import copy

import numpy as np
import pandas as pd
from alfpy.utils import distmatrix
from Bio import SeqIO

import orthodb_tools.sql_queries as sql_queries
import orthodb_tools.tools.alignment_tools as aln_tools
import orthodb_tools.tools.cli_wrappers as cli


def setup_df(seqrecord_dict_in: dict[str, SeqIO.SeqRecord]) -> pd.DataFrame:
    seqrecord_dict = copy.deepcopy(seqrecord_dict_in)
    seqrecord_list = [seq for seq in seqrecord_dict.values()]
    df = pd.DataFrame(columns=["id"], index=range(len(seqrecord_dict)))
    df["id"] = [seqrecord.id for seqrecord in seqrecord_list]
    df["organism"] = df["id"].apply(lambda x: sql_queries.odb_gene_id_2_species_id(x))
    df["sequence"] = df["id"].map(seqrecord_dict)
    return df


def addpid_by_msa(
    df_in: pd.DataFrame,
    query_seqrecord: SeqIO.SeqRecord,
    seqrecord_dict: dict[str, SeqIO.SeqRecord],
    n_align_threads: int = 8,
    **mafft_kwargs,
) -> pd.DataFrame:
    df = df_in.copy()
    seqrecord_list = [seq for seq in seqrecord_dict.values()]
    _, msa_seqrecord_dict = cli.mafft_align_wrapper(
        seqrecord_list, n_align_threads=n_align_threads, **mafft_kwargs
    )
    query_msa_seqrecord = msa_seqrecord_dict[query_seqrecord.id]  # type: ignore
    # query_msa_seqrecord = [i for i in msa_seqrecord_dict if i.id == query_seqrecord.id][ # type: ignore
    # 0
    # ]
    # could do this by applying a function to the `id` column but this seems a bit simpler
    pid_map_dict = {
        seq.id: aln_tools.compute_pairwise_percent_id_from_msa(query_msa_seqrecord, seq)
        for seq in msa_seqrecord_dict.values()
    }
    df["PID"] = df["id"].map(pid_map_dict)
    return df


def addpid_by_msa_by_organism(
    df_in: pd.DataFrame,
    query_seqrecord: SeqIO.SeqRecord,
    n_align_threads: int = 8,
    **mafft_kwargs,
) -> pd.DataFrame:
    df = df_in.copy()
    print("aligning sequences using mafft, one organism at a time")
    org_set = df["organism"].unique()
    pid_map_dict = {}
    # counter = 0
    for org_i in org_set:
        # counter += 1
        # get all the sequences for this organism
        seqs = list(df[df["organism"] == org_i]["sequence"].values)
        # add the query sequence to the list
        if query_seqrecord.id not in [seq.id for seq in seqs]:
            seqs.append(query_seqrecord)
        # align the sequences
        _, msa_i_dict = cli.mafft_align_wrapper(
            seqs, n_align_threads=n_align_threads, **mafft_kwargs
        )
        # compute the pairwise percent identity from the MSA
        for seq_id, seqrecord in msa_i_dict.items():
            pid_map_dict[seq_id] = aln_tools.compute_pairwise_percent_id_from_msa(
                seqrecord, msa_i_dict[query_seqrecord.id]
            )
        # if counter % 100 == 0:
        # print(f'finished {counter} of {len(org_set)} organisms')
    df["PID"] = df["id"].map(pid_map_dict)
    return df


def addpid_by_alfpy_google_distance(
    df_in: pd.DataFrame, query_seqrecord: SeqIO.SeqRecord
) -> pd.DataFrame:
    df = df_in.copy()
    print("comparing sequences using alignment free comparison (alfpy google distance)")
    org_set = df["organism"].unique()
    pid_map_dict = {}
    # counter = 0
    for org_i in org_set:
        # counter += 1
        # get all the sequences for this organism
        seqs = list(df[df["organism"] == org_i]["sequence"].values)
        # add the query sequence to the list
        if query_seqrecord.id not in [seq.id for seq in seqs]:
            seqs.append(query_seqrecord)
        matrix = aln_tools.alfpy_distance_matrix(seqs, word_size=2)
        id_list, query_row_similarity = aln_tools.query_alfpy_distance_matrix(
            query_seqrecord.id, matrix, similarity=True
        )
        for seq_id, similarity in zip(id_list, query_row_similarity):
            pid_map_dict[seq_id] = similarity
        # if counter % 100 == 0:
        # print(f'finished {counter} of {len(org_set)} organisms')
    df["PID"] = df["id"].map(pid_map_dict)
    return df


def addpid_by_pairwise(
    df_in: pd.DataFrame,
    query_seqrecord: SeqIO.SeqRecord,
    seqrecord_dict: dict[str, SeqIO.SeqRecord],
) -> pd.DataFrame:
    df = df_in.copy()
    seqrecord_list = [seq for seq in seqrecord_dict.values()]
    print("aligning sequences using pairwise alignment")
    pid_map_dict = {
        seq.id: aln_tools.align_and_get_PID(query_seqrecord, seq)
        for seq in seqrecord_list
    }
    df["PID"] = df["id"].map(pid_map_dict)
    return df


def get_LDOs_from_pids(df: pd.DataFrame, query_seqrecord: SeqIO.SeqRecord) -> list[str]:
    query_species_id = sql_queries.odb_gene_id_2_species_id(query_seqrecord.id)
    # remove sequences in the query organism that are not the query sequence
    df = df[(df["organism"] != query_species_id) | (df["id"] == query_seqrecord.id)]
    assert query_seqrecord.id in df["id"].values, "query sequence not found in df"
    # groupby select the closest sequence for each organism
    ldo_df = df.loc[df.groupby("organism")["PID"].idxmax()].copy()
    ldo_df = ldo_df.sort_values("PID", ascending=False)
    return list(ldo_df["id"].values)


def find_LDOs_main(
    seqrecord_dict: dict[str, SeqIO.SeqRecord],
    query_seqrecord: SeqIO.SeqRecord,
    pid_method: str = "alfpy_google_distance",
    n_align_threads: int = 8,
    **mafft_kwargs,
) -> tuple[pd.DataFrame, list[str]]:
    assert pid_method in [
        "msa_by_organism",
        "alfpy_google_distance",
        "pairwise",
        "msa",
    ], "LDO selection method not recognized. must be one of: msa_by_organism, alfpy_google_distance, pairwise, msa"
    df = setup_df(seqrecord_dict)
    if pid_method == "msa_by_organism":
        df = addpid_by_msa_by_organism(
            df, query_seqrecord, n_align_threads=n_align_threads, **mafft_kwargs
        )
    elif pid_method == "alfpy_google_distance":
        df = addpid_by_alfpy_google_distance(df, query_seqrecord)
    elif pid_method == "pairwise":
        df = addpid_by_pairwise(df, query_seqrecord, seqrecord_dict)
    elif pid_method == "msa":
        df = addpid_by_msa(
            df,
            query_seqrecord,
            seqrecord_dict,
            n_align_threads=n_align_threads,
            **mafft_kwargs,
        )
    return df, get_LDOs_from_pids(df, query_seqrecord)
