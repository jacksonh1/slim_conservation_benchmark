# %%
import json
import time
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO
from Bio.SeqRecord import SeqRecord

from local_conservation_scores.tools import capra_singh_2007_scores, general
from local_env_variables import matrices
from local_seqtools import alignment_tools as aln_tools
from local_seqtools import general_utils as tools
import local_seqtools.esm_tools as esm_tools
from collections import defaultdict
import torch
import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
import concurrent.futures

# %%


def make_empty_kmer_ortho_df(positions, ortholog_ids: list[str]):
    cols = ["reference_kmer"] + ortholog_ids
    df = pd.DataFrame(
        index=positions,
        columns=cols,
    )
    return df


def get_kmers(seq: str, k: int):
    """Generates kmers from a sequence"""
    k2 = k - 1
    kmers = []
    for i in range(len(seq) - k2):
        kmers.append(seq[i : i + k])
    return kmers


def get_subsequences(indices, input_string, k):
    """Get subsequences of ortho string based on best indices"""
    char_list = list(input_string)
    start_indices = indices.view(-1, 1)
    end_indices = start_indices + k
    max_length = len(char_list)
    subsequences = [
        "".join(char_list[start:end]) for start, end in zip(start_indices, end_indices)
    ]
    return subsequences


def run_pairwise_kmer_emb_aln(
    reference_id: str,
    embedding_dict: dict,
    k: int,
):
    # get the reference sequence and remove it from the embedding_dict
    ref_seq_str, ref_seq_embedding = embedding_dict.pop(reference_id)
    kmers = get_kmers(ref_seq_str, k)
    positions = list(range(len(ref_seq_str) - (k - 1)))
    score_df = make_empty_kmer_ortho_df(positions, list(embedding_dict.keys()))
    subseq_df = make_empty_kmer_ortho_df(positions, list(embedding_dict.keys()))
    pos_df = make_empty_kmer_ortho_df(positions, list(embedding_dict.keys()))
    score_df.loc[positions, "reference_kmer"] = kmers
    subseq_df.loc[positions, "reference_kmer"] = kmers
    pos_df.loc[positions, "reference_kmer"] = kmers

    # Make expanded ref tensor
    expand_inds_ref = torch.arange(k).view(1, -1) + torch.arange(
        ref_seq_embedding.shape[0]
    ).view(-1, 1)
    expand_inds_ref[ref_seq_embedding.shape[0] - (k - 1) :] = 0
    expand_inds_ref = (
        expand_inds_ref.unsqueeze(-1)
        .expand(-1, -1, ref_seq_embedding.shape[1])
        .to(dtype=torch.int64)
    )
    expand_ref = ref_seq_embedding.unsqueeze(1).expand(-1, expand_inds_ref.shape[1], -1)
    expand_ref = torch.gather(expand_ref, 0, expand_inds_ref)
    expand_ref = expand_ref[: ref_seq_embedding.shape[0] - (k - 1)].reshape(
        -1, k * expand_ref.shape[2]
    )

    # for each ortholog sequence
    for ortholog_id, v in embedding_dict.items():
        ortholog_seq = v[0]
        ortholog_embedding = v[1]
        if ortholog_seq is None or ortholog_embedding is None:
            continue
        if len(ortholog_seq) < k:
            continue

        # Make expanded ortholog tensor
        expand_inds_ortho = torch.arange(k).view(1, -1) + torch.arange(
            ortholog_embedding.shape[0]
        ).view(-1, 1)
        expand_inds_ortho[ortholog_embedding.shape[0] - (k - 1) :] = 0
        expand_inds_ortho = (
            expand_inds_ortho.unsqueeze(-1)
            .expand(-1, -1, ortholog_embedding.shape[1])
            .to(dtype=torch.int64)
        )
        expand_ortho = ortholog_embedding.unsqueeze(1).expand(
            -1, expand_inds_ortho.shape[1], -1
        )
        expand_ortho = torch.gather(expand_ortho, 0, expand_inds_ortho)
        expand_ortho = expand_ortho[: ortholog_embedding.shape[0] - (k - 1)].reshape(
            -1, k * expand_ortho.shape[2]
        )

        # Calculate pairwise distances and get stats
        pairwise_dists = torch.cdist(
            expand_ref, expand_ortho, p=2
        )  # Optional: compute_mode='donot_use_mm_for_euclid_dist'
        min_dists, min_dists_pos = torch.min(pairwise_dists, dim=-1)
        score_df.loc[positions, ortholog_id] = min_dists.cpu().numpy()
        pos_df.loc[positions, ortholog_id] = min_dists_pos.cpu().numpy()
        subseq_df.loc[positions, ortholog_id] = get_subsequences(
            min_dists_pos, ortholog_seq, k
        )
    return score_df, subseq_df, pos_df


def import_embedding_dict(
    embedding_dict_file: str | Path,
):
    embedding_dict = torch.load(embedding_dict_file)
    return embedding_dict


def run_embedding_alignment(
    reference_id: str,
    k: int,
    embedding_file: str | Path,
    **kwargs,
):
    embedding_dict = import_embedding_dict(embedding_file)
    score_df, subseq_df, pos_df = run_pairwise_kmer_emb_aln(
        reference_id,
        embedding_dict,
        k,
    )
    output_dict = {
        "score_dataframe": score_df.to_dict(orient="split"),
        "subseq_dataframe": subseq_df.to_dict(orient="split"),
        "position_dataframe": pos_df.to_dict(orient="split"),
        # "reciprocal_best_match_dataframe": rbm_df.to_dict(orient="split"),
    }
    return output_dict


def run_embedding_alignment_save_file(
    reference_id: str,
    k: int,
    embedding_file: str | Path,
    output_file: str | Path,
    **kwargs,
):
    output_dict = run_embedding_alignment(reference_id, k, embedding_file, **kwargs)
    with open(output_file, "w") as json_file:
        json.dump(output_dict, json_file)


# %%
# ==============================================================================
# // TITLE
# ==============================================================================


def process_row(row, level, score_key, output_folder, lflank, rflank):
    og = group_tools.ConserGene(row["json_file"])
    if level not in og.levels_passing_filters:
        return row
    lvlo = og.get_level_obj(level)
    embedding_file = row["embedding_file"]
    if embedding_file is None:
        return row
    hit_st_idr = og.hit_start_position - og.idr_start
    hit_end_idr = og.hit_end_position - og.idr_start
    query_idr = og.query_sequence[og.idr_start : og.idr_end + 1]
    flanked_hit_st_idr, flanked_hit_end_idr, flanked_hit = tools.pad_hit(
        query_idr, hit_st_idr, hit_end_idr, lflank, rflank
    )
    orig_hit_st_in_flanked_hit = hit_st_idr - flanked_hit_st_idr
    orig_hit_end_in_flanked_hit = hit_end_idr - flanked_hit_st_idr
    k = len(flanked_hit)
    # if score_key in og.info_dict["orthogroups"][level]["conservation_scores"]:
    #     print(f"WARNING: {score_key} already exists in {og.reference_index}-{level}")
    #     print(f"skipping {score_key} for {og.reference_index}-{level}")
    #     return row
    output_folder = Path(output_folder)
    output_file = output_folder / f"{og.reference_index}-{score_key}.json"
    output_file = output_file.resolve()
    run_embedding_alignment_save_file(
        reference_id=og.query_gene_id,
        k=k,
        embedding_file=embedding_file,
        output_file=output_file,
    )
    og.info_dict["orthogroups"][level]["conservation_scores"][f"{score_key}"] = {}
    score_dict = og.info_dict["orthogroups"][level]["conservation_scores"][
        f"{score_key}"
    ]
    score_dict["file"] = str(output_file)
    score_dict["flanked_hit"] = flanked_hit
    score_dict["flanked_hit_start_position_in_idr"] = flanked_hit_st_idr
    score_dict["original_hit_st_in_flanked_hit"] = orig_hit_st_in_flanked_hit
    score_dict["original_hit_end_in_flanked_hit"] = orig_hit_end_in_flanked_hit
    # score_dict["score_params"] = {k: v for k, v in score_params.items() if k != "mod"}
    score_dict["lflank"] = lflank
    score_dict["rflank"] = rflank
    og._overwrite_json()
    return row


def process_chunk(chunk, level, score_key, output_folder, lflank, rflank):
    return chunk.apply(
        process_row,
        level=level,
        score_key=score_key,
        output_folder=output_folder,
        lflank=lflank,
        rflank=rflank,
        axis=1,
    )


def main(
    input_table: str | Path,
    input_embedding_mapping: str | Path,
    output_folder: str | Path,
    level: str,
    score_key: str,
    lflank: int,
    rflank: int,
    # overwrite: bool = False,
    n_cores: int = 60,
):
    output_folder = Path(output_folder)
    output_folder.mkdir(exist_ok=True)
    df = pd.read_csv(input_table)
    embdf = pd.read_csv(
        input_embedding_mapping,
        header=None,
        names=["reference_index", "embedding_file"],
    )
    df = df.merge(embdf, on="reference_index", how="left")
    df_filtered = df[df["critical_error"].isna()].copy()
    # df_filtered = df_filtered[df_filtered['reference_index']==2]
    # for i, row in df_filtered.iterrows():
    #     process_row(row, level, score_key, output_folder, lflank, rflank)
    # return
    chunks = np.array_split(df_filtered, n_cores)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        processed_chunks = list(
            executor.map(
                process_chunk,
                chunks,
                [level] * n_cores,
                [score_key] * n_cores,
                [output_folder] * n_cores,
                [lflank] * n_cores,
                [rflank] * n_cores,
            )
        )
    df_result = pd.concat(processed_chunks)
    return df_result


N_CORES = 60
TABLE_FILE = (
    "../p3_conservation_analysis_pipeline/benchmark_table_renamed_ANNOTATED.csv"
)
EMBEDDING_MAPPING_FILE = "./emb_mappings.csv"
SCORE_KEY = "fragpair_gapless_embedding_lf5_rf5"
OUTPUT_FOLDER = "./embedding_scores"
LFLANK = 5
RFLANK = 5
LEVEL = "Vertebrata"
# OVERWRITE = True

if __name__ == "__main__":
    main(
        input_table=TABLE_FILE,
        input_embedding_mapping=EMBEDDING_MAPPING_FILE,
        output_folder=OUTPUT_FOLDER,
        level=LEVEL,
        score_key=SCORE_KEY,
        lflank=LFLANK,
        rflank=RFLANK,
        n_cores=N_CORES,
    )


# df_result = pd.concat(processed_chunks)
