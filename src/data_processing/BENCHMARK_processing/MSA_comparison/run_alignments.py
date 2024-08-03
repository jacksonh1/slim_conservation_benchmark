# %%
import json
import multiprocessing
from pathlib import Path

import pandas as pd
from Bio import SeqIO, Seq
import copy

import orthodb_tools.tools.cli_wrappers as cli
from orthodb_tools.tools import general_utils as tools
from Bio.SeqRecord import SeqRecord

LEVELS = ["Vertebrata", "Metazoa", "Tetrapoda"]
N_ALIGN_THREADS = 4
MUSCLE_OUTPUT_FOLDER = Path(
    "../../../../benchmark/benchmark_v4/p2_orthogroups/orthogroups/alignments_muscle"
)
MUSCLE_OUTPUT_FOLDER.mkdir(exist_ok=True, parents=True)
CLUSTALO_OUTPUT_FOLDER = Path(
    "../../../../benchmark/benchmark_v4/p2_orthogroups/orthogroups/alignments_clustalo"
)
CLUSTALO_OUTPUT_FOLDER.mkdir(exist_ok=True, parents=True)
DATABASE_KEY_FILE = (
    "../../../../benchmark/benchmark_v4/p2_orthogroups/orthogroups/database_key.json"
)
N_PROC = 15
# get true positives from benchmark
TABLE_FILE = (
    "../../../../benchmark/benchmark_v4/p3_conservation/benchmark_table_ANNOTATED.csv"
)

# %%


def strip_dashes_from_str(seq_str):
    """remove `-` characters from a sequence string

    Parameters
    ----------
    seq_str : str
        sequence string

    Returns
    -------
    str
        sequence string with `-` characters removed
    """
    return seq_str.replace("-", "")


def strip_dashes_from_sequences(sequences: list[SeqRecord]) -> list[SeqRecord]:
    """remove `-` characters from sequences in a list of SeqRecord objects

    Parameters
    ----------
    sequences : list
        list of SeqRecord objects

    Returns
    -------
    list
        list of SeqRecord objects with `-` characters removed from sequences
    """
    sequences = copy.deepcopy(sequences)
    for s in sequences:
        s.seq = Seq.Seq(strip_dashes_from_str(str(s.seq)))
    return sequences


def get_remaining_gene_ids(
    gene_ids, level, database_dict, muscle_output_folder, clustalo_output_folder
):
    remaining_gene_ids = []
    for gene_id in gene_ids:
        d = database_dict[gene_id]
        if level not in d:
            continue
        aln_file = Path(d[level]["alignment_file"])
        clustal_output_file = clustalo_output_folder / f"{aln_file.name}"
        muscle_output_file = muscle_output_folder / f"{aln_file.name}"
        if clustal_output_file.exists() and muscle_output_file.exists():
            continue
        remaining_gene_ids.append(gene_id)
    return remaining_gene_ids


# %%
def run_alignment_from_keydict(
    database_key_level_dict: dict,
    level: str,
    n_align_threads: int,
    muscle_output_folder: Path,
    clustalo_output_folder: Path,
):
    if level not in database_key_level_dict:
        return
    aln_file = Path(database_key_level_dict[level]["alignment_file"])
    faimporter = tools.FastaImporter(aln_file)
    aln_seqrecords = faimporter.import_as_list()
    stripped_seqrecords = strip_dashes_from_sequences(aln_seqrecords)

    clustal_output_file = clustalo_output_folder / f"{aln_file.name}"
    if clustal_output_file.exists():
        print(f"Clustal alignment file {clustal_output_file} already exists")
    else:
        clustal_aln_seqrecords = cli.clustal_align_wrapper(
            stripped_seqrecords,
            alignment_type="basic",
            output_type="list",
            n_align_threads=n_align_threads,
        )
        with open(clustal_output_file, "w") as f:
            SeqIO.write(clustal_aln_seqrecords, f, "fasta")
    muscle_output_file = muscle_output_folder / f"{aln_file.name}"
    if muscle_output_file.exists():
        print(f"Muscle alignment file {muscle_output_file} already exists")
    else:
        muscle_aln_seqrecords = cli.muscle_align_wrapper(
            stripped_seqrecords,
            muscle_binary="muscle",
            output_type="list",
            n_align_threads=n_align_threads,
        )
        with open(muscle_output_file, "w") as f:
            SeqIO.write(muscle_aln_seqrecords, f, "fasta")


def main(
    database_key_file,
    level,
    n_align_threads,
    muscle_output_folder,
    clustalo_output_folder,
    table_file,
    n_proc=10,
):
    df = pd.read_csv(table_file)
    df = df[df["critical_error"].isnull()]
    # filter out this really long gene (will be removed in the next benchmark)
    # it's ank3 which is not really a protein that I'm interested in anyway.
    # It has a lot of isoforms and the corresponding isoform in each organism
    # is probably all really underdetermined
    df = df[df["gene_id"] != "9606_0:0027f1"]
    gene_ids = list(df[df["verified interaction"]]["gene_id"].unique())
    with open(database_key_file, "r") as f:
        database_key = json.load(f)
    gene_ids = get_remaining_gene_ids(
        gene_ids, level, database_key, muscle_output_folder, clustalo_output_folder
    )
    p = multiprocessing.Pool(n_proc)
    f_args = [
        (
            database_key[i],
            level,
            n_align_threads,
            muscle_output_folder,
            clustalo_output_folder,
        )
        for i in gene_ids
    ]
    p.starmap(run_alignment_from_keydict, f_args)
    p.close()
    p.join()


if __name__ == "__main__":
    for level in LEVELS:
        main(
            database_key_file=DATABASE_KEY_FILE,
            level=level,
            n_align_threads=N_ALIGN_THREADS,
            muscle_output_folder=MUSCLE_OUTPUT_FOLDER,
            clustalo_output_folder=CLUSTALO_OUTPUT_FOLDER,
            n_proc=N_PROC,
            table_file=TABLE_FILE,
        )

# %%
