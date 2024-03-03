import copy
import os
from pathlib import Path

import dotenv
import pandas as pd
# from attrs import define, field
from attrs import frozen
from Bio import SeqIO

dotenv_path = os.path.join(os.path.dirname(__file__), '.env')
dotenv.load_dotenv(dotenv_path)
src_dir = Path(dotenv_path).parent.parent

# this should be absolute path
IUPRED_DIR = os.environ['IUPRED2A_LIB_DIR']

phylogeny_lvl_ordering_file = src_dir / os.environ['PHYLOGENY_LVL_ORDERING']
with open(phylogeny_lvl_ordering_file, 'r') as f:
    PHYLOGENY_LVL_ORDERING = [l.strip() for l in f.readlines()]

colormap_dir = src_dir / os.environ['COLORMAP_DIR']
COLOR_MAP_FILES = {
    'clustal': colormap_dir / 'clustal_hex_map.json',
}

orthodb_dir = Path(os.environ['ORTHODB_DATA_DIR'])
MAFFT_EXECUTABLE = os.environ['MAFFT_EXECUTABLE']
MAFFT_ADDITIONAL_ARGUMENTS = os.environ['MAFFT_ADDITIONAL_ARGUMENTS']
CD_HIT_EXECUTABLE = os.environ['CD_HIT_EXECUTABLE']
CD_HIT_ADDITIONAL_ARGUMENTS = os.environ['CD_HIT_ADDITIONAL_ARGUMENTS']

# ==============================================================================
# // getting odb filepaths
# ==============================================================================

@frozen
class orthoDB_files_object:
    all_seqs_fasta: str = str(orthodb_dir / "odb11v0_all_og_fasta.tab")
    all_seqs_sqlite: str = str(orthodb_dir / "odb11v0_all_og.sqlite")
    gene_refs_tsv: str = str(orthodb_dir / "odb11v0_genes.tab")
    gene_refs_sqlite: str = str(orthodb_dir / "odb11v0_genes.sqlite")
    gene_xrefs_tsv = str(orthodb_dir / "odb11v0_gene_xrefs.tab")
    gene_xrefs_sqlite = str(orthodb_dir / "odb11v0_gene_xrefs.sqlite")
    ogs_tsv: str = str(orthodb_dir / "odb11v0_OGs.tab")
    ogs_sqlite: str = str(orthodb_dir / "odb11v0_OGs.sqlite")
    OG2genes_tsv: str = str(orthodb_dir / "odb11v0_OG2genes.tab")
    OG2genes_sqlite: str = str(orthodb_dir / "odb11v0_OG2genes.sqlite")
    levels_tsv: str = str(orthodb_dir / "odb11v0_levels.tab")
    levels2species_tsv: str = str(orthodb_dir / "odb11v0_level2species.tab")
    species_tsv: str = str(orthodb_dir / "odb11v0_species.tab")

orthoDB_files = orthoDB_files_object()

# ==============================================================================
# // data loading functions
# ==============================================================================
def load_data_all_odb_seqs(database_files: orthoDB_files_object = orthoDB_files):
    data_all_seqrecords_dict = SeqIO.index_db(
        str(database_files.all_seqs_sqlite),
        str(database_files.all_seqs_fasta),
        "fasta",
    )
    return data_all_seqrecords_dict


def load_data_species_df(database_files: orthoDB_files_object = orthoDB_files):
    species_df = pd.read_csv(
        database_files.species_tsv,
        sep="\t",
        header=None,
        names=[
            "NCBI id",
            "species ID",
            "species name",
            "assembly ID",
            "n clustered genes",
            "n OGs",
            "mapping type",
        ],
    )
    return species_df

def load_data_levels_df(database_files: orthoDB_files_object = orthoDB_files):
    levels_df = pd.read_csv(
        database_files.levels_tsv,
        sep="\t",
        header=None,
        names=[
            "level NCBI tax id",
            "level name",
            "total non-redundant count of genes in all underneath clustered species",
            "total count of OGs built on it",
            "total non-redundant count of species underneath",
        ],
    )
    return levels_df


class orthoDB_database:
    '''
    main class that holds the orthoDB data
    '''    

    def __init__(self, database_files: orthoDB_files_object = orthoDB_files):
        self.datafiles = database_files
        self.data_all_seqrecords_dict = load_data_all_odb_seqs(self.datafiles)
        self.data_levels_df = load_data_levels_df(self.datafiles)
        self.data_species_df = load_data_species_df(self.datafiles)
        # special dictionaries that I want to have available for quick lookup
        self.data_species_dict = self._load_data_species_dict()
        self.data_levels_taxid_name_dict = self._load_data_levels_taxid_name_dict()

    def _load_data_species_dict(self):
        return (
            self.data_species_df[["species ID", "species name"]]
            .set_index("species ID")
            .to_dict()["species name"]
        )

    def _load_data_levels_taxid_name_dict(self):
        return (
            self.data_levels_df[["level NCBI tax id", "level name"]]
            .set_index("level NCBI tax id")
            .to_dict()["level name"]
        )

    def get_sequences_from_list_of_seq_ids(self, sequence_ids: list[str]) -> dict:
        og_seq_dict = {}
        for odb_gene_id in sequence_ids:
            og_seq_dict[odb_gene_id] = self.data_all_seqrecords_dict[odb_gene_id]
        return copy.deepcopy(og_seq_dict)


# make the orthoDB_database object available as a global environment variable
# ODB_DATABASE = orthoDB_database()

# Below is an attempt to deal with the fact that the orthoDB file names might change with different versions
# file_wildcards = {
#     "all_seqs_fasta": "*all_og.fasta",
#     "gene_refs_tsv": "*_genes.tab",
#     "gene_xrefs_tsv": "*gene_xrefs.tab",
#     "ogs_tsv": "*OGs.tab",
#     "OG2genes_tsv": "*OG2genes.tab",
#     "levels_tsv": "*levels.tab",
#     "levels2species_tsv": "*level2species.tab",
#     "species_tsv": "*_species.tab",
# }

# odb_data_files = {}
# for k, v in file_wildcards.items():
#     try:
#         odb_data_files[k] = next(orthodb_dir.glob(v))
#     except StopIteration:
#         raise FileNotFoundError(f"Could not find file for `{k}`. \nSearched with wildcard {orthodb_dir}/{v}")
#     if len(list(orthodb_dir.glob(v))) > 1:
#         raise FileNotFoundError(f"Found multiple files for `{k}`. \nSearched with wildcard {orthodb_dir}/{v}")

# # The sqlite files are created by the scripts in scripts-gen_SQLite_dbs and I 
# # want to import and use this library to retrieve those file names before they are created
# # so I am not checking for their existence here.
# odb_data_files["all_seqs_sqlite"] = odb_data_files["all_seqs_fasta"].with_suffix(".sqlite")
# odb_data_files["gene_refs_sqlite"] = odb_data_files["gene_refs_tsv"].with_suffix(".sqlite")
# odb_data_files["gene_xrefs_sqlite"] = odb_data_files["gene_xrefs_tsv"].with_suffix(".sqlite")
# odb_data_files["ogs_sqlite"] = odb_data_files["ogs_tsv"].with_suffix(".sqlite")
# odb_data_files["OG2genes_sqlite"] = odb_data_files["OG2genes_tsv"].with_suffix(".sqlite")
# # later - add a way to output the sqlite files to a different directory specified in the .env file
# odb_data_files = {k: str(v) for k, v in odb_data_files.items()}

# # I'm not really sure if this is an appropriate use of attrs, but I really like
# # having the autocomplete and type checking when working with these files throughout
# # the codebase
# @frozen
# class orthoDB_files_object:
#     all_seqs_fasta: str
#     all_seqs_sqlite: str
#     gene_refs_tsv: str
#     gene_refs_sqlite: str
#     gene_xrefs_tsv: str
#     gene_xrefs_sqlite: str
#     ogs_tsv: str
#     ogs_sqlite: str
#     OG2genes_tsv: str
#     OG2genes_sqlite: str
#     levels_tsv: str
#     levels2species_tsv: str
#     species_tsv: str

# orthoDB_files = orthoDB_files_object(**odb_data_files)
