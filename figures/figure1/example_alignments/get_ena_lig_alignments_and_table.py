import shutil
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
from pathlib import Path
import local_seqtools.general_utils as tools

table_file = "../../../benchmark/benchmark_v3/MSA_comparison/benchmark_table_renamed_MSA_comparison_ANNOTATED.csv"
table = pd.read_csv(table_file)
table = table[
    [
        "reference_index",
        "Organism",
        "Primary_Acc",
        "Accessions",
        "UniprotID",
        "regex",
        "hit_sequence",
        "gene_id",
        "hit start position",
        "hit end position",
        "verified interaction",
        "ELM_motif_class",
        "name",
        "critical_error",
        "json_file",
    ]
]

enah_df = table[(table["ELM_motif_class"]=="enah_LPPPP_FPPPP")&(table["verified interaction"])].copy()
OUTPUT_DIR = Path("./ena_motif_alignments")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# %%
# ==============================================================================
# // copy over alignments and slices
# ==============================================================================

LEVELS = ["Metazoa", "Vertebrata"]
for i, row in enah_df.iterrows():
    folder = OUTPUT_DIR / f"{row['reference_index']}_{row['name']}_{row['gene_id']}"
    folder.mkdir(exist_ok=True, parents=True)
    json_path = Path(row["json_file"])
    og = group_tools.ConserGene(json_path)
    for level in LEVELS:
        if level not in og.levels_passing_filters:
            continue
        lvlo = og.get_level_obj(level)
        # copy the alignment file
        shutil.copy(lvlo.alignment_file, folder / lvlo.alignment_file.name)
        shutil.copy(lvlo.info_dict["aln_slice_file"], folder / Path(lvlo.info_dict["aln_slice_file"]).name)
    enah_df.loc[i, "folder"] = str(folder)

# %%
# ==============================================================================
# // add flank to the table to make it easier to identify hits you want to look at
# ==============================================================================

import local_env_variables.env_variables as env
data_all_seqrecords_dict = env.load_data_all_odb_seqs()

enah_df['odb_seq'] = enah_df['gene_id'].apply(lambda x: str(data_all_seqrecords_dict[x].seq) if x in data_all_seqrecords_dict.keys() else False)
enah_df['flanked_hit'] = enah_df.apply(lambda x: tools.pad_with_aas_or_gaps(x['odb_seq'], x['hit start position'], x['hit end position']+1, flank=5), axis=1, result_type='expand')




enah_df.to_csv("enah_table.csv")