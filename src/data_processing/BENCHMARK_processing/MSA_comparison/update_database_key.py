# %%
import json
from pathlib import Path
import pandas as pd
import re
import copy

LEVELS = ["Vertebrata", "Metazoa", "Tetrapoda"]
MUSCLE_ALN_FOLDER = Path(
    "../../../../benchmark/benchmark_v4/p2_orthogroups/orthogroups/alignments_muscle"
)
CLUSTALO_ALN_FOLDER = Path(
    "../../../../benchmark/benchmark_v4/p2_orthogroups/orthogroups/alignments_clustalo"
)
DATABASE_KEY_FILE = Path(
    "../../../../benchmark/benchmark_v4/p2_orthogroups/orthogroups/database_key.json"
)
TABLE_FILE = Path("../../../../benchmark/benchmark_v4/p3_conservation/benchmark_table_ANNOTATED.csv")

def get_file_dict(aln_folder, level):
    p = re.compile(f"(\d+_\d)_([\d\w]+)_{level}")
    filelist=list(Path(aln_folder).rglob("*.fasta"))
    filedict={}
    for f in filelist:
        if level not in f.name:
            continue
        m=re.findall(p, f.name)
        assert len(m) == 1, f"wrong number of matches: {m}, {f.name}, {p}"
        assert len(m[0]) == 2, f"wrong number of groups in match: {m}"
        ref_id=f'{m[0][0]}:{m[0][1]}'
        filedict[ref_id]=f
    print(f"Found {len(filedict)} files for {level} in {aln_folder}")
    return filedict

def add_files_2_newkey(
    muscle_aln_folder: Path,
    clustalo_aln_folder: Path,
    level,
    full_database_key: dict,
    msa_comparison_database_key: dict
):
    muscle_alnfiledict = get_file_dict(muscle_aln_folder, level)
    clustalo_alnfiledict = get_file_dict(clustalo_aln_folder, level)

    for ref_id, aln_file in muscle_alnfiledict.items():
        if ref_id not in full_database_key:
            continue
        if ref_id not in msa_comparison_database_key:
            msa_comparison_database_key[ref_id] = {}
        msa_comparison_database_key[ref_id][f'{level}'] = {}
        msa_comparison_database_key[ref_id][f'{level}'] = copy.deepcopy(full_database_key[ref_id][f'{level}'])
        msa_comparison_database_key[ref_id][f'{level}-Muscle'] = {}
        msa_comparison_database_key[ref_id][f'{level}-Muscle']['alignment_file'] = str(aln_file.resolve())

    for ref_id, aln_file in clustalo_alnfiledict.items():
        if ref_id not in msa_comparison_database_key:
            raise ValueError(f"ref_id not found in database key: {ref_id}, {aln_file}")
        msa_comparison_database_key[ref_id][f'{level}-clustalo'] = {}
        msa_comparison_database_key[ref_id][f'{level}-clustalo']['alignment_file'] = str(aln_file.resolve())

    return msa_comparison_database_key

    
with open(DATABASE_KEY_FILE, "r") as f:
    database_key = json.load(f)

msa_comparison_database_key = {}
for level in LEVELS:
    msa_comparison_database_key = add_files_2_newkey(
        MUSCLE_ALN_FOLDER,
        CLUSTALO_ALN_FOLDER, 
        level, 
        database_key, 
        msa_comparison_database_key
    )


with open(str(DATABASE_KEY_FILE).replace(".json", "-MSA_comparison.json"), "w") as f:
    json.dump(msa_comparison_database_key, f, indent=4)


df = pd.read_csv(TABLE_FILE)
df=df[df['critical_error'].isnull()]
# filter out this really long gene (will be removed in the next benchmark)
# it's ank3 which is not really a protein that I'm interested in anyway. 
# It has a lot of isoforms and the corresponding isoform in each organism
# is probably really underdetermined
df = df[df["gene_id"] != "9606_0:0027f1"]
df = df[df['verified interaction']]
df = df.drop(columns=['critical_error', 'json_file'])
df.to_csv(Path("../../../../benchmark/benchmark_v4/MSA_comparison/") / TABLE_FILE.name.replace("ANNOTATED", "MSA_comparison"), index=False)



