# %%
import torch
import pairk
from slim_conservation_scoring.seqtools import esm_tools
from slim_conservation_scoring.seqtools import general_utils as tools
from pathlib import Path
import orthodb_tools.env_variables.env_variables as env
import copy

odbdatabase = env.orthoDBDatabase()

# %%

aln_dir = Path(
    "../../../../benchmark/benchmark_v4/p2_orthogroups/orthogroups/alignments/"
)
emb_output_dir = Path("/mnt/shared/jch/embeddings_tmp/")
fasta_files = list(aln_dir.glob("*.fasta"))

def fasta_to_ids(fasta_path):
    faimporter = tools.FastaImporter(fasta_path)
    return set(faimporter.import_as_dict().keys())

all_benchmark_seqs = set()
for f in fasta_files:
    all_benchmark_seqs.update(fasta_to_ids(f))
n = len(all_benchmark_seqs)
print(n)
counter = 0
mod = esm_tools.ESM_Model()
for i in all_benchmark_seqs:
    output_file = emb_output_dir / f"{i}.pt"
    if output_file.exists():
        counter += 1
        continue
    s = copy.deepcopy(str(odbdatabase.data_all_seqrecords_dict[i].seq))
    t = mod.encode(s, device="cuda")
    torch.save(t.to("cpu"), output_file)
    counter += 1
    if counter % 100 == 0:
        print(f"done {counter} of {n}")

# embdict = {}
# for i in ex1.full_length_dict:
#     embdict[i] = torch.load(f"./{i}.pt")
