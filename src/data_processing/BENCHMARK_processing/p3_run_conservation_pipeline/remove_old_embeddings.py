# %%
import pairk
from slim_conservation_scoring.seqtools import general_utils as tools
from pathlib import Path
import orthodb_tools.env_variables.env_variables as env
import copy
import shutil


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
# %%


# f2make = set(['a','b','c'])
# fhave = set(['a','b', 'd'])
# fhave.difference(f2make)

embedding_files = set([i.stem for i in list(emb_output_dir.glob("*.pt"))])
trash = emb_output_dir / "trash"
trash.mkdir(exist_ok=True)
# remove embeddings that are not in the benchmark
for i in embedding_files.difference(all_benchmark_seqs):
    x = emb_output_dir / f"{i}.pt"
    # shutil.rmtree(x)
    shutil.move(emb_output_dir / f"{i}.pt", trash / f"{i}.pt")

shutil.rmtree(trash)