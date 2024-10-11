# %%
import torch
import pairk
from transformers import AutoModel, AutoTokenizer
import torch
import numpy as np
from slim_conservation_scoring.seqtools import general_utils as tools
from pathlib import Path
import orthodb_tools.env_variables.env_variables as env
import copy


# %%

odbdatabase = env.orthoDBDatabase()
aln_dir = Path(
    "../../../../benchmark/benchmark_v4/p2_orthogroups/orthogroups/alignments/"
)
emb_output_dir = Path("/mnt/shared2/jch/dr_bert_embeddings_tmp/")
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


class DrBertEncoder:

    def __init__(
        self, checkpoint="../../../../data/DR-BERT/checkpoint-final/", threads: int = 1
    ):
        torch.set_num_threads(threads)
        self.tokenizer = AutoTokenizer.from_pretrained(
            checkpoint, clean_up_tokenization_spaces=True
        )
        self.model = AutoModel.from_pretrained(checkpoint)
        self.model = self.model.eval()
        # Get max input size from model config
        self.max_n_residues = self.model.config.max_position_embeddings - 4
        # self.max_n_residues = 20

    def embed(self, sequence, device="cuda"):
        inputs = self.tokenizer(
            sequence,
            return_tensors="pt",
        ).to(device)
        with torch.no_grad():
            outputs = self.model(**inputs)
            embeddings = outputs.last_hidden_state.to("cpu")[
                0
            ]  # Get embeddings and move to CPU
        return embeddings

    def encode(self, sequence, device="cuda", stride=50, padding=5):
        try:
            chunksize = self.max_n_residues - (padding * 2) - 2
            self.model = self.model.to(device)
            sequence_length = len(sequence)

            # If sequence is shorter than the model max length, encode it directly
            if sequence_length <= self.max_n_residues:
                return self.embed(sequence, device)

            # For longer sequences, split into chunks
            all_embeddings = []
            for i in range(
                0,
                sequence_length + 1,
                chunksize - (stride),
            ):
                # Define the window with overlap
                chunk = sequence[i : i + chunksize]
                if i == 0:
                    # pad the end with glycines
                    chunk = chunk + "G" * (padding)
                    chunk_embeddings = self.embed(chunk, device)
                    # remove end token and padding
                    all_embeddings.append(chunk_embeddings[: -(padding + 1)])
                elif i + chunksize >= sequence_length:
                    chunk = "G" * (padding) + chunk
                    chunk_embeddings = self.embed(chunk, device)
                    # remove the start token and padding
                    all_embeddings.append(chunk_embeddings[(padding + 1) :])
                    break
                else:
                    chunk = "G" * (padding) + chunk + "G" * (padding)
                    chunk_embeddings = self.embed(chunk, device)
                    # remove both start/end tokens and padding
                    all_embeddings.append(
                        chunk_embeddings[(padding + 1) : -(padding + 1)]
                    )
                # Combine embeddings, handling overlaps (e.g., average overlapping parts)
            full_embeddings = self._combine_chunks(all_embeddings, stride)
            assert (
                full_embeddings.shape[0] == sequence_length + 2
            ), f"sequence length mismatch: {full_embeddings.shape[0]} vs {sequence_length+2}"
            return full_embeddings

        except Exception as e:
            if device != "cpu":
                print(f"failed with {device}, trying cpu...")
                print(f"{e}")
                return self.encode(
                    sequence, device="cpu", stride=stride, padding=padding
                )
            else:
                print(f"failed on {sequence}")
                raise e

    def _combine_chunks(self, chunks, stride):
        for idx, chunk in enumerate(chunks):
            if idx == 0:
                combined_embeddings = chunk
            else:
                overlap_size = min(
                    combined_embeddings[-stride:].size(0), chunk[:stride].size(0)
                )
                overlap = torch.mean(
                    torch.stack(
                        [combined_embeddings[-overlap_size:], chunk[:overlap_size]]
                    ),
                    dim=0,
                )
                non_overlap = chunk[overlap_size:]
                combined_embeddings = torch.cat(
                    [combined_embeddings[:-overlap_size], overlap, non_overlap]
                )
        return combined_embeddings


# %%

counter = 0
mod = DrBertEncoder(threads=60)
for i in all_benchmark_seqs:
    output_file = emb_output_dir / f"{i}.pt"
    if output_file.exists():
        # print('file exists')
        counter += 1
        continue
    s = copy.deepcopy(str(odbdatabase.data_all_seqrecords_dict[i].seq))
    t = mod.encode(s, device="cuda", stride=50, padding=5)
    torch.save(t.to("cpu"), output_file)
    # print(f'processed {i}')
    counter += 1
    if counter % 100 == 0:
        print(f"done {counter} of {n}")
    # break


# %%
