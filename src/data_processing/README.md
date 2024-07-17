The directories here are for data downloading/processing and analysis of the benchmark.

I used different environments for this:
- slim_conservation_orthogroup_generation - [link](https://github.com/jacksonh1/slim_conservation_orthogroup_generation)
- slim_conservation_scoring - [link](https://github.com/jacksonh1/slim_conservation_scoring)
- kibby - [link](https://github.com/esbgkannan/kibby)

subdirectories (in order of execution):
- `./all_human_odb_sequences_removed_duplicates/`:
  - environment: slim_conservation_orthogroup_generation
  - get all human proteins in the orthodb, remove duplicates by clustering. This is the set of protein sequences that are used for the benchmark.
- `./ELM_instances/`:
  - environment: slim_conservation_orthogroup_generation
  - Downloads ELM instances from the ELM and verifies that the annotations match the given motifs. It then maps each ELM instance (position and sequence of motif) to the corresponding sequence in the orthodb.
- `./BENCHMARK_processing/`:
  - environment: All three
  - This is the main benchmark processing directory. It contains subdirectories for each step of the benchmark generation and scoring.
  - See the [README](./README.md) in this directory for more information.
