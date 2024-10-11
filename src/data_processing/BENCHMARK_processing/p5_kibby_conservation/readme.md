
The scripts/notebooks (in order of execution):
- `./build_fasta_for_kibby.ipynb`
    - environment: slim_conservation
    - build fasta file for all of the query protein sequences in the benchmark
- `./run_kibby.sh`
    - environment: kibby
    - run the kibby conservation scoring method on the benchmark sequences (uses the `conservation_from_fasta.py` script from the kibby repo)
