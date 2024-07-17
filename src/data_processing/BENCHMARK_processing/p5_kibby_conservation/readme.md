
scripts/notebooks were run with these environments:
- slim_conservation_orthogroup_generation - [link](https://github.com/jacksonh1/slim_conservation_orthogroup_generation)
- kibby - [link](https://github.com/esbgkannan/kibby)

The scripts/notebooks (in order of execution):
- `./build_fasta_for_kibby.ipynb`
    - environment: slim_conservation_orthogroup_generation
    - build fasta file for all of the query protein sequences in the benchmark
- `./run_kibby.sh`
    - environment: kibby
    - run the kibby conservation scoring method on the benchmark sequences (uses the `conservation_from_fasta.py` script from the kibby repo)
