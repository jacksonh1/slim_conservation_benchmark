
Run with both [slim_conservation_orthogroup_generation](https://github.com/jacksonh1/slim_conservation_orthogroup_generation) and [slim_conservation_scoring](https://github.com/jacksonh1/slim_conservation_scoring).

commands (in order):
- activate slim_conservation_orthogroup_generation conda environment
    - `mamba activate slim_conservation_orthogroup_generation`
- run alignments:
    - `nohup python "./run_alignments.py" > run_alignments.out &`
- update database key:
    - `python "./update_database_key.py"`
- activate slim_conservation_scoring conda environment
    - `mamba activate slim_conservation_scoring`
- run conservation scoring:
    - `nohup python "../../../local_scripts/conservation_analysis.py" -c "./params.yaml" -n 60`
