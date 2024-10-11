
commands (in order):
- activate slim_conservation conda environment
    - `mamba activate slim_conservation`
- run alignments:
    - `nohup python "./run_alignments.py" > run_alignments.out &`
- update database key:
    - `python "./update_database_key.py"`
- run conservation scoring:
    - `nohup bash ./score_conservation.sh > score_conservation.out &`
    - (contents of ./score_conservation.sh is just one line - `slim_conservation_scoring-pipeline -c "./params.yaml" -n 60`)
