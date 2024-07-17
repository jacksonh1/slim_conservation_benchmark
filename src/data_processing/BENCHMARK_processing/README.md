I used different environments to generate the benchmark and quantify conservation:
slim_conservation_orthogroup_generation - [link](https://github.com/jacksonh1/slim_conservation_orthogroup_generation)
slim_conservation_scoring - [link](https://github.com/jacksonh1/slim_conservation_scoring)
kibby - [link](https://github.com/esbgkannan/kibby)

The pipeline that was used is stated in readme files or the top of the notebooks

`./p1_build_table/`:<br>
- environment: slim_conservation_orthogroup_generation
- generate a table of known instances and background regex matches in the proteome for a selection of SLiMs

`./p2_generate_orthogroup_database/`:<br>
- environment: slim_conservation_orthogroup_generation
- run orthogroup generation pipeline to generate a database of orthogroups for just the benchmark table.

`./p3_run_conservation_pipeline/`:<br>
- environment: slim_conservation_scoring
- run conservation scoring pipeline to generate json files with motif information and add the filepaths to the benchmark table

`./p4_create_score_comparison_table/`:<br>
- environment: slim_conservation_scoring
- generate tables of scores for comparison

`./p5_kibby_conservation/`:<br>
- environment: kibby
- run the kibby conservation scoring method on the benchmark sequences

`./MSA_comparison/`:<br>
- environment: slim_conservation_orhtogroup_generation and slim_conservation_scoring
- generate multiple sequence alignments for the true positive subset of the benchmark
- see [readme.md](./MSA_comparison/readme.md) for commands

`./kibby_scores/`:<br>
- environment: kibby and slim_conservation_orhtogroup_generation
- generate combined fasta file for all of the query protein sequences in the benchmark. Calculate kibby scores for the sequences
- see [readme.md](./kibby_scores/readme.md) for more info



