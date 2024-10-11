I used different environments to generate the benchmark and quantify conservation:
- slim_conservation - [../../environment/](../../environment/)
- kibby - [link](https://github.com/esbgkannan/kibby)

The environment that was used is stated in readme files or the top of the notebooks

`./p1_build_table/`:<br>
- environment: slim_conservation
- generate a table of known instances and background regex matches in the proteome for a selection of SLiMs

`./p2_generate_orthogroup_database/`:<br>
- environment: slim_conservation
- run orthogroup generation pipeline to generate a database of orthogroups for just the benchmark table.

`./p3_run_conservation_pipeline/`:<br>
- environment: slim_conservation
- run conservation scoring pipeline to generate json files with motif information and add the filepaths to the benchmark table

`./p4_create_score_comparison_table/`:<br>
- environment: slim_conservation
- generate tables of scores for comparison

`./p5_kibby_conservation/`:<br>
- environment: slim_conservation and kibby
- generate fasta files for the benchmark query sequences
- run the kibby conservation scoring method on the benchmark sequences

`./MSA_comparison/`:<br>
- environment: slim_conservation
- generate multiple sequence alignments for the true positive subset of the benchmark
- see [readme.md](./MSA_comparison/readme.md) for commands

`./kibby_scores/`:<br>
- environment: slim_conservation
- generate combined fasta file for all of the query protein sequences in the benchmark. Calculate kibby scores for the sequences
- see [readme.md](./p5_kibby_conservation/readme.md) for more info



