`./p1_build_table/`:<br>
generate a table of known instances and background regex matches in the proteome for a selection of SLiMs

`./p2_generate_orthogroup_database/`:<br>
run pipeline 1 to generate a database of orthogroups for just the benchmark table.

`./p3_run_conservation_pipeline/`:<br>
run pipeline 2 to generate json files with motif information and add the filepaths to the benchmark table. also calculates fragment_pairwise_gapless matrices

`./p4_create_score_comparison_table/`:<br>
calculate a table to compare scores

`./p5_kibby_conservation/`:<br>
run the kibby conservation scoring method on the benchmark sequences. link - [Kibby]() 

`./MSA_comparison/`:<br>
run multiple sequence alignments on the true positive set of the benchmark

