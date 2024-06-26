benchmark_folder="../../../../benchmark/benchmark_v3/"
output_folder="$benchmark_folder""p2_orthogroups/"

mkdir $output_folder
python "../../../local_scripts/pipeline_input_table.py" -c "./params.yaml" --table "$benchmark_folder""p1_table/benchmark_table.csv" --odb_gene_id_column gene_id -l Vertebrata -n 60
python "../../../local_scripts/create_filemap.py" --main_output_folder "$output_folder""orthogroups" --output_file "$output_folder""orthogroups/database_key.json"
