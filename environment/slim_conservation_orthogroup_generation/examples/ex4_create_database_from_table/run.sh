python "../../orthodb_tools/scripts/map_uniprotid.py" -i "./table.csv" --uni_column "Uniprotid"
python "../../orthodb_tools/scripts/pipeline_input_table.py" -c "./params.yaml" --table "./table_mapped_odbgeneid.csv" --odb_gene_id_column gene_id
python "../../orthodb_tools/scripts/create_filemap.py" --main_output_folder "./output/" --output_file "./output/database_key.json"