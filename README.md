Written by Jackson Halpin <br>

This work was supported by the National Institutes of Health under Award Number R35GM149227. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

data download and processing scripts are in `./src/data_processing/`<br>

Building the benchmark requires my [pipeline 1](https://github.com/jacksonh1/orthogroup_generation) and IUPRED2A. You can download IUPRED2A from [here](https://iupred2a.elte.hu/download_new).<br>
You then have to add the path to the downloaded iupred folder to the file `./src/data_processing/benchmark_generation/build_table/iupred_tools.py` as the variable `IUPRED_DIR` variable.<br>


The scripts were run in the following order:
1. `./src/data_processing/all_human_odb_sequences_removed_duplicates/get_all_human_proteins_in_odb.ipynb`
2. `./src/data_processing/ELM_instances/build_instance_map.ipynb`
3. `./src/data_processing/benchmark_generation/build_table/benchmark_table.ipynb`
4. `./src/data_processing/benchmark_generation/generate_orthogroup_database/run.sh`


They are mostly included as jupyter notebooks for documentation reasons
