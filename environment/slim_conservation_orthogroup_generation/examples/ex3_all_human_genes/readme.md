# construct groups for all genes in an organism
This is the way that I've used the pipeline mostly. I generated groups of filtered, clustered, least divergent orthologs (LDOs) for all of the human genes and then accessed the output files for further conservation analysis later. In this way I used this output as a database of alignments to annotate tables with conservation scores.

The `run.sh` script runs two steps:
1. `pipeline_all_genes_in_species.py` to construct the groups
2. `create_filemap.py` to create a filemap that maps the odb_gene_id to the groups generated for that gene id.


## construct groups: `pipeline_all_genes_in_species.py`

to see the command line parameters:
```bash
$ python "../../orthodb_tools/scripts/pipeline_all_genes_in_species.py" --help
```
```
usage: pipeline_all_genes_in_species.py [-h] [-c <file>] [-n <int>] [-s <str>] [-l [<list> ...]] [-o]

run the pipeline for all genes in an organism

options:
  -h, --help            show this help message and exit
  -c <file>, --config <file>
                        path to config file (default: None)
  -n <int>, --n_cores <int>
                        number of cores to use (default: 6)
  -s <str>, --species_id <str>
                        species id to use (default: 9606_0)
  -l [<list> ...], --og_levels [<list> ...]
                        list of phylogenetic levels at which to construct ortholog groups (default: ['Eukaryota', 'Mammalia', 'Metazoa', 'Tetrapoda', 'Vertebrata'])
  -o, --overwrite       if flag is provided and the main_output_folder exists, it will be removed and overwritten by the new files. Otherwise, an error will be raised if the folder exists (default: False)
```
The default number of cores is the number of cores on your machine -2.

if you are running on a workstation, I would run it using nohup and save the output to a file like this:
`nohup python "../../orthodb_tools/scripts/pipeline_all_genes_in_species.py" -c "./params.yml" --species_id "9606_0" > pipeline_genes_ids_in_species.out &`

If using a cluster, I would specifiy the number of cores to use in the script with the `--n_cores` parameter

If you wanted to do this for another organism, you would just need to change the variable `--species_id` to the appropriate species code (e.g. `9606_0` for human).

You can find the species code for any organism in orthoDB in the `odb11v0_species.tab` file that you downloaded from the orthoDB website. The file in this repository is actually the full sized file and so you could just use that one as well (`../../data/orthoDB_sample_data/odb11v0_species.tab`).


Or you could do something like this:
```python
import orthodb_tools.env_variables.env_variables as env
ODB_DATABASE = env.orthoDB_database()

for k,v in ODB_DATABASE.data_species_dict.items():
    if 'sapiens' in v:
        print(k,v)
```
## construct database key: `create_filemap.py`
The script `../../orthodb_tools/scripts/create_filemap.py` creates a json file that maps the odb_gene_id to the files that were generated for that gene id. This is useful for accessing the files later. <br>

```bash
$ python "../../orthodb_tools/scripts/create_filemap.py" --help
```
```bash
usage: create_filemap.py [-h] --main_output_folder MAIN_OUTPUT_FOLDER --output_file OUTPUT_FILE

run this after the pipeline to create a single json file that maps the gene ids to the files output by the pipeline
creates a json file mapping odb gene ids to the files output by the pipeline
will create a file looking like:
{
    "odb_gene_id": {
        "Eukaryota": {
            "alignment_file": "path/to/file",
            "info_file": "path/to/file",
        },
        "Vertebrata": {
            "alignment_file": "path/to/file",
            "info_file": "path/to/file",
        },
        ...
    },
    ...
}

options:
  -h, --help            show this help message and exit
  --main_output_folder MAIN_OUTPUT_FOLDER
                        The main pipeline output folder. The folder should contain the folder `info_jsons/`, where it will look for the json files.
  --output_file OUTPUT_FILE
                        path of the output json file
```


## some more details

If you look at the script (`../../orthodb_tools/scripts/pipeline_all_genes_in_species.py`), this line retrieves all of the odb_gene_id's for an organism code:
`odbgeneid_list = sql_queries.get_all_odbids_from_species_id(species_id)`
The `odbgeneid_list` is then used to run the pipeline for each gene id in the list.<br>


