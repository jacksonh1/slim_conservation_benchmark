# examples of how to use the pipeline

## examples shown here:
In the examples shown here, the scripts are executed via a bash script called `run.sh`. Each folder (except `ex2_single_gene`) contains a `run.sh` script that can be run via `bash run.sh` to recreate the outputs in each folder.


### `./uniprotid_mapping/`
A simple script to map uniprot ids to odb_gene_ids. Takes a table with a column of uniprot ids and tries to map each uniprot id to a odb_gene_id. outputs a copy of the table with an added column of odb_gene_ids<br>
run `python ../orthodb_tools/scripts/map_uniprotid.py --help` to see more.

### `./ex1_single_gene/`:

In this example, the pipeline is run for a single gene using the script as a command line script (via a one-line bash script `./ex1_single_gene/run.sh`).

### `./ex2_table_with_uniprot_ids/`:

In this example, a table with a column of uniprot ids is imported and ortholog groups are constructed for each uniprot id. The main pipeline is imported and used in a simple script `./ex2_table_with_uniprot_ids/get_orthologs_from_table.py`. The resulting file paths are mapped back to the original table and the table is written to a new file. <br>

### `./ex3_all_human_genes/`:
In this example, there is a script that runs the pipeline for all of the human genes in the database. The aim is to construct a "database" of precomputed orthogroups.<br>
- pipeline is imported and used in the script `../orthodb_tools/scripts/pipeline_all_genes_in_species.py`
- uses multiprocessing
- A filemap is also created after the pipeline is run using the script `../orthodb_tools/scripts/create_filemap.py`

note: Running on the full dataset using 62 cores, this pipeline took a few hours to run and generated around 10 Gb of data. When I ran it with align=True, it generated around 100 Gbs of data and takes considerably longer (~1 day). <br>

### `./ex4_create_database_from_table/`:
In this example, I use the script `../orthodb_tools/scripts/pipeline_input_table.py` to create a database of orthogroups from a table. This is similar to the previous example, but the genes are retrieved from a table instead of being run for all genes in a species. <br>

### `./pipeline_walkthrough_for_single_gene/`:

This folder contains a notebook that walks through the pipeline step by step to show how the source code works.

# workflow example

1. database construction: construct groups for all human genes (using the script in `./ex3_all_human_genes/`)
    - generates a "database" of:
        - orthogroup info files `info_jsons/`
        - alignments: `alignments/`
        - file map: `filemap.json`
2. for a table of candidate motifs from human genes with a column of uniprot ids, map the uniprot ids to odb_gene_ids (using the script `../orthodb_tools/scripts/map_uniprotid.py`).
3. use the database filemap (e.g. `filemap.json`) to access the database files from the odb_gene_ids

