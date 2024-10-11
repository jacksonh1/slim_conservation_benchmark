
Written by Jackson Halpin <br>

This work was supported by the National Institutes of Health under Award Number R35GM149227. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

# Table of contents
- [Table of contents](#table-of-contents)
- [orthologous groups for conservation analysis](#orthologous-groups-for-conservation-analysis)
  - [Pipeline overview:](#pipeline-overview)
  - [note on future improvements](#note-on-future-improvements)
- [necessary background knowledge (beginner level)](#necessary-background-knowledge-beginner-level)
- [setup TL;DR:](#setup-tldr)
- [Usage](#usage)
  - [overview](#overview)
    - [using the tools as a command line script](#using-the-tools-as-a-command-line-script)
    - [using the tools as a module](#using-the-tools-as-a-module)
  - [useful scripts: `./orthodb_tools/scripts/`](#useful-scripts-orthodb_toolsscripts)
  - [pipeline parameters](#pipeline-parameters)
    - [pipeline parameters explained](#pipeline-parameters-explained)
- [comments](#comments)
- [source code](#source-code)
- [tools used (and links):](#tools-used-and-links)

# orthologous groups for conservation analysis
This repository contains tools to retrieve and process ortholog groups from a local copy of the orthoDB database files ([link](https://www.orthodb.org/)). The pipeline finds a protein of interest in the database, retrieves its homologous proteins as defined by the orthoDB, and processes the group of homologs in preparation for downstream [conservation analysis](https://github.com/jacksonh1/motif_conservation_in_IDRs). <br>
- Note: I refer  to these sequences as orthologs in many places throughout this repo but you could argue that it's more accurate to refer to them as homologs depending on how you define the term. I use the term orthologs because that is what orthoDB calls them. see orthoDB [terminology](https://www.ezlab.org/orthodb_userguide.html#terminology)


## Pipeline overview:
1. **Find a query protein in the orthoDB database** (retrieve the corresponding orthoDB ID)
   - you can find a protein by looking up its uniprot ID in the orthoDB database. You can also set the orthoDB ID manually if you already know it. 
   - *Note: If it can't find an orthoDB ID for a given UniprotID, it doesn't mean that the protein is absent from the orthoDB. It could still be present but was retrieved from a different database and a uniprot ID was not mapped to it. View OrthoDB documentation for more info on where the sequences come from. I have not solved this problem. In the future, it would be nice to develop a way to search for the actual full length sequence using blast or something, if it fails to find the uniprot ID*
2. **Retrieve the orthoDB-defined groups of homologous sequences** (orthogroup IDs) containing the query protein
3. **Select a group based on phylogenetic level**
    - OrthoDB constructs groups of homologs at different phylogenetic levels (e.g. Eukaryota, Metazoa, Vertebrata, etc.). Thus, a protein will probably be part of multiple groups at different levels. Look at the orthoDB page for pcare as an example (https://www.orthodb.org/?query=pcare)
    - For short linear motifs, we have found that the Vertebrata level typically works well for conservation analysis
4. **Filter** out sequences in the group that are too short (relative to the length of the query sequence) or that contain non amino acid characters
5. **Filter to least divergent orthologs (LDOs)**:
   - For each organism in the group, select the sequence that is most similar to the query sequence such that there is only one sequence per organism
6. **Cluster the filtered LDOs using CD-HIT**
7. **Align the clustered sequences** (uses MAFFT by default)
8. **output** the alignment and the ortholog group information in a directory structure that is compatible with the conservation analysis pipeline (link)
    - the group information is output in the form of a json file that can be imported into python

## note on future improvements
- I wanted to avoid using a workflow manager (e.g. snakemake or nextflow) and a more sophisticated database structure because I wanted to keep the pipeline easy for anyone to use, however, I think the result is potentially more confusing. If I had the time, I might refactor everything to use a workflow manager to handle the pipeline and a different database structure to store the data. It might make the pipeline more efficient and easier to use. <br>

# necessary background knowledge (beginner level)
- **basic use of a unix terminal** (i.e. navigating directories, running scripts, etc.). If you are unfamiliar with this, here's where to start:
  - Make sure you can open a terminal.
    - **windows** - I recommend using the windows subsystem for linux (WSL) and using the default ubuntu linux distribution (https://learn.microsoft.com/en-us/windows/wsl/install). Note, I know very little about windows, which is why I'm recommending using wsl, which is basically just linux within windows. People who are knowlegable about windows may not have to do this.
    - **Mac** - you can use the default terminal app under applications/utilities
    - **Linux** - you probably already know how to use the terminal
  - do a quick tutorial to get the basics of navigating directories and running scripts.
    - This one seems good: https://ubuntu.com/tutorials/command-line-for-beginners#1-overview
    - but there are tons of others out there
- **a way to manage python environments** (e.g. conda, virtualenv, etc.)
  - I recommend using conda because it is easy to install and use. To save space I would recommend using miniconda (https://docs.conda.io/projects/miniconda/en/latest/index.html). 
      - Windows WSL - use the linux version of miniconda, which you can install via the command line (See the "Quick command line install" section on that page.)
      - here's a quick tutorial on how to use conda: https://conda.io/projects/conda/en/latest/user-guide/getting-started.html
- **the ability to edit text files**
  - You can use any editor you want. A lot of people use VSCode. For a simple and lightweight editor, I recommend sublime text to people sometimes. You can also use the default text editor on your computer (TextEdit on Mac, Notepad on Windows, etc.).

# setup TL;DR:

1. download the orthoDB database files from here: [link](https://data.orthodb.org/download/)
2. download this repository - `git clone https://github.com/jacksonh1/slim_conservation_orthogroup_generation.git`
3. navigate to this repository in terminal (where this README file is located) - `cd slim_conservation_orthogroup_generation`
4. edit the file `./orthodb_tools/env_variables/.env` with the location of the downloaded orthoDB database files: `ORTHODB_DATA_DIR=/absolute/path/to/folder/with/orthodb_files/`
5. create a new python environment with the dependencies: 
   - Mac - `conda env create -f environment.yml` <br>
       - if you have an ARM64 mac (M1/M2) you have to create an x86 environment to install all of the packages at this point in time
       - to do so run the following commands:
          ```bash
          CONDA_SUBDIR=osx-64 conda create -n slim_conservation_orthogroup_generation
          conda activate slim_conservation_orthogroup_generation
          conda config --env --set subdir osx-64
          conda env update --file=envirenment.yml --name=slim_conservation_orthogroup_generation
          ```
   - Linux/windows WSL - `conda env create -f environment_linux.yml` <br>
6. activate the environment: `conda activate slim_conservation_orthogroup_generation` <br>
7. install the local package: `pip install .` <br>
8. generate the SQLite databases: `bash ./prepare_data.sh` <br>
   - *Note: This creates separate databases for each file. You could easily make one database with all of the tables, however I tried this and it was significantly slower to query. I don't know why.* <br>

If you have issues or need more information, check out the [detailed setup instructions](setup_instructions_detailed.md).

# Usage

There are different ways to use this pipeline depending on your use case:
- You may want to run the pipeline on a small number of genes from a table
- Or you may want to run the pipeline on all of the proteins in an organism's proteome at different phylogenetic levels and use the output as a database for downstream conservation analysis. <br>

examples of both of these use cases are shown in the `./examples/` directory. There are also command line scripts to run the pipeline for these use cases in `./orthodb_tools/scripts/`. (see [useful scripts](#useful-scripts-orthodb_toolsscripts)) <br>

## overview
The main pipeline can be executed via the script: `./orthodb_tools/scripts/orthogroup_pipeline.py`. <br>

The script `./orthodb_tools/scripts/orthogroup_pipeline.py` runs the pipeline for a single gene.<br>
  - The gene can be specified using a uniprot ID or an odb_gene_id (see *brief explanation of the orthodb data* in [advanced.md](./advanced.md) for more info on the ids used in orthoDB).
  - it outputs files with the results in json format and alignments in fasta format (if alignment is specified in the parameters). <br>
The output will look like this:
```bash
main_output_folder/
    ├── info_jsons/
    │   ├── 9606_0_001c7b_Eukaryota_5471876at2759_info.json
    │   ├── ...
    ├── ├── {odb_gene_id}_{og_level}_{og_id}_info.json
    │---alignments/
    │   ├── 9606_0_001c7b_Eukaryota_5471876at2759_clustered_ldos_aln.fasta
    │   ├── ...
    ├── ├── {odb_gene_id}_{og_level}_{og_id}_clustered_ldos_aln.fasta

```
The json files contain information about the ortholog groups<br><br>
See the [examples](./examples/) folder for example output. <br>


### using the tools as a command line script

The easiest way to explain how to use the script is via the command line script help message:
```bash
$ python ./orthodb_tools/scripts/orthogroup_pipeline.py --help
```
```
usage: orthogroup_pipeline.py [-h] (-unid <str> | -odbid <str>) [-c <file>]

run main orthoDB group generation pipeline for a single gene
    processing parameters should be provided in a config file. (-c/--config)
    if no config file is provided, default parameters will be used
    The default parameters are:
- filter_params: {'min_fraction_shorter_than_query': 0.5}
- og_select_params: {'OG_selection_method': 'level_name', 'OG_level_name': 'Vertebrata'}
- ldo_select_params: {'LDO_selection_method': 'alfpy_google_distance', 'LDO_mafft_threads': 8}
- align_params: {'align': False, 'n_align_threads': 8}
- main_output_folder: ./processed_odb_groups_output
- write_files: True

options:
  -h, --help            show this help message and exit
  -unid <str>, --uniprot_id <str>
                        the uniprot id of the gene of interest
  -odbid <str>, --odb_gene_id <str>
                        the odb gene id of the gene of interest (e.g. "9606_0:001c7b")
  -c <file>, --config <file>
                        path to config file, default=None
```
- Either `-odbid` or `-unid` must be specified. <br>
  - This is the "query" sequence that the analysis is centered on (provided by either uniprot id or odb gene id). It is the sequence that you are retrieving orthologs for.<br>
- The `-c` argument is optional and specifies the location of a configuration file which can be used to specify the pipeline parameters. If it is not specified, the default parameters will be used. <br>
  - Any of the parameters can be specified in the config file and they will overwrite the defaults. Any parameters absent from the config file will take on the default value. More on what the parameters actually do [below](#pipeline-parameters) <br>
- see [example 1](./examples/ex1_single_gene/) in the examples folder to see it used as a command line script. <br>

You would probably very rarely want to run the pipeline as a command line script, but it is useful for testing how the pipeline works or if you wanted to use a workflow manager like snakemake or nextflow. <br>

### using the tools as a module

The pipeline can be imported and used from within a script, see [example 2](./examples/ex2_table_with_uniprot_ids/).

In general you import the `orthodb_tools` package and then use the `orthogroup_pipeline` function. <br>:
```python
import orthodb_tools
```
This import should work from anywhere in your filesystem because you installed the orthodb_tools code as a package, as long as the environment is activated (`conda activate slim_conservation_orthogroup_generation`). <br><br>
The main function is `orthodb_tools.orthogroup_pipeline` but it requires a configuration file to be loaded first (if using any non-default options). <br>
This is done with the `orthodb_tools.load_config` function:
```python
config = orthodb_tools.load_config("./params.yml")
```
Then you can run the pipeline with the `orthodb_tools.orthogroup_pipeline` function, which takes the `config` object as an argument and either a uniprot id or odb gene id as another argument:
```python
orthodb_tools.orthogroup_pipeline(config, odb_gene_id="9606_0:002f40")
```
or
```python
orthodb_tools.orthogroup_pipeline(config, uniprot_id="Q8TC90")
```
In this way, you can run the pipeline for any number of genes in a script as is shown in example 2 <br>

## useful scripts: `./orthodb_tools/scripts/`
There are a few scripts in the `./orthodb_tools/scripts/` directory that are useful for running the pipeline in different scenarios or provide some other common use. They are explained below: <br>

- `orthogroup_pipeline.py`: runs the main pipeline. (described above) <br>
- `pipeline_all_genes_in_species.py`: runs the pipeline for all of the proteins in an organism (in the orthoDB) at different phylogenetic levels. The levels are "Eukaryota", "Mammalia", "Metazoa", "Tetrapoda", and "Vertebrata". But you can easily change this in the script if you wanted. The levels are currently hard coded but that could easily be changed. <br>
- `pipeline_input_table.py`: Runs the pipeline for all of the proteins in a table that has a column of uniprot ids or odb gene ids. The pipeline is run for each unique gene. This is useful if you want to create a starting point for conservation analysis for just a specific set of genes (can be from different organisms as well).<br>
- `create_filemap.py`: Intended to be run after the pipeline. It creates a json file that maps the odb_gene_ids to the generated files. This is useful if you are running the pipeline on a lot of genes and you want to keep track of the files. This also creates a "database key" for use in the [motif conservation pipeline](https://github.com/jacksonh1/motif_conservation_in_IDRs)<br>
- `map_uniprotid.py`: maps uniprot ids to orthoDB gene ids in an input table. <br>
    - outputs a new table with the orthoDB gene ids added as a new column

For any of the above, you can run `python <script_name>.py --help` to see the help message. <br>

In the current version of these tools, the scripts should be accessible in your path as long as the environment in which you pip installed the tools is activated. They are installed as scripts in the environment and are accessible with the `odb_groups-` prefix. For example:

| script                             | command                                    |
| ---------------------------------- | ------------------------------------------ |
| `orthogroup_pipeline.py`           | `odb_groups-orthogroup_pipeline`           |
| `pipeline_all_genes_in_species.py` | `odb_groups-pipeline_all_genes_in_species` |
| `pipeline_input_table.py`          | `odb_groups-pipeline_input_table`          |
| `create_filemap.py`                | `odb_groups-create_filemap`                |
| `map_uniprotid.py`                 | `odb_groups-map_uniprotid`                 |

<br>
If not, you can add the `./orthodb_tools/scripts/` directory to your PATH. <br>
example of how to add this directory to your path via your bashrc file:
```bash
echo 'export PATH=$PATH:/path/to/this/repo/orthodb_tools/scripts/' >> ~/.bashrc
```
restart your terminal or run `source ~/.bashrc` to make the changes take effect. <br>

## pipeline parameters

The pipeline parameters are specified in a yaml file. <br>
If you are not familiar with yaml, it is fairly easy to understand by just looking at an example. Here is an example yaml file:<br>
```yaml
overwrite: false

filter_params:
  min_fraction_shorter_than_query: 0.5

og_select_params:
  OG_selection_method: level_name
  OG_level_name: Vertebrata

ldo_select_params:
  LDO_selection_method: alfpy_google_distance

align_params:
  align: false
  n_align_threads: 8

main_output_folder: processed_odb_groups_output
write_files: true
```
More configuration files are used in the examples here so you can also just copy and edit those for your own use (*e.g.* `./examples/ex1_single_gene/params.yml`).

### pipeline parameters explained
Here is an explanation of the parameters and what they do:
- `overwrite`: whether or not to overwrite the output files if they already exist. Can be one of:
  - true: overwrite the output files
  - false: do not overwrite the output files
- `filter_params`:
  - `min_fraction_shorter_than_query`: A number between 0 and 1. Sequences that are shorter than `min_fraction_shorter_than_query`*(`length of query sequence`) will be removed from the group of sequences. <br>Default=0.5
- `og_select_params`:
  - `OG_selection_method`: the method for selecting the ortholog groups. Can be one of:
    - `level_name`: (Default) selects the ortholog groups at a specified taxonomic level
    - there are no other options currently. I used to have options to select the group based on the number of species in the group, but I didn't think it was very useful. I've left it open to add more options in the future.
- `ldo_select_params`:
  - `LDO_selection_method`: the method for selecting the LDOs. Can be one of:
    - `alfpy_google_distance`: (Default) uses the alfpy package to calculate the google distance between the query sequence and each ortholog sequence. The ortholog sequence with the smallest google distance is selected as the LDO.
    - `msa`: uses the mafft package to construct 1 multiple sequence alignment all of the ortholog sequences. Then the paralog in each organism with the highest percent identity (PID) to the query sequence is chosen as the LDO.
    - `msa_by_organism`: Uses the mafft package to construct a separate alignment for each organism in the group, composed of the paralogs in each organism and the query sequence. Then the paralog in each organism with the highest PID to the query sequence is chosen as the LDO.
    - `pairwise`: Performs a pairwise alignment between the query sequence and each individual ortholog in the group using BioPython. Selects the paralog in each organism with the highest PID to the query sequence is chosen as the LDO.
    - `LDO_msa_threads`: the number of threads to use for the msa. Default=8. Only used if `LDO_selection_method` is `msa` or `msa_by_organism`.
- `align_params`:
  - `align`: whether or not to align the clustered LDOs. Can be one of:
    - true: (Default) align the LDOs
    - false: do not align the LDOs
  - `n_align_threads`: the number of threads to use for alignment. Default=8
- `main_output_folder`: the folder to write the output files to. Default=`./processed_odb_groups_output`
- `write_files`: whether or not to write the output files. Can be one of:
  - true: (Default) write the output files
  - false: do not write the output files

See the [advanced.md](./advanced.md) file for advanced configuration options.

# comments
The main advantages of using these tools:
- the ability to query the orthoDB files quickly (using SQLite databases)
    - The orthoDB files are very large and it is impractical to access them using normal methods (e.g. loading table into python and filtering based on a column value). The key is to pre-index the columns that you are going to filter on, however, it can't be done using normal text files. So I used SQLite to create databases with indexed columns. It makes the lookup time faster by orders of magnitude. The scripts to create the databases here index the columns that are used for lookups in the pipeline.
- go straight from a query protein to a set of orthologs
    - The OrthoDB database can be used for many different analyses. My specific use case is to just retrieve orthologs for a query protein. 
    - Retrieving ortholog groups from the downloaded files requires cross-referencing multiple large files and linking different tables with ids. 
    - I've basically figured all of that out for you so that it's easy to go from a query protein to a set of orthologs.
- reducing ortholog groups to least divergent orthologs (LDOs). 
    - This simplifies a conservation analysis by removing paralogs. It crudely assumes that the most similar sequence to the query sequence is the appropriate sequence in that organism to compare with.
    - You arguably lose information doing this (paralogs), but it greatly simplifies the analysis when you are looking at conservation for a lot of proteins (too many to manually inspect all of the sequences/alignments)
- clustering
    - Clustering reduces sequence redundancy and results in a more even distribution of sequence diversity over the group, which is very helpful in a conservation analysis.
        - Imagine you had a group of 100 homologous sequences: 60 from mammals, all with > 95 % identity and 40 distributed across more distant Vertebrates. The conservation analysis would be dominated by the mammals and would not be very informative. Preclustering would collapse highly similar sequences (>90% identity by default) into one sequence, which would make the analysis more informative.

# source code
- see [orthodb_tools/readme.md](./orthodb_tools/readme.md) for more info on how the source code is structured

# tools used (and links):
- [orthoDB database](https://www.orthodb.org/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [CD-HIT](https://sites.google.com/view/cd-hit)
- [alfpy](https://github.com/aziele/alfpy)
- [BioPython](https://biopython.org/)
- and other python packages ... see `./environment.yml` for a full list of dependencies

