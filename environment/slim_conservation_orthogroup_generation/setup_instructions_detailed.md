# setup (more detailed):

## Download orthoDB database files
- download the orthoDB database files from here: [link](https://data.orthodb.org/download/)
  - I used version 11.0 for my project
  - The files you need are:
    - `odb11v0_all_og_fasta.tab.gz`
    - `odb11v0_levels.tab.gz`
    - `odb11v0_species.tab.gz`
    - `odb11v0_level2species.tab.gz`
    - `odb11v0_genes.tab.gz`
    - `odb11v0_gene_xrefs.tab.gz`
    - `odb11v0_OGs.tab.gz`
    - `odb11v0_OG2genes.tab.gz`
- Unzip the files (todo: see if you can work with the compressed files instead)
- save the path to this directory so that you can use it in the `.env` file (see below)

## Download this repository and edit the `.env` file
```bash
git clone https://github.com/jacksonh1/orthogroup_generation.git
```
- Edit the `./orthodb_tools/env_variables/.env` file
  - add a `ORTHODB_DATA_DIR` variable with the absolute path to the directory containing the downloaded orthoDB files
  - example `.env` file: 
```
ORTHODB_DATA_DIR=/Users/username/project/data/orthodb/odb11v0/
```
- I have included example files in the `./data/orthoDB_sample_data/` directory. These files are from orthoDB version 11.0 and are just subsets of the full files<br>
- Note: *if you use a different version of orthoDB, you will need to change the file names in `./orthodb_tools/env_variables/env_variables.py`. File names are hard coded in the `OrthoDBFiles` class*

## conda environment
Create a conda environment with the necessary packages using the environment file. It includes CD-HIT and MAFFT from bioconda so you do not have to install these separately. <br>
`conda env create -f environment.yml` <br>

### Mac
if you have an ARM64 mac (M1/M2) you have to create an x86 environment to install all of the packages at this point in time <br>
to do so run the following commands:<br>
```bash
CONDA_SUBDIR=osx-64 conda create -n slim_conservation_orthogroup_generation
conda activate slim_conservation_orthogroup_generation
conda config --env --set subdir osx-64
conda update -f=envirenment.yml --name=slim_conservation_orthogroup_generation
```
if you have an older intel mac, you shouldn't have to use the above commands, you should be able to just run the following command: <br>
`conda env create -f environment.yml` <br>


## install local tools in environment

**TL;DR**: 
- activate the environment: `conda activate slim_conservation_orthogroup_generation`
- run `pip install .` in this directory (where `setup.py` file is located) <br>

I wrote this pipeline (code in `./orthodb_tools/`) for it to be installed in the environment as a local package. <br>
I took inspiration from this example: https://github.com/ericmjl/Network-Analysis-Made-Simple <br>
combined with information from here: https://setuptools.pypa.io/en/latest/userguide/package_discovery.html <br>

Install the modules/tools as a local package using setuptools. <br>
first activate the environment: <br>
`conda activate slim_conservation_orthogroup_generation` <br>

Make sure that you are in this directory in terminal (where `pyproject.toml` is located) and run the following command to install: <br>
`pip install .` <br>

If you want to make modifications to the orthodb_tools code, you can install it as an editable install. Then you can make changes to the orthodb_tools files and they will be updated in the environment. <br>
If you want to do this, run the following command instead: <br>
`pip install -e .` <br><br>

## generate SQLite databases for orthoDB files
In this directory, run the following command to generate SQLite databases for the orthoDB files: <br>
```bash
bash ./prepare_data.sh
```

Warning - this will take a while to run and I don't think it can be parallelized (Let me know if this is possible because I would love to know how if so). <br>

*Note: This creates separate databases for each file. You could easily make one database with all of the tables, however I tried this and it was significantly slower to query.* <br>
