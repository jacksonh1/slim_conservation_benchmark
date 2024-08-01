# source code structure


## `env_variables/`
`env_variables.py`:
- reads the environment variables from the `.env` file
- has some data importing functions and stores the main data files as objects
- This is imported in many of the other orthodb_tools files

## `config/`
- contains code for structuring the pipeline parameters and setting their default and allowed values

## `orthogroup_processing/`
- contains the main code for each individual step of the pipeline, stored in separate files. the `pipeline.py` file contains the main functions for running the entire pipeline and importing the config file.

## `tools/`
- contains some basic tools for working with sequences which are used throughout the pipeline

## `scripts/`
- contains scripts with CLIs for running the pipeline from various inputs
- also contains some scripts for generating a file that maps the output files to the odb gene ids


