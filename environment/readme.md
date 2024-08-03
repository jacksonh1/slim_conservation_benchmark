
# slim conservation environment
this environment is composed of 2 custom pipelines. 

Both pipelines have their own repositories, however they are included in this directory for reproducibility. Such that the exact code used in the benchmark is available here. Both tools are installed in the same environment to run the benchmark and the makefile in this directory creates a conda environment with the dependencies installed.

**slim_conservation_orthogroup_generation:**<br>
This tool was used to generate the orthogroups. The tool is available here: [slim_conservation_orthogroup_generation](https://github.com/jacksonh1/orthogroup_generation). The exact version of the code used in the benchmark is available in [./slim_conservation_orthogroup_generation](./slim_conservation_orthogroup_generation).

**slim_conservation_scoring**:<br>
This tool was used to score the conservation using the generated orthogous sequences. The code is available here: [slim_conservation_scoring](https://github.com/jacksonh1/slim_conservation_scoring). The exact version of the code used in the benchmark is available in [./slim_conservation_scoring](./slim_conservation_scoring).


To install both tools in the same environment:
- download the OrthoDB database from the OrthoDB website as described in the [slim_conservation_orthogroup_generation](./slim_conservation_orthogroup_generation) readme.md file.
    - add the directory of the downloaded files to the `./slim_conservation_orthogroup_generation/orthodb_tools/env_variables/.env` file (as described in the readme.md file)
- download iupred from the IUPred website as described in the [slim_conservation_scoring](./slim_conservation_scoring) readme.md file.
    - add the path to the `./slim_conservation_scoring/slim_conservation_scoring/env_variables/.env` file (as described in the readme.md file)
- navigate to this directory (`slim_conservation_benchmark/environment/`)
- if you use mamba (instead of conda), run this command:
    ```bash
    make environment
    ```
- if you just want to use conda, change the value of `CONDA` in the [./Makefile](./Makefile) to `conda` and run this command:
    ```bash
    make environment
    ```
- if you are having issues with the installation, try changing the channel priority to flexible with this command:
    ```bash
    conda config --set channel_priority flexible
    ```
    - I am not sure why this is necessary, but it worked for me
- The makefile will create a conda environment called slim_conservation with the necessary dependencies. If you are using a mac with ARM64 architecture (M1/M2 etc cpus), the make file should detect that and install the everything in an x86 environment.
    - note: I cannot offer full support for the ARM64 architecture. If you want to get this working on a mac without using an x86 environment this should help: the reason why an x86 environment is required is just for cd-hit and mafft in the slim_conservation_orthogroup_generation tool. If you want to build an ARM64 version of these tools, you can probably just separately install cd-hit and mafft (outside the conda environment) and then add the path to the cd-hit/mafft executables to the `./slim_conservation_orthogroup_generation/orthodb_tools/env_variables/.env` file.
- activate the environment:
    ```bash
    conda activate slim_conservation
    ```
- install the tools:
    ```bash
    pip install ./slim_conservation_orthogroup_generation
    pip install ./slim_conservation_scoring
    ```
- generate the sqlite database as described in the [slim_conservation_orthogroup_generation](./slim_conservation_orthogroup_generation) README.md file.
    - `bash ./slim_conservation_orthogroup_generation/prepare_data.sh`
- If you want to use a GPU for the ESM2 embeddings, install the required tools to access and use cuda (cuda must be installed on your system). Follow the instructions from pytorch.
    - For me, I installed the tools using this command since I have cuda 12.1 installed on my machine:
        ```bash
        conda install pytorch torchvision torchaudio pytorch-cuda=12.1 -c pytorch -c nvidia
        ```




# kibby environment

Some of the benchmark uses the kibby tool, follow the instructions at the kibby [repository](https://github.com/esbgkannan/kibby) to set up the environment. the `conservation_from_fasta.py` script from that repo was executed to get conservation scores.


<!-- - run these commands: -->
<!--     ```bash -->
<!--     conda create --name slim_conservation -->
<!--     conda activate slim_conservation -->
<!--     conda env update --name slim_conservation --file ./slim_conservation_orthogroup_generation/environment.yml -->
<!--     pip install ./slim_conservation_orthogroup_generation -->
<!--     conda env update --name slim_conservation --file ./slim_conservation_scoring/environment.yml -->
<!--     pip install ./slim_conservation_scoring -->
<!--     ``` -->
<!-- - generate the sqlite database as described in the slim_conservation_orthogroup_generation [readme.md](./slim_conservation_orthogroup_generation/README.md) file. -->
<!--     - `bash ./slim_conservation_orthogroup_generation/prepare_data.sh` -->

