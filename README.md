Written by Jackson Halpin <br>

This work was supported by the National Institutes of Health under Award Number R35GM149227. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

---

This repo holds the code and results of the benchmark analysis for [pairk](https://github.com/jacksonh1/pairk). It is mainly for reproducibility purposes. If you would like to reproduce this analysis, note that it requires quite a bit of computational resources, uses 3 separate virtual environments and generates ~100gb of data.

below are instructions for reproducing the analysis


## environments

Building the benchmark requires IUPRED2A. You can download IUPRED2A from [here](https://iupred2a.elte.hu/download_new).<br>
You then have to add the path to the downloaded iupred folder to the file `./src/data_processing/benchmark_generation/build_table/iupred_tools.py` as the variable `IUPRED_DIR`.<br>
I cannot release iupred here because it is not open source (however it is free for academic users)
<br>

3 different python environments were used to perform the benchmark analysis. I decided to keep the tools separate so that they remain uncoupled, however, this makes reproducing the benchmark more annoying. Be careful to set up the environments correctly and read the readme files in the `src` subdirectories for more information.
They are as follows:

### slim_conservation_orthogroup_generation
This environment was used to generate the orthogroups. The code and instructions are available here: [slim_conservation_orthogroup_generation](https://github.com/jacksonh1/orthogroup_generation). The full OrthoDB (v11) data files were used.


### slim_conservation_scoring
This environment was used to score the conservation using the generated orthogroup sequences. The code and instructions are available here: [slim_conservation_scoring](https://github.com/jacksonh1/slim_conservation_scoring).

### kibby
For the kibby conservation scores, an environment created from the [kibby](https://github.com/esbgkannan/kibby) repo was used and the `conservation_from_fasta.py` script was executed to get conservation scores.
<br>
See the included links for more about how the tools work


## data processing
data download and processing scripts are in `./src/data_processing/`<br>
To get started, see [./src/data_processing/README.md](./src/data_processing/README.md) for a general outline of the steps.
