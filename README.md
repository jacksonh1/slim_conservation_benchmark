Written by Jackson Halpin <br>

This work was supported by the National Institutes of Health under Award Number R35GM149227. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.


This repo holds the code and results of the benchmark analysis for [pairk](https://github.com/jacksonh1/pairk). It is mainly for reproducibility purposes. If you would like to reproduce this analysis, note that it requires quite a bit of computational resources, uses 3 separate virtual environments and generates ~100gb of data.

below are instructions for reproducing the analysis

data download and processing scripts are in `./src/data_processing/`<br>

Building the benchmark requires my [slim_conservation_orthogroup_generation](https://github.com/jacksonh1/orthogroup_generation) tool and IUPRED2A. You can download IUPRED2A from [here](https://iupred2a.elte.hu/download_new).<br>
You then have to add the path to the downloaded iupred folder to the file `./src/data_processing/benchmark_generation/build_table/iupred_tools.py` as the variable `IUPRED_DIR`.<br>
I cannot release iupred here because it is not open source (however it is free for academic users)

For scoring the benchmark and further analysis, my [slim_conservation_scoring](https://github.com/jacksonh1/slim_conservation_scoring) was used. For the kibby scores, an environment created from the [kibby](https://github.com/esbgkannan/kibby) repo was used and the `conservation_from_fasta.py` was executed to get conservation scores.

See the included links for more about how the tools work


I decided to keep the tools separate so that they remain uncoupled, however, this makes reproducing the benchmark more annoying. Be careful to set up the environments correctly and read the readme files in the `src` subdirectories for more information. 
To get started, see [./src/data_processing/README.md](./src/data_processing/README.md) for a general outline of the steps.