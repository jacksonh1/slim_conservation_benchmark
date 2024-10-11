Written by Jackson Halpin <br>

This work was supported by the National Institutes of Health under Award Number R35GM149227. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

---

This repo holds the code and results of the benchmark analysis for [pairk](https://github.com/jacksonh1/pairk). It is mainly for reproducibility purposes. If you would like to reproduce this analysis, note that it requires quite a bit of computational resources, uses 2 separate virtual environments and requires ~300Gb of data.

below are instructions for reproducing the analysis

## environments

Building the benchmark requires IUPRED2A. You can download IUPRED2A from [here](https://iupred2a.elte.hu/download_new).<br>
You then have to add the path to the downloaded iupred folder to the file `./src/data_processing/benchmark_generation/build_table/iupred_tools.py` as the variable `IUPRED_DIR`.<br>
I cannot release iupred here because it is not open source (however it is free for academic users)
<br>

3 different python tools were used to perform the benchmark analysis contained in 2 different environments. Be careful to set up the environments correctly and read the readme files in the `src/data_processing/` subdirectories for more information. to set up the environments, follow the instructions in the [./environment](./environment) readme file.

<!-- ### reproducing the benchmark
To reproduce the benchmark, you will need to set up the environments described above. To make this as reproducible as possible, I've included copies of the `slim_conservation` source code in this repo (in the [tools](./tools) directory) so that the exact version of the code used in the benchmark is available. The `kibby` code is not included here, but you can find it at the link above.
To create the exact environments used:
- `git clone http....`
- `cd [benchmark repo]/tools/slim_conservation_orthogroup_generation`
- follow the instructions in [./tools/slim_conservation_orthogroup_generation/README.md](./tools/slim_conservation_orthogroup_generation/README.md) to set up the environment and install the src code in that environment
- `cd [benchmark repo]/tools/slim_conservation_scoring`
- follow the instructions in [./tools/slim_conservation_scoring/README.md](./tools/slim_conservation_scoring/README.md) to set up the environment and install the src code in that environment -->


## data processing
data download and processing scripts are in `./src/data_processing/`<br>
To get started, see [./src/data_processing/README.md](./src/data_processing/README.md) for a general outline of the steps.
