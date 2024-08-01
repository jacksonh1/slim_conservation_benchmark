# Advanced stuff

## brief explanation of the orthoDB data
You can view the readme file that comes with the orthoDB download for more information. <br>
The data is organized using some internal id numbers. <br>
Here is a brief explanation of the ids and how I've refered to them in the code:
- **odb_gene_id**: An internal orthoDB id. It defines a specific protein sequence in the database. The sequences in the fasta file have the orthoDB id as the sequence id. It is not consistent across versions of orthoDB.
  - example: `9606_0:001c7b`
  - odb_gene_id's are mapped to a variety of other ids corresponding to outside databases (e.g. uniprot, ensembl, etc.) or they were downloaded from some database and have a corresponding id. This information is stored in the `odb11v0_gene_xrefs.tab`/`odb11v0_genes.tab` files.
- **og_id**: An internal orthoDB id. It is probably not consistent across versions of orthoDB. It defines a group of homologous sequences. Each phylogenetic level of homologs has a unique og_id. So a single odb_gene_id probably belongs to multiple og_ids.
  - example og_id: `1567973at7742`
  - example - all of the og_id's that contain the odb_gene_id `9606_0:001c7b`:
    - `605262at9347`, `70995at314146`, `5821at9604`, `1742826at33208`, `5821at314295`, `4349at40674`, `1005199at32523`, `1567973at7742`, `5471876at2759`, `70995at9443`
  - Each of the og_id's is associated with a phylogenetic level at which it was constructed:
      | OG id          | level name       | total non-redundant count of species underneath |
      | :------------- | :--------------- | ----------------------------------------------: |
      | 5471876at2759  | Eukaryota        |                                            1952 |
      | 1742826at33208 | Metazoa          |                                             817 |
      | 1567973at7742  | Vertebrata       |                                             470 |
      | 1005199at32523 | Tetrapoda        |                                             325 |
      | 4349at40674    | Mammalia         |                                             191 |
      | 605262at9347   | Eutheria         |                                             182 |
      | 70995at314146  | Euarchontoglires |                                              70 |
      | 70995at9443    | Primates         |                                              30 |
      | 5821at314295   | Hominoidea       |                                               7 |
      | 5821at9604     | Hominidae        |                                               5 |
      - These correspond with the groups on the website: https://www.orthodb.org/?query=9606_0%3A001c7b
      - note that the `total non-redundant count of species underneath` is the number of species in the database underneath that phylogenetic level, NOT the number of sequences in that orthogroup. OrthoDB doesn't necessarily find an ortholog in each species (nor should it, some species don't have a homolog for a given gene).

## advanced configuration

### clustering and alignment parameters
if you have MAFFT and/or CD-HIT installed somewhere else and would prefer to use versions, you can change the `MAFFT_EXECUTABLE` and `CD_HIT_EXECUTABLE` variables in the `.env` file to point to the executables on your computer. <br>

You can add also add additional command line arguments to the MAFFT and CD-HIT commands by editing the `MAFFT_ADDITIONAL_ARGUMENTS` and `CD_HIT_ADDITIONAL_ARGUMENTS` variables in the `.env` file. <br>
Those variables are set to empty strings by default, but they are inserted into the mafft/cd-hit commands where extra arguments would go:
- Mafft: `{MAFFT_EXECUTABLE} --thread {n_align_threads} --quiet --anysymbol {MAFFT_ADDITIONAL_ARGUMENTS} {input_file} > {output_alignment}`
- CD-hit: `{CD_HIT_EXECUTABLE} -i {input_file} -o {output_file} -M 0 -d 0 -g 1 {CD_HIT_ADDITIONAL_ARGUMENTS}` <br>
  - Note that for CD-HIT, the default of 90% sequence identity is used, therefore you can change the clustering % identity by providing it to the CD_HIT_ADDITIONAL_ARGUMENTS variable <br>
  - For example, if you wanted to change the clustering step to cluster the sequences to 80% identity, you would change the `CD_HIT_ADDITIONAL_ARGUMENTS` variable in the `.env` file to `-c 0.8` <br>

Setting these at the environment level is not really ideal if you want these parameters to be flexible.<br>
Therefore, you can also change the MAFFT and CD-HIT commands in the yaml config file (see `./examples/readme.md` for explanation) via some hidden parameters shown in this example:
```yaml 
ldo_select_params:
  _LDO_msa_exe: mafft
  _LDO_msa_exe_additional_args: --retree 1

align_params:
  _mafft_exe: mafft
  _mafft_exe_additional_args: --retree 1

_cd_hit_exe: cd-hit
_cd_hit_exe_additional_args: -c 0.8
```
These parameters are typically set via the variables in the .env file, but that behavior can be overwritten by the config file <br>

### incorporating different aligners

It would be fairly straightforward to incorporate different aligners into the pipeline. <br>
You would first have to add a new command line script wrapper for the aligner to `./orthodb_tools/tools/cli_wrappers.py` (There are already some unused functions in there for muscle/clustal). <br>
Then you could add a new configuration class to configure the aligner in `./orthodb_tools/config/orthodb_pipeline_parameters.py`, and alter the main pipeline script (`./orthodb_tools/orthogroup_processing/pipeline.py`) to use the new aligner. <br>
However, it might be easier to just use the pipeline to generate the ortholog groups using align=false (creates just json files) and then run an aligner outside of the pipeline. <br>
I would just write a script that imports the json file to get the sequence ids of the clustered ldo sequences and then runs the aligner on those sequences. <br>
