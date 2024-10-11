The purpose of these scripts are to create tables with the different conservation scores from the results of the pipeline. 

The pipeline is designed to add a few scores to a final table however, for the benchmark I want to compare all of the scores and processing parameters. So I wrote scripts that just assigns a numerical index to each score and associated parameters and saves a single table for each. I will then merge the tables together later for performance analysis and comparison. 

commands run:
```bash
conda activate slim_conservation
nohup python create_score_table_wide_form.py > create_score_table_wide_form.out &
```