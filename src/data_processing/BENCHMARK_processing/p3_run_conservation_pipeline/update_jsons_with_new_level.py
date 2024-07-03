# %%
import copy
from local_conservation_analysis_pipeline.conservation_pipeline import load_config
import json
import pandas as pd
from pathlib import Path


def update_json(output_dict: dict, database_key: dict, query_gene_id: str):
    # if query_gene_id not in database_key, this will fail and stop the whole analysis
    # we want this to happen, since this is only supposed to be run if you are updating the
    # database and don't want to redo the whole analysis. 
    database_files = database_key[query_gene_id]
    levels = [k for k in database_files.keys() if k != "query_uniprot_id"]
    for level in levels:
        if level not in output_dict["orthogroups"]:
            output_dict["orthogroups"][level] = copy.deepcopy(database_files[level])
            if "conservation_scores" not in output_dict["orthogroups"][level]:
                output_dict["orthogroups"][level]["conservation_scores"] = {}
    return output_dict


def add_level_to_analysis_from_database_key(indexed_hits_file, database_key_file, output_folder):
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    hits_df = pd.read_csv(indexed_hits_file)
    with open(database_key_file, "r") as f:
        database_key = json.load(f)

    for i, row in hits_df.iterrows():
        ref_ind = int(row["reference_index"])
        query_gene_id = row["gene_id"]
        analysis_folder = output_folder / f'{ref_ind}-{str(query_gene_id).replace(":","_")}'
        analysis_folder.mkdir(parents=True, exist_ok=True)
        output_file = (
            analysis_folder / f'{ref_ind}-{str(query_gene_id).replace(":","_")}.json'
        )
        if output_file.exists():
            with open(output_file, "r") as f:
                output_dict = json.load(f)
            if 'critical_error' in output_dict:
                print(f"Critical error in {output_file}")
                continue
            try:
                output_dict = update_json(output_dict, database_key, query_gene_id)
            except Exception as e:
                print(e)
                continue
            with open(output_file, "w") as f:
                json.dump(output_dict, f, indent=4)


def main():
    config = load_config(config_file='./params.yaml')
    table_filename = Path(config.table_file).name
    reindexed_table_file = Path(config.output_folder) / table_filename.replace(
        ".csv", "_original_reindexed.csv"
    )
    add_level_to_analysis_from_database_key(reindexed_table_file, config.database_filekey, config.output_folder)



if __name__ == "__main__":
    main()