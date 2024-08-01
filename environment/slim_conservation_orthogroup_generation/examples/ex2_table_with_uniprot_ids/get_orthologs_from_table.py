import orthodb_tools
import pandas as pd


def main(config, uniprot_id_list):
    for uniprot_id in uniprot_id_list:
        try:
            orthodb_tools.orthogroup_pipeline(config, uniprot_id=uniprot_id)
        except ValueError as err:
            print(f"{uniprot_id} - {err}")


if __name__ == "__main__":
    config = orthodb_tools.load_config("./params.yml")
    df = pd.read_csv("./table.csv")
    uniprot_id_list = df["Uniprotid"].unique().tolist()
    main(config, uniprot_id_list)
