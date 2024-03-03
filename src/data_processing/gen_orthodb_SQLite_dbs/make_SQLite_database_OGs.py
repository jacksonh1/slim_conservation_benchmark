import sqlite3
from pathlib import Path

import sqlite3_db_tools as sqltools

import local_env_variables.env_variables as env

db_file_name = env.orthoDB_files.ogs_sqlite
Path(db_file_name).touch()
connection = sqlite3.connect(db_file_name)
cursor = connection.cursor()
sqltools.create_sqlitedb_from_csv(
    cursor,
    connection,
    csv_file_name=env.orthoDB_files.ogs_tsv,
    table_name="OGs",
    column_names=["OG_id", "level_NCBI_tax_id", "OG_name"],
    index_column_name="OG_id",
)
connection.close()
