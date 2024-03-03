import sqlite3
from pathlib import Path

import sqlite3_db_tools as sqltools

import local_env_variables.env_variables as env

db_file_name = env.orthoDB_files.OG2genes_sqlite
Path(db_file_name).touch()
connection = sqlite3.connect(db_file_name)
cursor = connection.cursor()
sqltools.create_sqlitedb_from_csv(
    cursor,
    connection,
    csv_file_name = env.orthoDB_files.OG2genes_tsv,
    table_name = 'OG2genes',
    column_names = ['OG_id', 'odb_gene_id'],
    index_column_name=['OG_id', 'odb_gene_id']
)

connection.close()
