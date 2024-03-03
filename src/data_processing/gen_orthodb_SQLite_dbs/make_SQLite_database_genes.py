import sqlite3
from pathlib import Path

import sqlite3_db_tools as sqltools

import local_env_variables.env_variables as env

db_file_name = env.orthoDB_files.gene_refs_sqlite
Path(db_file_name).touch()
connection = sqlite3.connect(db_file_name)
cursor = connection.cursor()
sqltools.create_sqlitedb_from_csv(
    cursor,
    connection,
    csv_file_name = env.orthoDB_files.gene_refs_tsv,
    table_name = 'gene_refs',
    column_names = [
        'odb_gene_id',
        'species_id',
        'source_id',
        'synonyms',
        'Uniprotid',
        'Ensemble',
        'NCBI_id',
        'description',
    ],
    index_column_name=['odb_gene_id', 'Uniprotid']
)

connection.close()
