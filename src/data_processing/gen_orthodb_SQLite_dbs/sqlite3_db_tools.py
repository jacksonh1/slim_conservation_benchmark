import csv
import sqlite3
from pathlib import Path


def create_sqlitedb_from_csv(
    cursor: sqlite3.Cursor,
    connection: sqlite3.Connection,
    csv_file_name: str | Path,
    table_name: str,
    column_names: list[str],
    index_column_name: str | list[str] | None = None,
):
    """create a sqlite table from a csv file

    Parameters
    ----------
    cursor : sqlite3.Cursor
        cursor to the database
    connection : sqlite3.Connection
        connection to the database
    csv_file_name : str | Path
        csv file to import
    table_name : str
        name of the table to create
    column_names : list[str]
        column names of the table
    index_column_name : str | list[str] | None, optional
        columns to index for faster queries, by default None
    """    
    print(f"creating sqlite database from: {csv_file_name}")
    print(f"table name: {table_name}")
    print(f"column names: {column_names}")
    print(f"indexed columns: {index_column_name}")
    print(f"text sent to sqlite3:")
    # import the csv file
    csv_file = open(csv_file_name, "r")
    csv_reader = csv.reader(csv_file, delimiter="\t")
    # Table Definition
    create_table = f"CREATE TABLE IF NOT EXISTS {table_name} (ind INTEGER PRIMARY KEY AUTOINCREMENT,"
    for column_name in column_names:
        create_table += f" {column_name} TEXT,"
    create_table = create_table[:-1] + ")"
    print(create_table)
    cursor.execute(create_table)
    insert_line = f"INSERT INTO {table_name} ("
    for column_name in column_names:
        insert_line += f"{column_name},"
    insert_line = insert_line[:-1] + ") VALUES ("
    for column_name in column_names:
        insert_line += "?,"
    insert_line = insert_line[:-1] + ")"
    print(insert_line)
    for line in csv_reader:
        cursor.execute(insert_line, line)
    csv_file.close()
    if index_column_name is not None:
        # if index_column_name is not a list, make it a list
        if isinstance(index_column_name, str):
            index_column_name = [index_column_name]
        for column_name in index_column_name:
            create_index = f"CREATE INDEX IF NOT EXISTS idx_{column_name} ON {table_name} ({column_name})"
            print(create_index)
            cursor.execute(create_index)
    connection.commit()
