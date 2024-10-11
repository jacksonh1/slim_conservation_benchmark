#!/usr/bin/env python

import argparse
from pathlib import Path

import pandas as pd

from orthodb_tools.orthogroup_processing import uniprotid_search


def map_uniprot_id(
    table: pd.DataFrame, uniprotid_column_name: str = "uniprot_id"
) -> pd.DataFrame:
    """maps a column of uniprot ids in a table to orthoDB gene ids.
    adds a new column containing the orthoDB gene ids named `gene_id`.
    if the uniprot ids have isoform information, e.g. "P12345-1", the isoform
    information will be stripped before mapping, e.g. "P12345-1" -> "P12345"."""
    assert (
        uniprotid_column_name in table.columns
    ), f"column {uniprotid_column_name} not found in input table"
    # strip leading and trailing whitespaces from uniprot ids.
    mapping_colname = uniprotid_column_name + "_stripped"
    table[mapping_colname] = table[uniprotid_column_name].str.strip()
    # strip isoform information from uniprot ids.
    table[mapping_colname] = table[mapping_colname].str.replace(
        r"-\d+$", "", regex=True
    )
    uniprot_ids = list(table[mapping_colname].unique())
    id_map = {}
    for uniprot_id in uniprot_ids:
        try:
            id_map[uniprot_id] = uniprotid_search.uniprotid_2_odb_gene_id(uniprot_id)
        except ValueError:
            print(f"uniprot id {uniprot_id} not found in orthoDB. Skipping...")
            continue
    table["gene_id"] = table[mapping_colname].map(id_map)
    table = table.drop(columns=[mapping_colname])
    return table


def main(
    input_file: str | Path,
    uniprotid_column_name: str = "uniprot_id",
    output_file: str | Path | None = None,
):
    """maps a column of uniprot ids in a table to orthoDB gene ids.
    exports a copy of the table with a new column containing the orthoDB gene ids.
    if the uniprot ids have isoform information, e.g. "P12345-1", the isoform
    information will be stripped before mapping, e.g. "P12345-1" -> "P12345"."""
    input_file = Path(input_file)
    if output_file is None:
        output_file = (
            input_file.parent / f"{input_file.stem}_mapped_odbgeneid{input_file.suffix}"
        )
    else:
        output_file = Path(output_file)
    table = pd.read_csv(input_file)
    table = map_uniprot_id(table, uniprotid_column_name)
    table.to_csv(output_file, index=False)


def main_cli():
    parser = argparse.ArgumentParser(
        description="""maps a column of uniprot ids in a table to orthoDB gene ids.
exports a copy of the table with a new column containing the orthoDB gene ids.\n
if the uniprot ids have isoform information, e.g. "P12345-1", the isoform 
information will be stripped before mapping, e.g. "P12345-1" -> "P12345".""",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        metavar="<file>",
        help="""input table with a column of uniprot ids. In csv format.""",
    )
    parser.add_argument(
        "--uni_column",
        type=str,
        default="uniprot_id",
        metavar="<str>",
        help="""name of the column containing the uniprot ids. Default: uniprot_id""",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        metavar="<file>",
        help='''output file name. Default: input file name + "_mapped_odbgeneid"''',
    )
    args = parser.parse_args()
    main(args.input, args.uni_column, args.output)


if __name__ == "__main__":
    main_cli()
