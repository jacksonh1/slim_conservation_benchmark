#!/usr/bin/env python

import argparse
from attrs import asdict

from orthodb_tools.config import orthodb_pipeline_parameters


import orthodb_tools.orthogroup_processing.pipeline as pipeline


if __name__ == "__main__":
    # get the default parameters just to print them in the help message
    # this is a bit hacky but it works
    # filter out the private attributes (those that start with '_') because they
    # are more advanced and probably won't be used by most users
    d_params = ""
    for k, v in asdict(
        orthodb_pipeline_parameters.PipelineParams(),
        filter=lambda attr, value: not str(attr.name).startswith("_"),
    ).items():
        d_params += f"- {k}: {v}\n"

    parser = argparse.ArgumentParser(
        description=f"""run main orthoDB group generation pipeline for a single gene
    processing parameters should be provided in a config file. (-c/--config)
    if no config file is provided, default parameters will be used
    The default parameters are:
{d_params}""",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument(
        "-unid",
        "--uniprot_id",
        type=str,
        metavar="<str>",
        help="the uniprot id of the gene of interest",
    )
    group.add_argument(
        "-odbid",
        "--odb_gene_id",
        type=str,
        metavar="<str>",
        help='the odb gene id of the gene of interest (e.g. "9606_0:001c7b")',
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        metavar="<file>",
        default=None,
        help="""path to config file, default=None""",
    )
    args = parser.parse_args()
    config = pipeline.load_config(args.config)
    pipeline.orthogroup_pipeline(config, args.uniprot_id, args.odb_gene_id)
