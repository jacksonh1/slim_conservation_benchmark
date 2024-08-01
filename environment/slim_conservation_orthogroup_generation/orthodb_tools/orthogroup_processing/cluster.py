import copy
from typing import Union

import pandas as pd
from Bio import SeqIO

import orthodb_tools.tools.cdhit_tools as cdhit_tools
import orthodb_tools.tools.cli_wrappers as cli


def cdhit_clstr_retrieve_representative_sequences(
    clstr_dict: dict, seqrecord_dict: dict[str, SeqIO.SeqRecord]
) -> dict[str, SeqIO.SeqRecord]:
    """
    pull out representative seqs defined in cdhit clstr_dict from full seqrecord_dict
    """
    clustered_seq_dict = {}
    for cluster_id in clstr_dict.keys():
        id_i = clstr_dict[cluster_id]["representative_seq"]
        # id_i = re.findall(r'\d+\_\d\:.+$', rep_i)[0]
        clustered_seq_dict[id_i] = copy.deepcopy(seqrecord_dict[id_i])
    return clustered_seq_dict


def cdhit_main(
    seqrecord_dict: dict[str, SeqIO.SeqRecord],
    query_odb_gene_id: str,
    repr_id_keywords: list[str] | None = None,
    **kwargs,
) -> tuple[str, dict[str, SeqIO.SeqRecord]]:
    """ """
    if repr_id_keywords is None:
        repr_id_keywords = []
    if query_odb_gene_id not in repr_id_keywords:
        repr_id_keywords.append(query_odb_gene_id)

    cdhit_command, _, cdhit_clstr_dict = cli.cd_hit_wrapper(
        list(seqrecord_dict.values()), **kwargs
    )
    cdhit_clstr_dict = (
        cdhit_tools.cd_hit_clstr_redefine_cluster_representative_by_keywords(
            cdhit_clstr_dict, repr_id_keywords
        )
    )

    clustered_seqrecord_dict = cdhit_clstr_retrieve_representative_sequences(
        cdhit_clstr_dict, seqrecord_dict
    )
    return cdhit_command, clustered_seqrecord_dict
