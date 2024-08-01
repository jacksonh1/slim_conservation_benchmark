import orthodb_tools.env_variables.env_variables as env
import orthodb_tools.sql_queries as sql_queries


def uniprotid_2_odb_gene_id(
    uniprotid: str,
    duplicate_action: str = "longest",
) -> str:
    """Given a uniprot id, return the gene id (orthodb id). If multiple gene ids are found, return the `duplicate_action` one (first or longest)

    Parameters
    ----------
    uniprotid : str
        uniprot id of the protein
    duplicate_action : str, optional
        What to do when multiple gene ids (orthodb ids) are found for a given uniprot id, by default "longest"
            if "longest", return the gene id with the longest sequence
            if "first", return the first gene id found in the database

    Returns
    -------
    str
        the gene id (orthodb id) corresponding to the provided uniprot id

    Raises
    ------
    ValueError
        if duplicate_action is not "first" or "longest"
    ValueError
        if the uniprot id is not found in the database
    """
    if duplicate_action not in ["first", "longest"]:
        raise ValueError(
            f"duplicate_action must be 'first' or 'longest', not {duplicate_action}"
        )
    # search the gene refs table
    odb_gene_ids = sql_queries.uniprotid_2_odb_gene_id_refs(uniprotid)
    if len(odb_gene_ids) == 0:
        print(f"{uniprotid} not found in gene key table, searching in xref table")
        # search the gene xref table
        odb_gene_ids = sql_queries.uniprotid_2_odb_gene_id_xrefs(uniprotid)
    if len(odb_gene_ids) == 0:
        print("not found in xref key table")
        raise ValueError(
            f"Uniprot id `{uniprotid}` not found in human gene key or xref tables"
        )
    if len(odb_gene_ids) > 1:
        print(f"Multiple matches for `{uniprotid}` in table")
        print(odb_gene_ids)
        print(f'choosing a single id from multiple matches by "{duplicate_action}"')
        if duplicate_action == "first":
            query_odb_gene_id = odb_gene_ids[0]
            return query_odb_gene_id
        else:
            data_all_seqrecords_dict = env.load_data_all_odb_seqs()
            seq_list = []
            for odb_gene_id in odb_gene_ids:
                if odb_gene_id not in data_all_seqrecords_dict:
                    print(
                        f"{odb_gene_id} not found in fasta file. Probably doesn't have an orthogroup. Skipping."
                    )
                    continue
                seq = data_all_seqrecords_dict[odb_gene_id]
                seq_list.append(seq)
            sorted_seq_list = sorted(seq_list, key=lambda x: len(x.seq), reverse=True)
            query_odb_gene_id = sorted_seq_list[0].id
            print(f"choosing {query_odb_gene_id} as the longest sequence")
            return query_odb_gene_id
    else:
        query_odb_gene_id = odb_gene_ids[0]
        return query_odb_gene_id
