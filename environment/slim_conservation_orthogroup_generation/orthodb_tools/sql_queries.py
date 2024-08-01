import sqlite3
from pathlib import Path

from orthodb_tools.env_variables import env_variables as env


def uniprotid_2_odb_gene_id_refs(
    uniprotid, db_path: str | Path = env.orthoDB_files.gene_refs_sqlite
) -> list[str]:
    """return the odb_gene_id from a uniprot ID"""
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()
    res = cursor.execute(
        f"SELECT odb_gene_id FROM gene_refs WHERE Uniprotid='{uniprotid}'"
    )
    odb_gene_ids = res.fetchall()
    odb_gene_ids = [x[0] for x in odb_gene_ids]
    connection.close()
    return odb_gene_ids


def uniprotid_2_odb_gene_id_xrefs(
    uniprotid, db_path: str | Path = env.orthoDB_files.gene_xrefs_sqlite
) -> list[str]:
    """return the odb_gene_id from a uniprot ID"""
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()
    # res = cursor.execute(f"SELECT odb_gene_id FROM gene_xrefs WHERE xref_id='{uniprotid}' AND DB_name='UniProt'")
    # This is much faster if you don't specify the DB_name and then filter the results
    res = cursor.execute(f"SELECT * FROM gene_xrefs WHERE xref_id='{uniprotid}'")
    odb_gene_ids = res.fetchall()
    odb_gene_ids = [x[1] for x in odb_gene_ids if x[3] == "UniProt"]
    connection.close()
    return odb_gene_ids


def odb_gene_id_2_species_id(
    odb_gene_id, db_path: str | Path = env.orthoDB_files.gene_refs_sqlite
) -> str:
    """return the species ID from an orthodb ID"""
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()
    res = cursor.execute(
        f"SELECT species_id FROM gene_refs WHERE odb_gene_id='{odb_gene_id}'"
    )
    species_id = res.fetchall()[0][0]
    connection.close()
    return species_id


def odb_gene_id_2_ogid_list(
    odb_gene_id, db_path: str | Path = env.orthoDB_files.OG2genes_sqlite
) -> list[str]:
    """Given a odb_gene_id, return the list of OGs it belongs to

    Parameters
    ----------
    odb_gene_id : str
        gene id (orthodb id)
    db_path : str|Path, optional
        path to the OG2genes sqlite database, by default env.orthoDB_files.OG2genes_sqlite

    Returns
    -------
    list[str]
        list of OGs the gene id (orthodb id) belongs to
    """
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()
    res = cursor.execute(
        f"SELECT OG_id FROM OG2genes WHERE odb_gene_id='{odb_gene_id}'"
    )
    og_ids = res.fetchall()
    og_ids = [og_id[0] for og_id in og_ids]
    og_ids = list(set(og_ids))
    connection.close()
    if len(og_ids) == 0:
        raise ValueError(f"no OGs found for gene id {odb_gene_id}")
    return og_ids


def get_ogid_info(
    ogid, db_path: str | Path = env.orthoDB_files.ogs_sqlite
) -> tuple[str]:
    """query the OGs sqlite database for the info about the OG with the given ogid

    Returns
    -------
    tuple[str]
        returns a tuple composed of (ogid, level NCBI tax id, and OG name)
    """
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()
    res = cursor.execute(f"SELECT * FROM OGs WHERE OG_id='{ogid}'")
    og_info = res.fetchall()[0]
    # raise error if no results found?
    connection.close()
    return og_info[1:]


def odb_gene_id_2_uniprotid(
    odb_gene_id, db_path: str | Path = env.orthoDB_files.gene_refs_sqlite
) -> str:
    """return the uniprot ID from an orthodb ID"""
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()
    res = cursor.execute(
        f"SELECT Uniprotid FROM gene_refs WHERE odb_gene_id='{odb_gene_id}'"
    )
    result = res.fetchall()
    connection.close()
    if len(result) == 0:
        # raise ValueError(f"no uniprot id found for gene id {odb_gene_id}")
        return ""
    uniprot_id = result[0][0]
    return uniprot_id


def ogid_2_odb_gene_id_list(
    ogid, db_path: str | Path = env.orthoDB_files.OG2genes_sqlite
) -> list[str]:
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()
    res = cursor.execute(f"SELECT odb_gene_id FROM OG2genes WHERE OG_id='{ogid}'")
    odb_gene_ids = res.fetchall()
    odb_gene_ids = [x[0] for x in odb_gene_ids]
    connection.close()
    return odb_gene_ids


def get_all_odb_gene_ids_from_species_id(
    species_id: str, db_path: str | Path = env.orthoDB_files.gene_refs_sqlite
) -> list[str]:
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()
    res = cursor.execute(f"SELECT * FROM gene_refs WHERE species_id='{species_id}'")
    results = res.fetchall()
    gene_list = list(set([i[1] for i in results]))
    connection.close()
    return gene_list
