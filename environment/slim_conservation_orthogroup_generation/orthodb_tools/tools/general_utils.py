import re
from pathlib import Path

from Bio import Align, AlignIO, Seq, SeqIO
from Bio.SeqRecord import SeqRecord


def import_fasta(
    fasta_path, output_format="list"
) -> list[SeqRecord] | dict[str, SeqRecord]:
    """import fasta file into a list or dictionary of SeqRecord objects

    Parameters
    ----------
    fasta_path : str
        file path to fasta file
    output_format : str, optional
        output_format of output. Either 'list' or 'dict'. by default 'list'

    Returns
    -------
    list or dictionary
        list or dictionary of SeqRecord objects for each sequence in the fasta file
    """
    allowed_formats = ["list", "dict"]
    with open(fasta_path) as handle:
        if output_format == "list":
            seqs = list(SeqIO.parse(handle, "fasta"))
        elif output_format == "dict":
            seqs = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        else:
            raise ValueError(
                f"Invalid output format - {output_format}. Expected one of: {allowed_formats}"
            )
    return seqs


class FastaImporter:
    """import fasta file and return seqrecord objects in various formats

    Parameters
    ----------
    fasta_path : str
        file path to fasta file
    """

    def __init__(self, fasta_path: str | Path):
        self.fasta_path = fasta_path

    def import_as_list(self) -> list[SeqRecord]:
        """return list of SeqRecord objects for each sequence in the fasta file

        Returns
        -------
        List[SeqRecord]
            list of SeqRecord objects
        """
        with open(self.fasta_path) as handle:
            return list(SeqIO.parse(handle, "fasta"))

    def import_as_dict(self) -> dict[str, SeqRecord]:
        """return dictionary of SeqRecord objects for each sequence in the fasta file

        Returns
        -------
        dict[str, SeqRecord]
            dictionary of SeqRecord objects, keys are the sequence ids and values are the SeqRecord objects
        """
        with open(self.fasta_path) as handle:
            return SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    def import_as_alignment(self) -> Align.MultipleSeqAlignment:
        """return multiple sequence alignment object

        Returns
        -------
        Align.MultipleSeqAlignment
            multiple sequence alignment object
        """
        with open(self.fasta_path) as handle:
            return AlignIO.read(handle, "fasta")
