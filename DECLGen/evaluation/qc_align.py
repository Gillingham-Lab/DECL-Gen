from typing import Tuple, Optional, List
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from .metadata import ReadfileMetadata, ReadfileWorkerMetadata
from .codon import extract
from ..template import get_codon_coordinates

def _qc_helper(read: Seq, r: ReadfileMetadata, f: float = 0.2) -> Tuple[bool, Optional[List[Seq]]]:
    """
    Helper method for qc_align()
    :param read:
    :param r:
    :param f:
    :return:
    """
    alignment = pairwise2.align.localms(read, r.template, 5, -3, -4, -4, one_alignment_only=True)[0]
    max_score = min(len(read), len(r.template)) * 5

    codons = extract(alignment[0], get_codon_coordinates(alignment[1]), r.is_reverse())

    if alignment[2] > max_score * f:
        return True, codons
    else:
        return False, None

def qc(
        reads: Tuple[SeqRecord, SeqRecord],
        metadata: ReadfileWorkerMetadata
) -> Tuple[bool, bool, Tuple[List[Seq], Optional[List[Seq]]]]:
    """
    Quality control with true alignment
    :param reads:
    :param metadata:
    :return:
    """
    read_1, read_2 = reads
    r1 = metadata.r1
    r2 = metadata.r2

    # Check quality and fetch adjusted coordinates
    r1_pass, r1_codons = _qc_helper(read_1.seq, r1)
    r2_pass, r2_codons = True, None

    if metadata.is_paired():
        r2_pass, r2_codons = _qc_helper(read_2.seq, r2)

    # Extract codons
    codons = (
        r1_codons,
        r2_codons
    )

    return r1_pass, r2_pass, codons

qc.__name__ = "align"
