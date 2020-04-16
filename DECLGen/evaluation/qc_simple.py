from typing import Tuple, Optional, List
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .codon import extract
from .metadata import ReadfileMetadata, ReadfileWorkerMetadata

def compare_sequence(seq1: Seq, seq2: Seq, start: int = 0, stop: Optional[int] = None):
    return seq1[start:stop] == seq2[start:stop]


def _qc_helper(read: Seq, r: ReadfileMetadata, n: int) -> bool:
    has_passed = True

    # Check beginning of read
    #if compare_sequence(read, r.template, 0, n) is False:
    #    has_passed = False

    # Check around codon +- n
    for positions in r.coordinates:
        # Check upstream of codon if enough positions
        if positions[0] - n >= 0 and compare_sequence(read, r.template, positions[0] - n, positions[0]) is False:
            has_passed = False

        # Check downstream of codon
        if positions[1] + n < len(r.template) and compare_sequence(read, r.template, positions[1], positions[1] + n) is False:
            has_passed = False

    return has_passed


def qc(
        reads: Tuple[SeqRecord, SeqRecord],
        metadata: ReadfileWorkerMetadata
) -> Tuple[bool, bool, Tuple[List[Seq], Optional[List[Seq]]]]:
    """
    "Simple" quality control.
    :param reads:
    :param metadata:
    :return:
    """
    read_1, read_2 = reads
    r1 = metadata.r1
    r2 = metadata.r2

    if metadata.quality is None:
        quality = 3
        print("Set quality number to 3")
    else:
        try:
            quality = int(round(metadata.quality)) or 3
        except TypeError:
            quality = 3

    # Check quality
    r1_pass = _qc_helper(read_1.seq, r1, quality)
    r2_pass = True

    if metadata.is_paired():
        r2_pass = _qc_helper(read_2.seq, r2, quality)

    # Extract codons
    codons = (
        extract(read_1.seq, r1.coordinates, r1.is_reverse()),
        extract(read_2.seq, r2.coordinates, r2.is_reverse()) if metadata.is_paired() else None
    )

    return r1_pass, r2_pass, codons

qc.__name__ = "simple"