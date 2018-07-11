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
    if compare_sequence(read, r.template, 0, n) is False:
        has_passed = False

    # Check around codon +- N for codon length L
    for positions in r.coordinates:
        L = positions[1] - positions[0]

        # Check upstream of codon
        if compare_sequence(read, r.template, positions[0] - L, positions[1] - L) is False:
            has_passed = False

        # Check downstream of codon
        if compare_sequence(read, r.template, positions[0] + L, positions[1] + L) is False:
            has_passed = False

    return has_passed


def qc(
        reads: Tuple[SeqRecord, SeqRecord],
        metadata: ReadfileWorkerMetadata
) -> Tuple[bool, bool, Tuple[List[Seq], Optional[List[Seq]]]]:
    read_1, read_2 = reads
    r1 = metadata.r1
    r2 = metadata.r2

    # Check quality
    r1_pass = _qc_helper(read_1.seq, r1, metadata.compare_n)
    r2_pass = True

    if metadata.is_paired():
        r2_pass = _qc_helper(read_2.seq, r2, metadata.compare_n)

    # Extract codons
    codons = (
        extract(read_1.seq, r1.coordinates, r1.is_reverse()),
        extract(read_2.seq, r2.coordinates, r2.is_reverse()) if metadata.is_paired() else None
    )

    return r1_pass, r2_pass, codons

qc.__name__ = "simple"