from typing import Tuple, Optional, List
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from .metadata import ReadfileMetadata, ReadfileWorkerMetadata
from .codon import extract
from ..template import get_codon_coordinates

def _qc_helper(read: Seq, r: ReadfileMetadata, f: float = 0.5) -> bool:
    score = pairwise2.align.localms(read, r.template, 5, -3, -4, -4, score_only=True)
    max_score = min(len(read), len(r.template)) * 5

    if score > max_score * f:
        return True

    return False

def qc(
        reads: Tuple[SeqRecord, SeqRecord],
        metadata: ReadfileWorkerMetadata
) -> Tuple[bool, bool, Tuple[List[Seq], Optional[List[Seq]]]]:
    read_1, read_2 = reads
    r1 = metadata.r1
    r2 = metadata.r2

    # Check quality
    r1_pass = _qc_helper(read_1.seq, r1)
    r2_pass = True

    if metadata.is_paired():
        r2_pass = _qc_helper(read_2.seq, r2)

    # Extract codons
    codons = (
        extract(read_1.seq, r1.coordinates, r1.is_reverse()),
        extract(read_2.seq, r2.coordinates, r2.is_reverse()) if metadata.is_paired() else None
    )

    return r1_pass, r2_pass, codons

qc.__name__ = "score"