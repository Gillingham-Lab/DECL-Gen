from typing import Tuple, Optional, List
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from .metadata import ReadfileMetadata, ReadfileWorkerMetadata
from .codon import extract
from ..template import get_codon_coordinates


def _match_score(x, y):
    if y == "N":
        return 1
    elif y == x:
        return 5
    else:
        return -2

class Tmp:
    N = None
    Max = None
    Length = None

def _qc_helper(read: Seq, r: ReadfileMetadata, f: float = 0.3) -> Tuple[bool, Optional[List[Seq]]]:
    """
    Helper method for qc_align()
    :param read:
    :param r:
    :param f:
    :return:
    """
    #alignment = pairwise2.align.localms(read, r.template, 5, -3, -4, -4, one_alignment_only=True)[0]exit

    if Tmp.N is None:
        Tmp.N = r.template.count("N")

    if Tmp.Max is None:
        Tmp.Max = (len(r.template)-Tmp.N)*5 + Tmp.N

    if Tmp.Length is None:
        Tmp.Length = len(r.template)

    read = read[0:Tmp.Length]

    alignment = pairwise2.align.localcs(read, r.template, _match_score, -10, -1, one_alignment_only=True)[0]

    codons = extract(alignment[0], get_codon_coordinates(alignment[1]), r.is_reverse())
    indels = sum([1 for codon in codons if codon.count("-") > 0])

    if indels > 0:
        return False, None

    if alignment[2] > (Tmp.Max - 5*7 - 1):
        return True, codons
    else:
        return False, None

def qc(
        reads: Tuple[SeqRecord, SeqRecord],
        metadata: ReadfileWorkerMetadata
) -> Tuple[bool, bool, Tuple[List[Seq], Optional[List[Seq]]]]:
    """
    Alignment-based quality control and codon extraction.
    :param reads:
    :param metadata:
    :return:
    """
    read_1, read_2 = reads
    r1 = metadata.r1
    r2 = metadata.r2
    quality = 0.3
    if metadata.quality is not None:
        quality = float(abs(metadata.quality))
    quality = quality if quality < 1 else 1

    # Check quality and fetch adjusted coordinates
    r1_pass, r1_codons = _qc_helper(read_1.seq, r1, quality)
    r2_pass, r2_codons = True, None

    if metadata.is_paired():
        r2_pass, r2_codons = _qc_helper(read_2.seq, r2, quality)

    # Extract codons
    codons = (
        r1_codons,
        r2_codons
    )

    return r1_pass, r2_pass, codons

qc.__name__ = "align"
