from typing import Tuple, Optional, List
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .codon import extract
from .metadata import ReadfileMetadata, ReadfileWorkerMetadata
from .qc_simple import qc as qc_simple


def compare_sequence(seq1: Seq, seq2: Seq, start: int = 0, stop: Optional[int] = None, shift: int = 0):
    return seq1[start+shift:stop+shift] == seq2[start:stop]


def _qc_helper(read: Seq, r: ReadfileMetadata, n: int, shifts) -> bool:
    has_passed = True

    # Check beginning of read
    # if compare_sequence(read, r.template, 0, n) is False:
    #    has_passed = False

    # Check around codon +- n
    a = 0
    for positions in r.coordinates:
        shift = shifts[a]

        # Check upstream of codon
        if compare_sequence(read, r.template, positions[0] - n, positions[0], shift) is False:
            has_passed = False

        # Check downstream of codon
        if compare_sequence(read, r.template, positions[1], positions[1] + n, shift) is False:
            has_passed = False

        a += 1

    return has_passed


def _codons_match(set1, set2):
    if set1 == set2:
        return True

    return False


def _get_shift_factor(x, y) -> Optional[int]:
    a = x >> (3 * y)
    if a & 0b100:
        f = -1
    else:
        f = 1
    if a & 0b010 and a & 0b001:
        f = None
    elif a & 0b010:
        f = f * 2
    elif a & 0b001:
        f = f * 1
    elif a & 0b100:
        f = None
    else:
        f = 0

    return f

def _adjust_coordinates(coordinates, shifts):
    if len(coordinates) != len(shifts):
        raise Exception("Coordinates must have as many items as shifts.")

    new_coordinates = []
    for i in range(len(coordinates)):
        new_coordinates.append((coordinates[i][0]+shifts[i], coordinates[i][1]+shifts[i]))
    return new_coordinates

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

    if metadata.is_single():
        raise Exception("The advanced qc is only possible with paired end reads.")

    read_1, read_2 = reads
    r1 = metadata.r1
    r2 = metadata.r2

    quality = int(round(metadata.quality)) or 5

    r1_pass, r2_pass, codons = qc_simple(reads, metadata)

    if _codons_match(codons[0], codons[1]) and r1_pass is True and r2_pass is True:
        return r1_pass, r2_pass, codons

    codons = len(codons[0])
    size = 3

    # First pass did not work. Lets get creativily looking for indels
    codons1 = []
    for x in range((2**size)**codons):
        shifts = []
        for y in range(codons):
            shifts.append(_get_shift_factor(x, y))

        if None in shifts:
            continue

        r1_pass = _qc_helper(read_1.seq, r1, quality, shifts)
        if r1_pass is True:
            coordinates = _adjust_coordinates(r1.coordinates, shifts)
            codons1.append(extract(read_1.seq, coordinates, r1.is_reverse()))

    codons2 = []
    for x in range((2**size)**codons):
        shifts = []
        for y in range(codons):
            shifts.append(_get_shift_factor(x, y))

        if None in shifts:
            continue

        r2_pass = _qc_helper(read_2.seq, r2, quality, shifts)
        if r2_pass is True:
            coordinates = _adjust_coordinates(r2.coordinates, shifts)
            codons2.append(extract(read_2.seq, coordinates, r2.is_reverse()))

    # Find pairs
    found = []
    for pair in codons1:
        if pair in codons2:
            found.append(pair)

    if len(found) == 1:
        return True, True, (found[0], found[0])
    elif len(found) > 1:
        return True, False, (found[0], found[-1])

    return False, False, ([], [])




qc.__name__ = "advanced"