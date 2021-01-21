from typing import Union, List, Tuple
from Bio.Seq import Seq


def extract(read: Union[Seq, str], coordinates: List[Tuple[int, int]], invert=False) -> List[Seq]:
    """
    Extracts codons from a given Bio.Seq.Seq object depending on the coordinates given.

    :param read: The read to extract the codons from.
    :param coordinates: A list of coordinate-tuples for all codons.
    :param invert: If true, returns the coordinates as their reverse complement.
    :return: A list of codons as strings.
    """
    codons = []
    for start, stop in coordinates:
        codon = read[start:stop]

        if type(codon) == str:
            codon = Seq(codon)

        if invert is True:
            codons.append(codon.reverse_complement())
        else:
            codons.append(codon)

    return codons if invert is False else codons[::-1]