from typing import NewType, Union

class CodonConfig:
    """ Defines which base translate to which digit."""
    bases = "ATGC"

# Define a few types for type hinting.
bases = CodonConfig.bases
CodonInt = NewType("CodonInt", int)
CodonStr = NewType("CodonStr", str)
CodonType = NewType("CodonType", Union[CodonStr, CodonInt])

def normalize(index: CodonType) -> CodonInt:
    """
    Normalizes a CodonType to a CodonInt.

    This is a helper method that can be utilized to make methods accept both CodonInt and CodonStr. It always decodes
    a CodonStr to a CodonInt, but does not return a CodonStr. Additionally, it also tries to convert string digits to
    real digits.
    :param index: An input to convert.
    :return: The decoded codon integer.
    """
    try:
        # Try to convert to integer
        if type(index) == str and index.isdigit():
            index = int(index)
    except ValueError:
        # If it fails, we leave index alone since its a string.
        pass

    # If index is a string, we assume its a CodonStr and decode it to an integer.
    if type(index) == str:
        index = decode(index)

    return index

def reverse(codon: str) -> str:
    """
    Returns the reverse complement of a given Sequence.
    :param codon: Input sequence
    :return: Reverse complement of the input sequence.
    """
    codon = codon.upper()
    codon = codon.replace("A", "t")
    codon = codon.replace("C", "g")
    codon = codon.replace("G", "c")
    codon = codon.replace("T", "a")
    return codon.upper()[::-1]


def encode(n: int, length: int) -> str:
    """
    Encodes a number as a codon.
    :param n: int, Number to encode
    :param length: int, Target length of the encoding DNA string
    :return: str, The encoded DNA string
    """
    if n < 0:
        raise ValueError("n must be >= 0")
    if length <= 0:
        raise ValueError("length must be > 0")
    if n >= len(bases)**length:
        raise ValueError("length must be sufficient for n.")

    # Build-up codon, lowest digit first
    l = len(bases)
    codon = ""
    for i in range(length):
        codon += bases[n%l]
        n = n//l

    # Return codon, highest digit first
    return codon[::-1]


def decode(codon: str) -> int:
    """
    Decodes a codon back to the according number
    :param codon:
    :return:
    """
    if len(codon) == 0:
        ValueError("codon must not be empty")
    l = len(codon)

    # Set lowest digit first
    codon = codon[::-1]
    value = 0

    for i in range(l):
        letter = codon[i]
        val = bases.find(letter)

        if val == -1:
            raise ValueError("codon must only contain the set bases ({}, but found {})".format(bases, letter))

        value += val*4**i

    return value
