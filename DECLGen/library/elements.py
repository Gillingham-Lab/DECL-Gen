import DECLGen.library

from DECLGen import template, codon
from DECLGen.exceptions import LibraryElementException


class Element:
    """
    Represents a diversity element. A element always belongs to a category.
    """
    cat: 'DECLGen.library.Category'
    raw_smiles: str
    parsed_smiles: str
    index: int

    def __init__(self, cat: 'DECLGen.library.Category', smiles: str, index: int):
        if not isinstance(cat, DECLGen.library.Category):
            raise ValueError("category must be DECLGen.library.Category")

        if not isinstance(index, int):
            raise ValueError("index must be an integer")

        if index < 0:
            raise ValueError("index must be >= 0")

        self.cat = cat
        self.user_smiles = smiles
        self.raw_smiles = template.sanitize(smiles)
        self.parsed_smiles = template.parse(self.raw_smiles)

        if not template.is_valid(self.parsed_smiles, template.get_anchors(self.raw_smiles)):
            raise LibraryElementException("Smiles for this element is not valid.")

        self.index = index

    @property
    def codon(self):
        return codon.encode(self.index, self.cat.get_codon_length())