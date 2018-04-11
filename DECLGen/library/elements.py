from DECLGen import template, codon
from DECLGen.exceptions import LibraryElementException


class Element:
    """Represents a diversity element. A element always belongs to a category"""
    cat = None
    raw_smiles = None
    parsed_smiles = None
    index = None

    def __init__(self, cat, smiles: str, index: int):
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