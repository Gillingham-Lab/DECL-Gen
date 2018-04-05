from DECLGen import template, codon
from DECLGen.exceptions import LibraryElementException


class Element:
    cat = None
    smiles = None
    index = None

    def __init__(self, cat, smiles: str, index: int):
        self.cat = cat
        self.smiles = template.sanitize(smiles)
        self.parsed_smiles = template.parse(self.smiles)

        if not template.is_valid(self.parsed_smiles, template.get_anchors(self.smiles)):
            raise LibraryElementException("Smiles for this element is not valid.")

        self.index = index

    @property
    def codon(self):
        return codon.encode(self.index, self.cat.get_codon_length())
