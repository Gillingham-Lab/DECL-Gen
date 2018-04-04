import math
from typing import List, Dict
from DECLGen.exceptions import LibraryCategoryException

class Category:
    id = None
    name = None
    anchors = None
    codon_length = None
    elements = None

    def __init__(self, id: str, name: str, anchors: List[str], codon_length: int = 0):
        self.id = id
        self.name = name
        self.anchors = anchors
        self.set_codon_length(codon_length)
        self.clear()

    def get_anchors(self):
        return self.anchors

    def set_name(self, name: str) -> None:
        self.name = name

    def set_codon_length(self, codon_length: int = 0) -> None:
        if codon_length < 0:
            raise LibraryCategoryException("Codon length must be a positive integer.")

        self.codon_length = codon_length

    def clear(self):
        self.elements = []

    def describe(self) -> Dict[str, str]:
        description = {
            "id": self.id,
            "name": self.name,
            "anchors": ", ".join(self.anchors),
            "codon_length": "variable" if self.codon_length == 0 else str(self.codon_length),
            "elements": len(self.elements),
            "max_elements": "unlimited" if self.codon_length == 0 else str(4**self.codon_length),
        }

        if self.codon_length == 0:
            description["min_codon_length"] = math.ceil(math.log(len(self.elements), 4))

        return description
