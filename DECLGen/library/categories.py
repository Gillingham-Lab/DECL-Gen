import math
from typing import List, Dict, Union
from DECLGen import codon, template
from DECLGen.exceptions import \
    LibraryCategoryException, \
    LibraryElementException, \
    LibraryElementNotFoundException, \
    LibraryElementExistsException
from .elements import Element


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

    def __iter__(self):
        for index in self.elements:
            yield self.elements[index]

    def get_anchors(self):
        return self.anchors

    def set_name(self, name: str) -> None:
        self.name = name

    def set_codon_length(self, codon_length: int = 0) -> None:
        if codon_length < 0:
            raise LibraryCategoryException("Codon length must be a positive integer.")

        self.codon_length = codon_length

    def get_codon_length(self):
        """ Returns the set codon length if set > 0 or the minimum codon length required if set = 0"""
        if self.codon_length > 0:
            return self.codon_length
        else:
            return math.ceil(
                math.log(len(self.elements), len(codon.CodonConfig.bases))
            )

    def clear(self):
        self.elements = {}
        self.element_smiles = []

    def describe(self) -> Dict[str, str]:
        description = {
            "id": self.id,
            "name": self.name,
            "anchors": ", ".join(self.anchors),
            "codon_length": "variable" if self.codon_length == 0 else str(self.codon_length),
            "elements": len(self.elements),
            "max_elements": "unlimited" if self.codon_length == 0 else str(len(codon.CodonConfig.bases)**self.codon_length),
        }

        if self.codon_length == 0:
            description["min_codon_length"] = math.ceil(math.log(len(self.elements), len(codon.CodonConfig.bases)))

        return description

    def has_index(self, index: Union[str, int]):
        if index in self.elements:
            return True
        else:
            return False

    def add_element(self, elm_smiles: str, index: Union[str, int] = None):
        if type(index) == str:
            index = codon.decode(index)
        if index is None:
            if len(self.elements) > 0:
                index = max(list(self.elements.keys()))+1
            else:
                index = 0

        if self.has_index(index):
            raise LibraryElementExistsException(
                "An element with the index <{index}> already exists in category <{cat}>".format(
                    index=index,
                    cat=self.id
                )
            )


        # Check if R-group requirement is met
        anchors = template.get_anchors(elm_smiles)
        not_found_in_elm = []
        num_found_in_elm = 0
        for cat_anchor in self.anchors:
            if cat_anchor not in anchors:
                not_found_in_elm.append(cat_anchor)
            else:
                num_found_in_elm += 1

        if len(not_found_in_elm) > 0:
            raise LibraryElementException("R-Group requirement not met: <{}> not found".format(", ".join(not_found_in_elm)))
        if len(anchors) - num_found_in_elm > 0:
            raise LibraryElementException("R-Group requirement not met: Anchors not belonging to this category found.")

        elm = Element(self, elm_smiles, index)

        # Check if the exact same R-group is already contained
        # This has the flaw that C[R1]CC and CCC[R1] are regarded as being different.
        if elm.smiles in self.element_smiles:
            print("Warning: The exact same R-Group is already contained in this category.")

        self.elements[index] = elm
        self.element_smiles.append(elm.smiles)

    def del_element(self, index: Union[str, int]):
        if type(index) == str:
            index = codon.decode(index)

        if not self.has_index(index):
            raise LibraryElementNotFoundException(
                "An element with the index <{index}> does not exist in category <{cat}>".format(
                    index=index,
                    cat=self.id
                )
            )

        elm = self.elements[index]
        del self.elements[index]
        try:
            self.element_smiles.remove(elm.smiles)
        except ValueError:
            pass
