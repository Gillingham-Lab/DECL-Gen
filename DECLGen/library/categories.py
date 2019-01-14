import math
import os
from typing import List, Dict, Union, Iterable, Optional
from DECLGen import codon, template
from DECLGen.exceptions import \
    LibraryCategoryException, \
    LibraryCategoryFullException, \
    LibraryCategoryEmptyException, \
    LibraryCategoryAnchorException, \
    LibraryElementException, \
    LibraryElementNotFoundException, \
    LibraryElementExistsException
from .elements import Element


class Category:
    """
        Represents a diversity element category.

        Attributes:
            id: str
                The user-given identifier. Must not contain spaces.
            name: str
                A human-readable name for the category. Can be any text.
            anchors: List[str]
                A list of anchors that this class requires (eg., ["R1", "R2"])
            codon_length: int
                The codon length set for this category. Can be anything if set to 0. Internal use only.
                Consider using get_codon_length() if the minimum codon length is required.
            elements: Dict[str, Element]
                A dictionary that associates encoded codon with the element.
            reverse_complement: Optional[bool]
                True if the codon should be read as it's reverse complement. Internal use only. Consider using is_reverse_complement()
        """
    id: str
    name: str
    anchors: List[str]
    codon_length: int
    elements: Dict[str, Element]
    reverse_complement: Optional[bool] = None

    def __init__(self, id: str, name: str, anchors: List[str], codon_length: int = 0):
        self.id = id
        self.name = name
        self.anchors = anchors
        self.set_codon_length(codon_length)
        self.clear()

    def __iter__(self) -> Iterable[Element]:
        for index in self.elements:
            yield self.elements[index]

    def __len__(self) -> int:
        return len(self.elements)

    def get_anchors(self) -> List[str]:
        """ Returns all anchors required within this category. """
        return self.anchors

    def set_name(self, name: str) -> None:
        """ Sets the descriptive name of this category """
        self.name = name

    def set_codon_length(self, codon_length: int = 0) -> None:
        """ Sets the maximum codon length. Setting to 0 makes the codon flexible."""
        if codon_length < 0:
            raise LibraryCategoryException("Codon length must be a positive integer.")

        self.codon_length = codon_length

    def get_codon_length(self) -> int:
        """ Returns the set codon length if set > 0 or the minimum codon length required if set = 0"""
        if self.codon_length > 0:
            return self.codon_length
        else:
            all_keys = list(self.elements.keys())
            if len(all_keys) == 0:
                return 1

            biggest_key = max(all_keys)
            return max(1, int(math.ceil(
                math.log(max(1, biggest_key+1), len(codon.CodonConfig.bases))
            )))

    def get_max_elements(self) -> int:
        """ Returns the maximum amount of elements available for the current set codon size """
        if self.codon_length == 0:
            return -1
        else:
            return len(codon.CodonConfig.bases)**self.codon_length

    def set_reverse_complement(self, reverse_complement: bool) -> None:
        """ Sets whether this category's codon is reverse complement or not. """
        self.reverse_complement = reverse_complement

    def is_reverse_complement(self) -> bool:
        """ Returns true of this category's codon is to be reverse complemented. """
        if self.reverse_complement is True:
            return True
        else:
            return False

    def clear(self):
        """ Removes all elements from the category. """
        self.elements = {}
        self.element_smiles = []

    def describe(self) -> Dict[str, str]:
        """ Returns a description of the category. """
        description = {
            "id": self.id,
            "name": self.name,
            "anchors": ", ".join(self.anchors),
            "codon_length": "variable" if self.codon_length == 0 else str(self.codon_length),
            "elements": len(self.elements),
            "max_elements": "unlimited" if self.codon_length == 0 else str(self.get_max_elements()),
            "reverse_complement": "Yes" if self.is_reverse_complement() else "No",
        }

        if self.codon_length == 0:
            description["min_codon_length"] = str(self.get_codon_length())

        return description

    def has_index(self, index: Union[str, int]):
        """ Returns true if an index is already in use. """
        if type(index) == str:
            index = codon.decode(index)

        if index in self.elements:
            return True
        else:
            return False

    def get_element(self, index: Union[str, int]) -> Element:
        """ Returns an element with a given index. """
        try:
            index = int(index)
        except ValueError:
            pass

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
        return elm

    def _get_element_by_list_index(self, list_index) -> Element:
        """ Internal helper method to access elements with an enumerated index. """
        elm = list(self.elements.values())[list_index]
        return elm

    def add_element(self, elm_smiles: str, index: Union[str, int] = None):
        """ Adds an element by providing a smiles and optionally an index. If the index is not given,
        DECL-Gen tries to automatically generate one. """
        if type(index) == str:
            old_index = index
            index = codon.decode(index)
        else:
            old_index = None
        if index is None:
            if len(self.elements) > 0:
                index = max(list(self.elements.keys()))+1
            else:
                index = 0

        if self.get_max_elements() > 0 and index >= self.get_max_elements():
            raise LibraryCategoryFullException(
                "This library category <{id}> can only store {max} compounds (codon size={codon}).".format(
                    id=self.id,
                    max=self.get_max_elements(),
                    codon=self.get_codon_length()
                )
            )

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
        if elm.raw_smiles in self.element_smiles:
            print("Warning: The exact same element is already contained in this category ({id}, {dna}, {smiles}).".format(
                id=index,
                dna=old_index,
                smiles=elm_smiles
            ))
            for k in self.elements:
                if self.elements[k].raw_smiles == elm.raw_smiles:
                    print("Other item at ", self.elements[k].index, codon.encode(self.elements[k].index, self.codon_length))

        self.elements[index] = elm
        self.element_smiles.append(elm.raw_smiles)

    def del_element(self, index: Union[str, int]):
        """ Deletes an element from a diversity element category"""
        elm = self.get_element(index)

        del self.elements[elm.index]
        try:
            self.element_smiles.remove(elm.raw_smiles)
            del elm
        except ValueError:
            pass

    def import_elements(self, filename: str, updateable = None) -> int:
        """ Imports elements from a tab separated file. """
        if not os.path.exists(filename):
            raise LibraryElementException("The file <{}> does not exist.".format(filename))

        with open(filename, "r") as file:
            # Count lines and auto-detect format
            line_count = 0
            for line in file:
                line_count += 1

            file.seek(0)

            c = 0
            for line in file:
                cols = line.strip().split("\t")
                codon = cols[0]
                smiles = cols[1]

                if self.has_index(codon):
                    self.del_element(codon)

                self.add_element(smiles, codon)

                updateable.update(c/line_count)
                c += 1

        return c

    def import_elements_from_cat(self, origin: 'Category', anchors: Dict[str, str], updateable=None) -> int:
        if len(origin) == 0:
            raise LibraryCategoryEmptyException("The origin category is empty, cannot import from this category.")

        anchor_origin = [*anchors.keys()]
        anchor_target = [*anchors.values()]

        # Check anchor consistency in origin
        found = 0
        for anchor in anchor_origin:
            if anchor in origin.anchors:
                found += 1

        if found < len(origin.anchors):
            raise LibraryCategoryAnchorException("Not all anchors set for category {} are given ({}).".format(origin.id, ",".join(origin.get_anchors())))

        # Check anchor consistency in target (self)
        found = 0
        for anchor in anchor_target:
            if anchor in self.get_anchors():
                found += 1

        if found < len(self.anchors):
            raise LibraryCategoryAnchorException("Not all anchors set for category {} are given ({}).".format(self.id, ",".join(self.get_anchors())))

        # Clear self elements
        self.clear()

        added = 0
        for element in origin:
            index = element.index
            smiles = element.raw_smiles

            # Translate table
            for k in range(len(anchors)):
                search = "[{}]".format(anchor_origin[k])
                replce = "[{}]".format(anchor_target[k])

                smiles = smiles.replace(search, replce)

            self.add_element(smiles, index)
            added += 1

        return added






