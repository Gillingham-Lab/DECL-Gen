import math
import os
from typing import List, Dict, Union, Iterable, Optional, Tuple
from DECLGen import codon, template
from DECLGen.exceptions import \
    LibraryCategoryException, \
    LibraryCategoryFullException, \
    LibraryCategoryEmptyException, \
    LibraryCategoryAnchorException, \
    LibraryCategoryExistsException, \
    LibraryCategoryNotFoundException, \
    LibraryElementException, \
    LibraryElementNotFoundException, \
    LibraryElementExistsException
from .elements import Element

def _get_new_index(objects: Dict[int, object]) -> codon.CodonInt:
    if len(objects) > 0:
        index = max(list(objects.keys())) + 1
    else:
        index = 0
    return index


class BaseCategory:
    @staticmethod
    def id_is_valid(id):
        if id.isalnum() and id[0].isalpha():
            return True
        return False

    @staticmethod
    def check_id_validity(id):
        if BaseCategory.id_is_valid(id) is False:
            raise LibraryCategoryException(
                "An id must start with a letter (a-z, A-Z). Any other position must be a letter or a number (0-9).")

    def __init__(self):
        pass

    def set_name(self, name: str) -> None:
        """ Sets the descriptive name of this category """
        self.name = name

    def get_anchors(self) -> List[str]:
        """ Returns all anchors required within this category. """
        return self.anchors

    def set_codon_length(self, codon_length: int = 0) -> None:
        """ Sets the maximum codon length. Setting to 0 makes the codon flexible."""
        if codon_length < 0:
            raise LibraryCategoryException("Codon length must be a positive integer.")

        self.codon_length = codon_length
        return

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

    def _clear(self):
        """ Removes all elements from the category. """
        self.elements = {}
        self.element_smiles = []

    def get_max_elements(self) -> int:
        """ Returns the maximum amount of elements available for the current set codon size """
        if self.codon_length == 0:
            return -1
        else:
            return len(codon.CodonConfig.bases)**self.codon_length

    def _raise_exception_if_full(self, index, what, of):
        if self.get_max_elements() > 0 and index >= self.get_max_elements():
            raise LibraryCategoryFullException(
                "This {what} <{id}> can only store {max} {of} (codon size={codon}).".format(
                    id=self.id,
                    max=self.get_max_elements(),
                    codon=self.get_codon_length(),
                    what=what,
                    of=of,
                )
            )

    def set_reverse_complement(self, reverse_complement: bool) -> None:
        self.reverse_complement = reverse_complement

    def is_reverse_complement(self) -> bool:
        """ Returns true of this category's codon is to be reverse complemented. """
        if self.reverse_complement is True:
            return True
        else:
            return False

    def is_subset(self) -> bool:
        return False

    def is_superset(self) -> bool:
        return False


class Category(BaseCategory):
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
        super().__init__()

        self.id = id
        self.name = name
        self.anchors = anchors
        self.set_codon_length(codon_length)
        self._clear()

    def __iter__(self) -> Iterable[Element]:
        for index in self.elements:
            yield self.elements[index]

    def __len__(self) -> int:
        return len(self.elements)

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
            "extra": "",
        }

        if self.codon_length == 0:
            description["min_codon_length"] = str(self.get_codon_length())

        return description

    def has_index(self, index: codon.CodonType):
        """ Returns true if an index is already in use. """
        index = codon.normalize(index)

        if index in self.elements:
            return True
        else:
            return False

    def clear(self):
        self.__clear()

    def get_element(self, index: codon.CodonType) -> Element:
        """ Returns an element with a given index. """
        index = codon.normalize(index)

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

    def add_element(self, elm_smiles: str, index: Optional[codon.CodonType] = None):
        """ Adds an element by providing a smiles and optionally an index. If the index is not given,
        DECL-Gen tries to automatically generate one. """
        if index is None:
            index = _get_new_index(self.elements)
        else:
            index = codon.normalize(index)

        self._raise_exception_if_full(index, "library category", "elements")

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
        # @ToDo: Change to add methyl group to all R positions; maybe save methylated canonical smiles, too.
        if elm.raw_smiles in self.element_smiles:
            print("Warning: The exact same element is already contained in this category ({id}, {dna}, {smiles}).".format(
                id=index,
                dna=index,
                smiles=elm_smiles
            ))
            for k in self.elements:
                if self.elements[k].raw_smiles == elm.raw_smiles:
                    print("Other item at ", self.elements[k].index, codon.encode(self.elements[k].index, self.codon_length))

        self.elements[index] = elm
        self.element_smiles.append(elm.raw_smiles)

    def del_element(self, index: codon.CodonType):
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

    def copy_elements_from(self, origin: 'Category', anchors: Dict[str, str], updateable=None) -> int:
        if len(origin) == 0:
            raise LibraryCategoryEmptyException("The origin category is empty, cannot import from this category.")

        if isinstance(origin, 'NoElementCategory'):
            raise LibraryCategoryException("You cannot import elements from the given category.")

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
                replce = "[::{}::]".format(anchor_target[k])

                smiles = smiles.replace(search, replce)

            smiles = smiles.replace("[::", "[")
            smiles = smiles.replace("::]", "]")

            self.add_element(smiles, index)
            added += 1

        return added


class NoElementCategory:
    def add_element(self, elm_smiles: str, index: Union[str, int] = None):
        raise LibraryCategoryException("You cannot add an element to this category.")

    def del_element(self, index: Union[str, int]):
        raise LibraryCategoryException("You cannot remove an element from this category.")

    def import_elements(self, filename: str, updateable=None) -> int:
        raise LibraryCategoryException("You cannot import elements into this category.")

    def copy_elements_from(self, origin: 'Category', anchors: Dict[str, str], updateable=None) -> int:
        raise LibraryCategoryException("You cannot import elements into this category.")


class SubCategory(Category):
    def codon(self, superset):
        return codon.encode(self.id, superset.get_codon_length())


class SupersetCategory(BaseCategory, NoElementCategory):
    id: str
    name: str
    anchors: List[str]
    codon_length: int
    elements: Dict[str, Element]
    reverse_complement: Optional[bool]
    subset: "SubsetCategory"
    categories: Dict[int, str]

    def __init__(self, id1: str, id2: str, name: str, anchors: List[str], codon1_length: int = 0, codon2_length: int = 0):
        super().__init__()

        self.id = id1
        self.name = name
        self.anchors = anchors
        self.set_codon_length(codon1_length)
        self._clear()
        self.clear_cats()

        self.subset = SubsetCategory(self, id2, codon2_length)

    def __iter__(self) -> Iterable[Element]:
        for index in self.categories:
            yield self.categories[index]

    def __len__(self) -> int:
        sum = 0
        for subcat in self.categories.values():
            sum += len(subcat)
        return sum

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
            "extra": "Superset-category (bound with {})".format(self.subset.id)
        }

        if self.codon_length == 0:
            description["min_codon_length"] = str(self.get_codon_length())

        return description

    def get_subset_category(self) -> "SubsetCategory":
        return self.subset

    def is_superset(self) -> bool:
        return True

    def clear(self):
        raise LibraryCategoryException(
            "You cannot clear a superset category with cat-clear. " + \
            "Instead, use ssc-clear-cats to remove all categories or ssc-clear-elms to remove all elements from a specific subcategory."
        )

    def clear_cats(self):
        self.categories = {}

    def has_index(self, index: codon.CodonInt):
        if index in self.categories:
            return True
        else:
            return False

    def get_next_index(self) -> codon.CodonInt:
        return _get_new_index(self.categories)

    def add_category(self, name, index: Optional[codon.CodonType]) -> codon.CodonInt:
        if index is None:
            index = self.get_next_index()
        else:
            index = codon.normalize(index)

        self._raise_exception_if_full(index, "superset category", "categories")

        if self.has_index(index):
            raise LibraryCategoryExistsException(
                "A category with the index <{index}> already exists in the superset category <{cat}>".format(
                    index=index,
                    cat=self.id
                )
            )

        # Create category
        cat = SubCategory(index, name, anchors=self.anchors, codon_length=self.subset.codon_length)
        cat.set_reverse_complement(self.subset.reverse_complement)

        # Add category to list.
        self.categories[index] = cat

        return index

    def del_category(self, index: codon.CodonType):
        index = codon.normalize(index)

        del self.categories[index]

    def get_category(self, index: codon.CodonType):
        index = codon.normalize(index)

        if self.has_index(index) is False:
            raise LibraryCategoryNotFoundException(
                "A category with the index <{index}> does not exist in the superset category <{cat}>".format(
                    index=index,
                    cat=self.id
                )
            )

        return self.categories[index]

    def _get_element_by_list_index(self, list_index) -> Tuple[Category, Element]:
        """ Internal helper method to access elements with an enumerated index. """
        i = list_index

        for subcat in self:
            if i >= len(subcat):
                i -= len(subcat)
            else:
                elm = subcat._get_element_by_list_index(i)
                break

        return subcat, elm



class SubsetCategory(BaseCategory, NoElementCategory):
    id: str
    codon_length: int
    elements: Dict[str, Element]
    reverse_complement: Optional[bool] = None
    superset: "SupersetCategory" = None

    def __init__(self, superset: SupersetCategory, id, codon_length):
        super().__init__()

        self.superset = superset
        self.id = id
        self.set_codon_length(codon_length)
        self._clear()

    def __len__(self) -> int:
        return 1

    def set_name(self, name: str) -> None:
        """ Sets the descriptive name of this category """
        raise LibraryCategoryException("You cannot change the name of a subset category. Please change the name of the corresponding superset.")

    def describe(self) -> Dict[str, str]:
        """ Returns a description of the category. """
        description = {
            "id": self.id,
            "name": self.superset.name,
            "anchors": ", ".join(self.superset.anchors),
            "codon_length": "variable" if self.codon_length == 0 else str(self.codon_length),
            "elements": len(self.elements),
            "max_elements": "unlimited" if self.codon_length == 0 else str(self.get_max_elements()),
            "reverse_complement": "Yes" if self.is_reverse_complement() else "No",
            "extra": "Subset-category (bound with {})".format(self.superset.id)
        }

        if self.codon_length == 0:
            description["min_codon_length"] = str(self.get_codon_length())

        return description

    def get_superset_category(self) -> SupersetCategory:
        return self.superset

    def is_subset(self) -> bool:
        return True

    def clear(self):
        raise LibraryCategoryException("You cannot clear a subset category.")

