from typing import List, Set, Dict, Generator, Tuple
from operator import mul, itemgetter
from functools import reduce
from DECLGen.exceptions import \
    LibraryCategoryException, \
    LibraryCategoryExistsException, \
    LibraryCategoryNotFoundException, \
    LibraryTemplateException
from .categories import Category
from DECLGen import template


def _check_anchor(anchor) -> bool:
    if not anchor.startswith("R"):
        raise LibraryCategoryException(
            "R-Group <{}> must start with R.".format(anchor)
        )

    if len(anchor) < 2:
        raise LibraryCategoryException(
            "R-Group <{}> must be Rn or Rnn, whereas n is a number from 0-9.".format(anchor)
        )

    if int(anchor[1:]) < 0:
        raise LibraryCategoryException(
            "R-Group <{}> must be a positive number.".format(anchor)
        )

    return True


class Library:
    storage = None
    anchors = None
    anchors_in_use = None
    categories = None

    def __init__(self, storage):
        self.storage = storage
        self.anchors = []
        self.anchors_in_use = {}
        self.categories = {}

    def set_anchors(self, anchors: List[str]):
        self.anchors = sorted(anchors)

    def change_template(self, new_template: str):
        # We must sanitize the template first in case the R group is at the beginning
        template_string = template.sanitize(new_template)

        raw_template = template_string
        smiles_template = template.parse(template_string)
        anchors = template.get_anchors(template_string)

        if len(self.anchors) > 0 and sorted(anchors) != sorted(self.anchors):
            raise LibraryTemplateException("Cannot change template with different anchors. You must remove the library first.")

        self.storage.raw_template = raw_template
        self.storage.smiles_template = smiles_template

    def describe(self) -> Dict[str, str]:
        description = {
            "Template": self.storage.raw_template,
            "R-Groups": ", ".join(self.anchors),
            "Library Shape": ", ".join([str(x) for x in self.shape]),
            "Library Size": self.get_size()
        }

        return description

    def get_size(self) -> int:
        return reduce(mul, self.shape, 1)

    def has_category(self, id: str):
        if id in self.categories:
            return True
        else:
            return False

    def get_categories(self) -> List[Category]:
        return list(self.categories.values())

    def get_category(self, id) -> Category:
        if not self.has_category(id):
            raise LibraryCategoryNotFoundException("Library category <{}> not found".format(id))

        return self.categories[id]

    @property
    def shape(self) -> Tuple[int]:
        shape = []
        for cat in self.categories.values():
            shape.append(len(cat))
        return tuple(shape)

    def add_category(self, id: str, name: str, anchors: List[str], codon_length: int = 0) -> bool:
        """
        Adds a new diversity element category to the library.
        :param id: The id of the category
        :param name: The name of the category
        :param anchors: A list of R-Group anchors (R1, R2 ...)
        :param codon_length: Codon length
        :return:
        """
        # Check if id is already in use
        if self.has_category(id):
            raise LibraryCategoryExistsException("A category with the id <{}> already exists".format(id))

        # Check validity of each anchor and if anchor is free
        if len(anchors) == 0:
            raise LibraryCategoryException("You must at least give 1 anchor, none were given.")

        anchors_used = []
        for anchor in anchors:
            _check_anchor(anchor)

            if anchor in self.anchors_in_use:
                anchors_used.append(anchor)

        if len(anchors_used) > 0:
            raise LibraryCategoryException(
                "Anchors <{}> are already in use.".format(", ".join(anchors_used))
            )

        # Create the new category
        self.categories[id] = Category(id, name, anchors, codon_length)
        # Marks anchors as used by cross-referencing Category
        for anchor in anchors:
            self.anchors_in_use[anchor] = self.categories[id]

        return True

    def del_category(self, id: str) -> bool:
        # Check if id is already in use
        if not self.has_category(id):
            raise LibraryCategoryNotFoundException("A category with the id <{}> does not exist exists".format(id))

        # Remove used anchors
        anchors = self.categories[id].get_anchors()
        for anchor in anchors:
            del self.anchors_in_use[anchor]

        # Remove category
        del self.categories[id]

        return True

    def get_molecule_data_by_index(self, elements: Dict[str, int]) -> str:
        fragments = []
        codons = []
        for cat_id in elements:
            fragments.append(self.categories[cat_id].get_element_by_index(elements[cat_id]).parsed_smiles)
            codons.append(self.categories[cat_id].get_element_by_index(elements[cat_id]).codon)

        fragments.append(self.storage.smiles_template)
        smiles = ".".join(fragments)

        return smiles, codons


    def generate_molecule_queue(self) -> Tuple[list, list]:
        categories_by_size = []
        for cat_id in self.categories:
            categories_by_size.append((len(self.categories[cat_id]), cat_id))

        categories_by_size.sort()

        queue = []
        for element1 in range(categories_by_size[0][0]):
            gen = ([categories_by_size[0][1], element1], categories_by_size[1:])
            queue.append(gen)

        return queue, [x[1] for x in categories_by_size]