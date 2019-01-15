from . import template
from .exceptions import LibraryTemplateException
from .library import Library


class Storage:
    # This sequence defines the number-value of each base
    bases = "ATGC"
    # The raw template as given by the user, sanitized
    raw_template = None
    # The smiles-converted raw template
    smiles_template = None
    # The raw dna string as given by the user
    dna_template = None
    # The library
    library = None

    def __init__(self, advanced_anchors: bool, superset_categories: bool):
        self.raw_template = None
        self.smiles_template = None
        self.dna_template = None

        self.library = Library(
            self,
            advanced_anchors=advanced_anchors,
            superset_categories=superset_categories
        )

    def set_template(self, template_string):
        if self.raw_template is not None:
            raise LibraryTemplateException("Cannot re-set the template via Storage.set_template() after it has been initialized.")

        if template.count_anchors(template_string) <= 0:
            raise LibraryTemplateException("The number of R-Groups ([R1], [R2], etc) must be >= 1")

        self.library.change_template(template_string)

        # Sets the libraries R-Groups
        self.library.set_anchors(template.get_anchors(template_string))