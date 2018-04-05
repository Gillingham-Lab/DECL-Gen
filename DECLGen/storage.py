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

    def __init__(self):
        self.raw_template = None
        self.smiles_template = None
        self.dna_template = None

        self.library = Library(self)

    def set_template(self, template_string):
        if template.count_anchors(template_string) <= 0:
            raise LibraryTemplateException("The number of R-Groups ([R1], [R2], etc) must be >= 1")

        # We must sanitize the template first in case the R group is at the beginning
        template_string = template.sanitize(template_string)

        self.raw_template = template_string
        self.smiles_template = template.parse(template_string)

        # Sets the libraries R-Groups
        self.library.set_anchors(template.get_anchors(template_string))