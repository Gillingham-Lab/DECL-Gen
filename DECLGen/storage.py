from . import template
from .exceptions import LibraryTemplateException
from .library import Library


class Storage:
    """
    This class stores all information from a library and is the class that gets pickled.
    """

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
        """
        Constructor.
        :param advanced_anchors: True if advanced-anchor-mode should be utilized.
        :param superset_categories: True of superset categories should be enabled.
        """
        # Initialize empty values.
        self.raw_template = None
        self.smiles_template = None
        self.dna_template = None

        # Create the (empty) library
        self.library = Library(
            self,
            advanced_anchors=advanced_anchors,
            superset_categories=superset_categories
        )

    def set_template(self, template_string):
        """
        Sets the template.

        This method sets the molecular scaffold if it has not yet been set. It also sets the R-Groups.
        :param template_string:
        :return:
        """
        if self.raw_template is not None:
            raise LibraryTemplateException("Cannot re-set the template via Storage.set_template() after it has been initialized.")

        if template.count_anchors(template_string) <= 0:
            raise LibraryTemplateException("The number of R-Groups ([R1], [R2], etc) must be >= 1")

        self.library.change_template(template_string)

        # Sets the libraries R-Groups. This cannot ever be changed.
        self.library.set_anchors(template.get_anchors(template_string))