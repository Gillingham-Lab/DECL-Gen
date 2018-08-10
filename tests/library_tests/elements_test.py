import unittest
from rdkit import RDLogger

from DECLGen.library import Category, Element
from DECLGen.exceptions import LibraryElementException

class ElementsTestCase(unittest.TestCase):
    def _get_dummy_cat(self):
        cat = Category("dummy", "Dummy", ["R1"], 3)
        return cat

    def test_if_element_creation_works_with_valid_parameters(self):
        cat = self._get_dummy_cat()
        element = Element(cat, "CCO[R1]", 0)
        self.assertIsNotNone(element)

    def test_if_exception_is_raised_if_smiles_is_invalid(self):
        # Silence RDKit logger
        logger = RDLogger.logger()
        logger.setLevel(RDLogger.CRITICAL)

        cat = self._get_dummy_cat()

        with self.assertRaises(LibraryElementException):
            element = Element(cat, "C(C)(C)(C)(C)(C)", 0)

    def test_if_exception_is_raised_if_category_is_not_a_category(self):
        with self.assertRaises(ValueError):
            element = Element(None, "CCO[R1]", 0)

    def test_if_exception_is_raised_if_index_is_negative(self):
        cat = self._get_dummy_cat()

        with self.assertRaises(ValueError):
            element = Element(cat, "CCO[R1]", -1)

    def test_if_exception_is_raised_if_index_is_not_integer(self):
        cat = self._get_dummy_cat()

        with self.assertRaises(ValueError):
            element = Element(cat, "CCO[R1]", "Hallo")


    def test_if_element_returns_codon(self):
        cat = self._get_dummy_cat()
        element = Element(cat, "CCO[R1]", 0)

        self.assertEqual(element.codon, "AAA")
