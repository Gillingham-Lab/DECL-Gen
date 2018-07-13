import unittest
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from DECLGen.evaluation import codon

class CodonTestCase(unittest.TestCase):
    def test_codon_extraction_from_string(self):
        template = "AAACCCGGGTTT"
        coordinates = [(0,3), (3, 6), (2, 8)]

        result_expected = ["AAA", "CCC", "ACCCGG"]
        result = codon.extract(template, coordinates)
        self.assertListEqual(result_expected, result)

        result_expected = ["CCGGGT", "GGG", "TTT"]
        result = codon.extract(template, coordinates, True)
        self.assertListEqual(result_expected, result)

    def test_codon_extraction_from_bioseq(self):
        template = Seq("AAACCCGGGTTT", IUPAC.unambiguous_dna)
        coordinates = [(0,3), (3, 6), (2, 8)]

        result_expected = ["AAA", "CCC", "ACCCGG"]
        result = codon.extract(template, coordinates)
        self.assertListEqual(result_expected, result)

        result_expected = ["CCGGGT", "GGG", "TTT"]
        result = codon.extract(template, coordinates, True)
        self.assertListEqual(result_expected, result)