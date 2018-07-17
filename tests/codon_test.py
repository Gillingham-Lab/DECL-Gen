import unittest

from DECLGen import codon

class CodonTestCase(unittest.TestCase):
    def test_codon_reversal(self):
        tests = [
            ("AAA", "TTT"),
            ("GGG", "CCC"),
            ("GAG", "CTC"),
            ("AGA", "TCT"),
            ("ACT", "AGT"),
            ("CTG", "CAG"),
            ("ACGT", "ACGT"),
        ]

        for a, b in tests:
            with self.subTest(i=(a, b)):
                self.assertEqual(b, codon.reverse(a))
                self.assertEqual(a, codon.reverse(b))

    def test_codon_encoding(self):
        tests = [
            (0, 3, "AAA"),
            (0, 7, "AAAAAAA"),
            (1, 7, "AAAAAAT"),
            (2, 7, "AAAAAAG"),
            (3, 7, "AAAAAAC"),
            (13, 7, "AAAAACT"),
        ]

        for i, l, result in tests:
            with self.subTest(i=(i, l, result)):
                self.assertEqual(result, codon.encode(i, length=l))

    def test_if_encoding_raises_exception_if_codon_size_not_sufficient(self):
        for i in range(10):
            with self.subTest(i=i):
                with self.assertRaises(ValueError):
                    codon.encode(4**i, i)

    def test_codon_decoding(self):
        tests = [
            ("AAA", 0),
            ("AAAA", 0),
            ("AAAAA", 0),
            ("AAAAAA", 0),
            ("TTT", 21),
            ("TTTT", 85),
            ("GTAC", 147),
            ("CTGG", 218),
            ("AGA", 8),
        ]

        for cod, result in tests:
            with self.subTest(i=(cod, result)):
                self.assertEqual(result, codon.decode(cod))

    def test_if_decode_raises_exception_if_codon_size_contains_unsupported_bases(self):
        tests = [
            "ACTN",
            "UTGA",
            "ABCD",
            "ZXUV"
        ]

        for test in tests:
            with self.subTest(i=test):
                with self.assertRaises(ValueError):
                    codon.decode(test)