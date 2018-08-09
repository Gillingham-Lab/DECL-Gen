import unittest

from DECLGen import template

class TemplateTestCase(unittest.TestCase):
    def test_count_anchors(self):
        tests = [
            ("CCC([R1])C(O)=O", 1),
            ("CCC([A])C(O)=O", 0),
            ("CCC([R7A])C(O)=O", 0),
            ("CCC([R88])C(O)=O", 1),
            ("CCC(N)C([R2])=O", 1),
            ("CCC([R1])C([R2])=O", 2),
            ("CCC(N[R2])([R3])C([R1])=O", 3),
        ]

        for a, b in tests:
            with self.subTest(i=(a, b)):
                self.assertEqual(b, template.count_anchors(a))

    def test_get_anchors(self):
        tests = [
            ("CCC([R1])C(O)=O", ["R1"]),
            ("CCC([A])C(O)=O", []),
            ("CCC([R7A])C(O)=O", []),
            ("CCC([R88])C(O)=O", ["R88"]),
            ("CCC(N)C([R2])=O", ["R2"]),
            ("CCC([R1])C([R2])=O", ["R1", "R2"]),
            ("CCC([R2])C([R1])=O", ["R1", "R2"]),
            ("CCC(N[R2])([R3])C([R1])=O", ["R1", "R2", "R3"]),
        ]

        for a, b in tests:
            with self.subTest(i=(a, b)):
                self.assertListEqual(b, template.get_anchors(a))

    def test_scaffold_sanitation(self):
        tests = [
            ("CCC([R1])C(O)=O", "CCC[R1]C(O)=O"),
            ("CCC([A])C(O)=O", "CCC([A])C(O)=O"),
            ("CCC([R7A])C(O)=O", "CCC([R7A])C(O)=O"),
            ("CCC([R88])C(O)=O", "CCC[R88]C(O)=O"),
            ("CCC(N)C([R2])=O", "CCC(N)C[R2]=O"),
            ("CCC([R1])C([R2])=O", "CCC[R1]C[R2]=O"),
            ("CCC([R2])C([R1])=O", "CCC[R2]C[R1]=O"),
            ("CCC(N[R2])C([R1])=O", "CCC(N[R2])C[R1]=O"),
            ("CCC(N[R2])([R3])C([R1])=O", "CCC[R3](N[R2])C[R1]=O"),
        ]

        for a, b in tests:
            with self.subTest(i=(a, b)):
                self.assertEqual(b, template.sanitize(a))

    def test_if_parsing_correct_scaffolds_returns_correct_string(self):
        tests = [
            ("CCC([R1])C(O)=O", "CCC%98C(O)=O"),
            ("CCC([A])C(O)=O", "CCC([A])C(O)=O"),
            ("CCC([R7A])C(O)=O", "CCC([R7A])C(O)=O"),
            ("CCC([R88])C(O)=O", "CCC%11C(O)=O"),
            ("CCC(N)C([R2])=O", "CCC(N)C%97=O"),
            ("CCC([R1])C([R2])=O", "CCC%98C%97=O"),
            ("CCC([R2])C([R1])=O", "CCC%97C%98=O"),
            ("CCC(N[R2])C([R1])=O", "CCC(N%97)C%98=O"),
            ("CCC(N[R2])([R3])C([R1])=O", "CCC%96(N%97)C%98=O"),
        ]

        for a, b in tests:
            with self.subTest(i=(a, b)):
                self.assertEqual(b, template.parse(template.sanitize(a)))

    def test_if_parsing_scaffolds_with_incorrect_anchors_raises_exception(self):
        tests = [
            ("CCC([R90])C(O)=O"),
        ]

        for a in tests:
            with self.subTest(i=(a)):
                with self.assertRaises(ValueError):
                    b = template.parse(template.sanitize(a))

    def test_if_is_valid_returns_are_true_if_smiles_is_good(self):
        tests = [
            "CCC([R1])C(O)=O",
            "CCC([R88])C(O)=O",
            "CCC(N)C([R2])=O",
            "CCC([R1])C([R2])=O",
            "CCC(N[R2])([R3])C([R1])=O",
        ]

        for a in tests:
            with self.subTest(i=(a)):
                self.assertTrue(template.is_valid(template.parse(template.sanitize(a)), template.get_anchors(a)))

    def test_if_is_valid_returns_are_false_if_smiles_is_not_good(self):
        tests = [
            "C([A])CC([R1])C(O)=O",
            "CCC([U7])C(O)=O",
        ]

        for a in tests:
            with self.subTest(i=(a)):
                self.assertFalse(template.is_valid(template.parse(template.sanitize(a)), template.get_anchors(a)))