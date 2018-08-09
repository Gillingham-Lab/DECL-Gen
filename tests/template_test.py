import unittest

from DECLGen import template

class TemplateTestCase(unittest.TestCase):
    def test_count_anchors(self):
        tests = [
            ("CCC([R1])C(O)=O", 1),
            ("CCC([A])C(O)=O", 0),
            ("CCC([R7A])C(O)=O", 0),
            ("CCC([R88])C(O)=O", 1),
            ("CCC(N)C([R2])=0", 1),
            ("CCC([R1])C([R2])=0", 2),
            ("CCC(N[R2])([R3])C([R1])=0", 3),
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
            ("CCC(N)C([R2])=0", ["R2"]),
            ("CCC([R1])C([R2])=0", ["R1", "R2"]),
            ("CCC([R2])C([R1])=0", ["R1", "R2"]),
            ("CCC(N[R2])([R3])C([R1])=0", ["R1", "R2", "R3"]),
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
            ("CCC(N)C([R2])=0", "CCC(N)C[R2]=0"),
            ("CCC([R1])C([R2])=0", "CCC[R1]C[R2]=0"),
            ("CCC([R2])C([R1])=0", "CCC[R2]C[R1]=0"),
            ("CCC(N[R2])C([R1])=0", "CCC(N[R2])C[R1]=0"),
            ("CCC(N[R2])([R3])C([R1])=0", "CCC[R3](N[R2])C[R1]=0"),
        ]

        for a, b in tests:
            with self.subTest(i=(a, b)):
                self.assertEqual(b, template.sanitize(a))

    def test_if_parsing_correct_scaffolds_returns_correct_string(self):
        tests = [
            ("CCC([R1])C(O)=O", "CCC%99C(O)=O"),
            ("CCC([A])C(O)=O", "CCC([A])C(O)=O"),
            ("CCC([R7A])C(O)=O", "CCC([R7A])C(O)=O"),
            ("CCC([R88])C(O)=O", "CCC%11C(O)=O"),
            ("CCC(N)C([R2])=0", "CCC(N)C%98=0"),
            ("CCC([R1])C([R2])=0", "CCC%99C%98=0"),
            ("CCC([R2])C([R1])=0", "CCC%98C%99=0"),
            ("CCC(N[R2])C([R1])=0", "CCC(N%98)C%99=0"),
            ("CCC(N[R2])([R3])C([R1])=0", "CCC%97(N%98)C%99=0"),
        ]