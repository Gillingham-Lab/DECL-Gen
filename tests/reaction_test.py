import unittest
from rdkit import Chem, RDLogger

from DECLGen.molecule import Molecule
from DECLGen.reactions import Reaction, predefinedReactions
from DECLGen.rtools import unmask

class ReactionTestCase(unittest.TestCase):
    def test_if_custom_reaction_gives_product(self):
        m = Molecule("CCO")
        rxn = Reaction("[C:1]-[O:2]>>[C:1]-[O:2]-C(C)(C)C")

        rxn.react(m)

        new_smiles = Chem.MolToSmiles(m._mol)

        self.assertNotEqual("CCO", new_smiles)
        self.assertEqual("CCOC(C)(C)C", new_smiles)

    def test_if_AmideCoupling_Amine_gives_correct_products(self):
        tests = [
            # Primary amines
            ("CN", "CN[R1]"),
            # secondary amines
            ("CCNCC", "CCN([R1])CC"),
            # Amines, not Amides
            ("CC(=O)NCCN", "CC(=O)NCCN[R1]"),
            # Amines, not anilines
            ("Nc1ccc(CN)cc1", "Nc1ccc(CN[R1])cc1"),
        ]

        for i in range(len(tests)):
            with self.subTest(i=i):
                input, expected = tests[i]

                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["AmideCoupling_Amine"], rGroup="R1")

                rxn.react(m)
                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertNotEqual(input, result)
                self.assertEqual(expected, result)