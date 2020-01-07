import unittest
from rdkit import Chem, RDLogger

from DECLGen.molecule import Molecule
from DECLGen.reactions import Reaction, predefinedReactions
from DECLGen.rtools import unmask
from DECLGen.exceptions import ReactionNoProductException

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

    def test_if_AmideCoupling_Amine_is_specific_enough(self):
        tests = [
            # Prevent Amide
            ("CCNC(C)=O"),
            # Prevent Anilines
            ("Cc1ccc(N)cc1"),
            # Prevent protected amines
            ("CCNC(=O)OC(C)(C)C")
        ]

        for smiles in tests:
            with self.subTest(smiles=smiles):
                input = smiles

                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["AmideCoupling_Amine"], rGroup="R1")

                with self.assertRaises(ReactionNoProductException):
                    rxn.react(m)

                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertEqual(input, result)

    def test_if_AmideCoupling_Acids_gives_correct_products(self):
        tests = [
            # Formic acid
            ("O=CO", "O=C[R1]"),
            # Acetic Acid
            ("CC(=O)O", "CC(=O)[R1]"),
            # Propionic Acid
            ("CCC(=O)O", "CCC(=O)[R1]"),
            # 2-Methylpropionic acid
            ("CC(C)C(=O)O", "CC(C)C(=O)[R1]"),
            # 2,2-Dimethylpropionic acid
            ("CC(C)(C)C(=O)O", "CC(C)(C)C(=O)[R1]"),
            # Benzoic acid
            ("c1cc(C(=O)O)ccc1", "O=C([R1])c1ccccc1"),
        ]

        for i in range(len(tests)):
            with self.subTest(i=i):
                input, expected = tests[i]

                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["AmideCoupling_Acid"], rGroup="R1")

                rxn.react(m)
                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertNotEqual(input, result)
                self.assertEqual(expected, result)

    def test_if_AmideCoupling_Acid_is_specific_enough(self):
        tests = [
            # Prevent Methyl ester
            ("CCC(=O)OC"),
            # Prevent Aminohydroxyl(?)
            ("CCC(=O)ON"),
            # Prevent Anhydrides
            ("CCC(=O)OC(C)=O")
        ]

        for smiles in tests:
            with self.subTest(smiles=smiles):
                input = smiles

                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["AmideCoupling_Amine"], rGroup="R1")

                with self.assertRaises(ReactionNoProductException):
                    rxn.react(m)

                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertEqual(input, result)

    def test_if_CuAAC_Alkyne_gives_correct_products(self):
        tests = [
            # Propyne
            ("CC#C", "Cc1cn([R1])nn1"),
            # Phenylethyne
            ("c1ccccc1C#C", "[R1]n1cc(-c2ccccc2)nn1"),
        ]

        for i in range(len(tests)):
            with self.subTest(i=i):
                input, expected = tests[i]

                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["CuAAC_Alkyne"], rGroup="R1")

                rxn.react(m)
                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertNotEqual(input, result)
                self.assertEqual(expected, result)

    def test_if_CuAAC_Alkyne_is_specific_enough(self):
        tests = [
            # Prevent non-terminal alkynes
            ("CC#CC"),
        ]

        for smiles in tests:
            with self.subTest(smiles=smiles):
                input = smiles

                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["CuAAC_Alkyne"], rGroup="R1")

                with self.assertRaises(ReactionNoProductException):
                    rxn.react(m)

                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertEqual(input, result)

    def test_if_CuAAC_Azide_gives_correct_products(self):
        tests = [
            # 3-Azidopentane
            ("CCC(N=[N+]=[N-])CC", "CCC([R1])CC"),
            # Phenylazide
            ("c1ccccc1(N=[N+]=[N-])", "[R1]c1ccccc1"),
        ]

        for i in range(len(tests)):
            with self.subTest(i=i):
                input, expected = tests[i]

                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["CuAAC_Azide"], rGroup="R1")

                rxn.react(m)
                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertNotEqual(input, result)
                self.assertEqual(expected, result)

    def test_if_CuAAC_Azide_is_specific_enough(self):
        tests = [
            # Prevent triazenes
            ("CN=NNC")
        ]

        for smiles in tests:
            with self.subTest(smiles=smiles):
                input = smiles

                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["CuAAC_Azide"], rGroup="R1")

                with self.assertRaises(ReactionNoProductException):
                    rxn.react(m)

                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertEqual(input, result)

    def test_if_NAr_Phenol_gives_correct_products(self):
        tests = [
            # Phenol
            ("Oc1ccccc1", "[R1]Oc1ccccc1"),
        ]

        for i in range(len(tests)):
            with self.subTest(i=i):
                input, expected = tests[i]

                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["NAr_Phenol"], rGroup="R1")

                rxn.react(m)
                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertNotEqual(input, result)
                self.assertEqual(expected, result)

    def test_if_NAr_Phenol_is_specific_enough(self):
        tests = [
            # Prevent triazenes
            ("CCO")
        ]

        for smiles in tests:
            with self.subTest(smiles=smiles):
                input = smiles

                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["NAr_Phenol"], rGroup="R1")

                with self.assertRaises(ReactionNoProductException):
                    rxn.react(m)

                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertEqual(input, result)

    def test_if_NAr_Aniline_gives_correct_products(self):
        tests = [
            # Aniline
            ("Nc1ccccc1", "[R1]Nc1ccccc1"),
        ]

        for i in range(len(tests)):
            with self.subTest(i=i):
                input, expected = tests[i]

                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["NAr_Aniline"], rGroup="R1")

                rxn.react(m)
                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertNotEqual(input, result)
                self.assertEqual(expected, result)

    def test_if_NAr_Aniline_is_specific_enough(self):
        tests = [
            # Prevent amines
            ("CCN"),
            # Prevent protected anilines
            ("O=CNc1ccccc1"),
        ]

        for smiles in tests:
            with self.subTest(smiles=smiles):
                input = smiles

                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["NAr_Aniline"], rGroup="R1")

                with self.assertRaises(ReactionNoProductException):
                    rxn.react(m)

                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertEqual(input, result)

    def test_if_NAr_Amine_gives_correct_products(self):
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
                rxn = Reaction(**predefinedReactions["NAr_Amine"], rGroup="R1")

                rxn.react(m)
                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertNotEqual(input, result)
                self.assertEqual(expected, result)

    def test_if_NAr_Amine_is_specific_enough(self):
        tests = [
            # Prevent Amide
            ("CCNC(C)=O"),
            # Prevent Anilines
            ("Cc1ccc(N)cc1"),
            # Prevent protected amines
            ("CCNC(=O)OC(C)(C)C")
        ]

        for smiles in tests:
            with self.subTest(smiles=smiles):
                input = smiles

                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["NAr_Amine"], rGroup="R1")

                with self.assertRaises(ReactionNoProductException):
                    rxn.react(m)

                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertEqual(input, result)

    def test_if_NAr_Acceptor_gives_correct_products(self):
        tests = [
            # Fluorbenzene
            ("Fc1ccccc1", "[R1]c1ccccc1"),
            # 2-Chloropyridine
            ("Clc1ncccc1", "[R1]c1ccccn1"),
            # 4-Chloropyridine
            ("Clc1ccncc1", "[R1]c1ccncc1"),
            # 2-Chloro-2-nitrobenzene
            ("Clc1c([N+](=O)[O-])cccc1", "O=[N+]([O-])c1ccccc1[R1]"),
            # 4-Chloro-2-nitrobenzene
            ("Clc1ccc([N+](=O)[O-])cc1", "O=[N+]([O-])c1ccc([R1])cc1"),
        ]

        for input, expected in tests:
            with self.subTest(input=input, expected=expected):
                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["NAr_Acceptor"], rGroup="R1")

                rxn.react(m)
                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertNotEqual(input, result)
                self.assertEqual(expected, result)

    def test_if_NAr_Acceptor_is_specific_enough(self):
        tests = [
            # Prevent Chlorbenzene
            ("Clc1ccccc1"),
        ]

        for smiles in tests:
            with self.subTest(smiles=smiles):
                input = smiles

                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["NAr_Acceptor"], rGroup="R1")

                with self.assertRaises(ReactionNoProductException):
                    rxn.react(m)

                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertEqual(input, result)

    def test_if_Suzuki_Halogen_gives_correct_products(self):
        tests = [
            # Fluorbenzene
            ("Clc1ccccc1", "[R1]c1ccccc1"),
            ("Brc1ccccc1", "[R1]c1ccccc1"),
            ("Ic1ccccc1", "[R1]c1ccccc1"),
        ]

        for input, expected in tests:
            with self.subTest(input=input, expected=expected):
                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["Suzuki_Halogen"], rGroup="R1")

                rxn.react(m)
                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertNotEqual(input, result)
                self.assertEqual(expected, result)

    def test_if_Suzuki_Halogen_is_specific_enough(self):
        tests = [
            # Prevent Fluorbenzene
            ("Fc1ccccc1"),
        ]

        for smiles in tests:
            with self.subTest(smiles=smiles):
                input = smiles

                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["Suzuki_Halogen"], rGroup="R1")

                with self.assertRaises(ReactionNoProductException):
                    rxn.react(m)

                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertEqual(input, result)

    def test_if_Suzuki_Boron_gives_correct_products(self):
        tests = [
            # Phenylboronic acid
            ("OB(O)c1ccccc1", "[R1]c1ccccc1"),
            # Phenylboronic acid pinacol ester
            ("CC1(C)C(C)(C)OB(O1)c1ccccc1", "[R1]c1ccccc1"),
        ]

        for input, expected in tests:
            with self.subTest(input=input, expected=expected):
                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["Suzuki_Boron"], rGroup="R1")

                rxn.react(m)
                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertNotEqual(input, result)
                self.assertEqual(expected, result)

    def test_if_Suzuki_Boron_is_specific_enough(self):
        tests = [
            # Prevent nothing?
        ]

        for smiles in tests:
            with self.subTest(smiles=smiles):
                input = smiles

                m = Molecule(input)
                rxn = Reaction(**predefinedReactions["Suzuki_Boron"], rGroup="R1")

                with self.assertRaises(ReactionNoProductException):
                    rxn.react(m)

                result = unmask(Chem.MolToSmiles(m._mol))

                self.assertEqual(input, result)