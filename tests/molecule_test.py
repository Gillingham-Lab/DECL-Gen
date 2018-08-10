import unittest
from rdkit import RDLogger

from DECLGen.molecule import Molecule, property_title
from DECLGen.exceptions import MoleculeInvalidSmilesException

class MoleculeTestCase(unittest.TestCase):
    def test_if_creating_molecule_with_valid_smiles_creates_a_molecule(self):
        m = Molecule("CCO")

        self.assertIsNotNone(m)

    def test_if_creating_molecule_with_invalid_smiles_raises_an_exception(self):
        # Silence RDKit logger
        logger = RDLogger.logger()
        logger.setLevel(RDLogger.CRITICAL)

        with self.assertRaises(MoleculeInvalidSmilesException):
            m = Molecule("CC[R1]")

    def test_molecule_canonical_smiles(self):
        m = Molecule("CCO")
        self.assertEqual(m.canonical_smiles(), "CCO")

    def test_molecule_molecular_weight(self):
        self.assertIsNotNone(Molecule.mw.title)

        m = Molecule("CCO")
        r = m.mw()
        self.assertIsInstance(r, float)
        self.assertAlmostEqual(46, r, 0)

    def test_molecule_qed(self):
        self.assertIsNotNone(Molecule.qed.title)

        m = Molecule("CCO")
        r = m.qed()
        self.assertIsInstance(r, float)
        self.assertAlmostEqual(0.4068, r, 3)

    def test_molecule_tpsa(self):
        self.assertIsNotNone(Molecule.tpsa.title)

        m = Molecule("CCO")
        r =  m.tpsa()
        self.assertIsInstance(r, float)
        self.assertAlmostEqual(20, r, 0)

    def test_molecule_labute_asa(self):
        self.assertIsNotNone(Molecule.labute_asa.title)

        m = Molecule("CCO")
        r = m.labute_asa()
        self.assertIsInstance(r, float)
        self.assertAlmostEqual(20, r, 0)

    def test_molecule_tpsapermw(self):
        self.assertIsNotNone(Molecule.tpsapermw.title)

        m = Molecule("CCO")
        r = m.tpsapermw()
        self.assertIsInstance(r, float)
        self.assertAlmostEqual(m.tpsa()/m.mw(), r, 5)

    def test_molecule_alogp(self):
        self.assertIsNotNone(Molecule.alogp.title)

        m = Molecule("CCO")
        r = m.alogp()
        self.assertIsInstance(r, float)
        self.assertAlmostEqual(0, r, 0)

    def test_molecule_hdonors(self):
        self.assertIsNotNone(Molecule.hdonors.title)

        m = Molecule("CCO")
        r = m.hdonors()
        self.assertIsInstance(r, int)
        self.assertEqual(1, r)

    def test_molecule_hacceptors(self):
        self.assertIsNotNone(Molecule.hacceptors.title)

        m = Molecule("CCO")
        r = m.hacceptors()
        self.assertIsInstance(r, int)
        self.assertEqual(1, r)

    def test_molecule_nhetero(self):
        self.assertIsNotNone(Molecule.nhetero.title)

        m = Molecule("CCO")
        r = m.nhetero()
        self.assertIsInstance(r, int)
        self.assertEqual(1, r)

    def test_molecule_rotatable(self):
        self.assertIsNotNone(Molecule.rotatable.title)

        m = Molecule("CC=O")
        r = m.rotatable()
        self.assertIsInstance(r, int)
        self.assertEqual(0, r)

    def test_molecule_no(self):
        self.assertIsNotNone(Molecule.no.title)

        m = Molecule("CC=O")
        r = m.no()
        self.assertIsInstance(r, int)
        self.assertEqual(1, r)

    def test_molecule_nhoh(self):
        self.assertIsNotNone(Molecule.nhoh.title)

        m = Molecule("CC=O")
        r = m.nhoh()
        self.assertIsInstance(r, int)
        self.assertEqual(0, r)

    def test_molecule_rings(self):
        self.assertIsNotNone(Molecule.rings.title)

        m = Molecule("CC=O")
        r = m.rings()
        self.assertIsInstance(r, int)
        self.assertEqual(0, r)

        self.assertEqual(1, Molecule("C1CC1").rings())
        self.assertEqual(4, Molecule("C1CC12CCC23CCCC34CCCCCC4").rings())

    def test_molecule_maxRingSize(self):
        self.assertIsNotNone(Molecule.maxRingSize.title)

        m = Molecule("CC=O")
        r = m.maxRingSize()
        self.assertIsInstance(r, int)
        self.assertEqual(0, r)

        self.assertEqual(3, Molecule("C1CC1").maxRingSize())
        self.assertEqual(7, Molecule("C1CC12CCC23CCCC34CCCCCC4").maxRingSize())

    def test_molecule_csp3(self):
        self.assertIsNotNone(Molecule.csp3.title)

        m = Molecule("CC=O")
        r = m.csp3()
        self.assertIsInstance(r, float)
        self.assertAlmostEqual(0.50001, r, 4)

    def test_molecule_heavyAtoms(self):
        self.assertIsNotNone(Molecule.heavyAtoms.title)

        m = Molecule("CC=O")
        r = m.heavyAtoms()
        self.assertIsInstance(r, int)
        self.assertEqual(3, r)

    def test_molecule_get_data_returns_specific_fields(self):
        fields = [
            "canonical_smiles",
            "mw",
            "qed",
            "tpsa",
            "tpsapermw",
            "labute_asa",
            "alogp",
            "hdonors",
            "hacceptors",
            "nhetero",
            "rotatable",
            "no",
            "nhoh",
            "rings",
            "maxRingSize",
            "csp3",
            "heavyAtoms",
        ]
        n = len(fields)

        m = Molecule("CC=O")

        for x in range(len(fields)):
            for i in range(n-x):
                subFields = fields[x:x+i+1]

                with self.subTest(i=(i,x+i+1,subFields)):
                    data = m.get_data({x: True for x in subFields})
                    headers = m.get_data_headers({x: True for x in subFields})

                    self.assertEqual(len(data), len(headers))
                    self.assertEqual(len(data), i+1)

    def test_if_decorator_decorates_properly(self):
        @property_title("A")
        def test(a):
            return 0

        self.assertEqual("A", test.title)
