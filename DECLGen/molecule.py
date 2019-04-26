from typing import Dict, List, Union
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Lipinski, Descriptors, QED, MolSurf

from .exceptions import MoleculeInvalidSmilesException


def property_title(name):
    """
    Decorator to add a column header to a molecular property.
    :param name:
    :return:
    """
    def property_title_decorator(func):
        func.title = name
        return func

    return property_title_decorator


class Molecule:
    """
    A container for a Molecule.

    Abstracts the molecular graph from rdkit and the corresponding smiles structure and enables easy access
    to common rdkit properties.
    """
    _mol = None
    _smiles = None

    def __init__(self, smiles: str):
        """
        Constructor. Accepts a smiles and creates a molecule.
        :param smiles:
        :raises MoleculeInvalidSmilesException if the smiles is invalid.
        """
        self._smiles = smiles
        self._mol = Chem.MolFromSmiles(smiles)

        if self._mol is None:
            raise MoleculeInvalidSmilesException

    @staticmethod
    def get_data_headers(data_fields: Dict[str, bool]) -> List[str]:
        """
        Static methods to convert a dict with method names to a list of column headers.

        This method accepts a dictionary in the style of {"method_name": True} and converts it to a list of
        column headers based on the name given the method with the property_title decorator.

        *Important:*
        This method relies on the dictionary keeping the sequence constant. This requires at least Python 3.6 or an
        ordered dict instead.

        :param data_fields:
        :return:
        """
        return_list = []
        for key in data_fields:
            if data_fields[key] is True:
                method = getattr(Molecule, key)
                return_list.append(method.title)

        return return_list

    def get_data(self, data_fields: Dict[str, bool]) -> List[str]:
        """
        Method to convert a dict with method names to a list of the actual properties.

        This methods accepts a dictionary in the style of {"method_name": True} and converts it to a list of
        the calculates properties corresponding to the method name.

        *Important:*
        This method relies on the dictionary keeping the sequence constant. This requires at least Python 3.6 or an
        ordered dict instead.

        :param data_fields:
        :return:
        """
        return_list = []
        for key in data_fields:
            if data_fields[key] is True:
                value = getattr(self, key)()
                return_list.append(value)

        return return_list

    @property_title("Canonical smiles")
    def canonical_smiles(self) -> str:
        """
        :return: Canonical SMILES.
        """

        # Note: Newer rdkit versions do not require the isomericSmiles-parameter as it True by default
        # Older versions have if False by default, requiring us to specify it here.
        return Chem.MolToSmiles(self._mol, isomericSmiles=True)

    @property_title("MW")
    def mw(self) -> float:
        """
        :return: Molecular Weight.
        """
        return Descriptors.MolWt(self._mol)

    @property_title("QED")
    def qed(self) -> float:
        """
        :return: Quantification of the Estimation of Druglike properties.
        """
        return QED.qed(self._mol)

    @property_title("TPSA")
    def tpsa(self) -> Union[int, float]:
        """
        :return: Topological Surface Area
        """
        return AllChem.CalcTPSA(self._mol)

    @property_title("LabuteASA")
    def labute_asa(self):
        """
        :return: Labute Accessible Surface Area
        """
        return MolSurf.pyLabuteASA(self._mol, 1)

    @property_title("TPSApMW")
    def tpsapermw(self) -> float:
        """
        Calculates the ratio of the TPSA to molecular weight to normalize TPSA for big molecules.

        This is a measure sometimes use to describe macrocycles.
        :return: Ratio of TPSA to MW
        """
        return self.tpsa()/self.mw()

    @property_title("ALogP")
    def alogp(self) -> Union[int, float]:
        """
        Calculates the ALogP, an atom based estimation of the LogP.
        :return: Crippen Atomic LogP
        """
        return Crippen.MolLogP(self._mol)

    @property_title("NumHDonors")
    def hdonors(self) -> Union[int, float]:
        """
        :return: Number of hydrogen bond donors.
        """
        return Lipinski.NumHDonors(self._mol)

    @property_title("NumHAcceptors")
    def hacceptors(self) -> Union[int, float]:
        """
        :return: Number of hydrogen bond acceptors.
        """
        return Lipinski.NumHAcceptors(self._mol)

    @property_title("NumHeteroatoms")
    def nhetero(self) -> Union[int, float]:
        """
        :return: Number of heteroatoms (all atoms except H, C)
        """
        return Lipinski.NumHeteroatoms(self._mol)

    @property_title("NumRotatableBonds")
    def rotatable(self) -> Union[int, float]:
        """
        :return: Number of bonds that can freely rotate.
        """
        return Lipinski.NumRotatableBonds(self._mol)

    @property_title("NOCount")
    def no(self) -> Union[int, float]:
        """
        :return: Number if nitrogen and oxygen (hydrogen bond acceptors).
        """
        return Lipinski.NOCount(self._mol)

    @property_title("NHOHCount")
    def nhoh(self) -> Union[int, float]:
        """
        :return: Number of nitrogen and oxygen bound with hydrogen (hydrogen bond donors).
        """
        return Lipinski.NHOHCount(self._mol)

    @property_title("RingCount")
    def rings(self) -> Union[int, float]:
        """
        :return: Number of unique rings in a molecule.
        """
        ringInfo = self._mol.GetRingInfo().AtomRings()
        return len(ringInfo)

    @property_title("BiggestRingSize")
    def maxRingSize(self) -> Union[int, float]:
        """
        :return: The size of the biggest ring. Important measure for macrocycles.
        """
        ringInfo = self._mol.GetRingInfo().AtomRings()
        if len(ringInfo) > 0:
            return max([len(t) for t in ringInfo])
        else:
            return 0

    @property_title("FractionCsp3")
    def csp3(self) -> Union[int, float]:
        """
        :return: Fraction of carbon atoms that are sp3 hybridizes.
        """
        return Lipinski.FractionCSP3(self._mol)

    @property_title("NumHeavyAtoms")
    def heavyAtoms(self) -> Union[int, float]:
        """
        :return: Number of heavy atoms (atoms bigger than hydrogen).
        """
        return Lipinski.HeavyAtomCount(self._mol)