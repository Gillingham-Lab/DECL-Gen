from typing import Dict, List, Union
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Crippen, Lipinski, Descriptors, QED

def property_title(name):
    def property_title_decorator(func):
        func.title = name
        return func

    return property_title_decorator

class Molecule:
    _mol = None
    _smiles = None

    def __init__(self, smiles):
        self._smiles = smiles
        self._mol = Chem.MolFromSmiles(smiles)

    @staticmethod
    def get_data_headers(data_fields: Dict[str, bool]) -> List[str]:
        return_list = []
        for key in data_fields:
            if data_fields[key] is True:
                method = getattr(Molecule, key)
                return_list.append(method.title)

        return return_list

    def get_data(self, data_fields: Dict[str, bool]) -> List[str]:
        return_list = []
        for key in data_fields:
            if data_fields[key] is True:
                value = getattr(self, key)()
                return_list.append(value)

        return return_list

    @property_title("Canonical smiles")
    def canonical_smiles(self) -> str:
        return Chem.MolToSmiles(self._mol, isomericSmiles=True)

    @property_title("MW")
    def mw(self) -> float:
        return Descriptors.MolWt(self._mol)

    @property_title("QED")
    def qed(self) -> float:
        return QED.qed(self._mol)

    @property_title("TPSA")
    def tpsa(self) -> Union[int, float]:
        return AllChem.CalcTPSA(self._mol)

    @property_title("TPSApMW")
    def tpsapermw(self) -> float:
        return self.tpsa()/self.mw()

    @property_title("ALogP")
    def alogp(self) -> Union[int, float]:
        return Crippen.MolLogP(self._mol)

    @property_title("NumHDonors")
    def hdonors(self) -> Union[int, float]:
        return Lipinski.NumHDonors(self._mol)

    @property_title("NumHAcceptors")
    def hacceptors(self) -> Union[int, float]:
        return Lipinski.NumHAcceptors(self._mol)

    @property_title("NumHeteroatoms")
    def nhetero(self) -> Union[int, float]:
        return Lipinski.NumHeteroatoms(self._mol)

    @property_title("NumRotatableBonds")
    def rotatable(self) -> Union[int, float]:
        return Lipinski.NumRotatableBonds(self._mol)

    @property_title("NOCount")
    def no(self) -> Union[int, float]:
        return Lipinski.NOCount(self._mol)

    @property_title("NHOHCount")
    def nhoh(self) -> Union[int, float]:
        return Lipinski.NHOHCount(self._mol)

    @property_title("RingCount")
    def rings(self) -> Union[int, float]:
        ringInfo = self._mol.GetRingInfo().AtomRings()
        return len(ringInfo)

    @property_title("BiggestRingSize")
    def maxRingSize(self) -> Union[int, float]:
        ringInfo = self._mol.GetRingInfo().AtomRings()
        if len(ringInfo) > 0:
            return max([len(t) for t in ringInfo])
        else:
            return 0

    @property_title("FractionCsp3")
    def csp3(self) -> Union[int, float]:
        return Lipinski.FractionCSP3(self._mol)

    @property_title("NumHeavyAtoms")
    def heavyAtoms(self) -> Union[int, float]:
        return Lipinski.HeavyAtomCount(self._mol)