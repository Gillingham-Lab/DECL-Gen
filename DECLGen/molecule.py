from typing import Dict, List
from rdkit import Chem

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