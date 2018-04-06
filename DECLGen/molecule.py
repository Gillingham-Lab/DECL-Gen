from rdkit import Chem

class Molecule:
    mol = None
    smiles = None

    def __init__(self, smiles):
        self.smiles = smiles
        self.mol = Chem.MolFromSmiles(smiles)