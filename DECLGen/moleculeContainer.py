from rdkit import Chem
from typing import List

from .molecule import Molecule
from .rtools import mask, unmask

class MoleculeContainer:
    molecules: List[Molecule]

    def __init__(self):
        self.molecules = []

    def __len__(self):
        return len(self.molecules)

    def __iter__(self):
        for molecule in self.molecules:
            yield molecule

    @classmethod
    def fromSmilesList(cls, source: str, doMask: bool=True) -> "MoleculeContainer":
        """ Reads a file and constructs molecules out of it. """
        container = cls()

        with open(source, "r") as file:
            for row in file:
                row = row.strip()

                if doMask:
                    row = mask(row)

                mol = Molecule(row)

                container.add(mol)

        return container

    def count(self):
        """ Returns the amount of molecules in the list. """
        return len(self)

    def add(self, molecule: Molecule):
        """ Adds a molecule to the container. """
        self.molecules.append(molecule)

    def export(self, outFile, doUnmask=True):
        with open(outFile, "w") as file:
            for m in self.molecules:
                smiles = m.canonical_smiles()

                if doUnmask:
                    smiles = unmask(smiles)

                file.write(f"{smiles}\n")