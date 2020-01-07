from typing import List
from .molecule import Molecule

class MoleculeContainer:
    molecules: List[Molecule]

    def __init__(self):
        self.molecules = []

    def __len__(self):
        return len(self.molecules)

    @classmethod
    def fromSmilesList(cls, source) -> "MoleculeContainer":
        """ Reads a file and constructs molecules out of it. """
        container = cls()

        with open(source, "r") as file:
            for row in file:
                row = row.strip()

                mol = Molecule(row)

                container.add(mol)

        return container

    def count(self):
        """ Returns the amount of molecules in the list. """
        return len(self)

    def add(self, molecule: Molecule):
        """ Adds a molecule to the container. """
        self.molecules.append(Molecule)