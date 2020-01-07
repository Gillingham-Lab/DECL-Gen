from typing import List
from rdkit import Chem
from rdkit.Chem import AllChem

from .exceptions import ReactionInvalidSmartsException, ReactionNoProductException
from .molecule import Molecule
from .rtools import RGroup, unmask, mask

class Reaction:
    smartsStrings: List[str]
    protectionStrings: List[str]
    rGroup: str

    smarts: List["rdkit.Chem.rdChemReactions.ChemicalReaction"]
    protections: List["rdkit.Chem.rdChemReactions.ChemicalReaction"]

    def __init__(self, smarts, protections=None, rGroup=None):
        if protections is None:
            protections = []

        if type(smarts) == str:
            smarts = [smarts]

        self.smartsStrings = smarts
        self.protectionStrings = protections

        if rGroup is not None:
            self.rGroup = RGroup.fromR(rGroup)

            for i in range(len(self.smartsStrings)):
                S = self.smartsStrings[i]
                self.smartsStrings[i] = S.replace("[R]", self.rGroup.Au)

        self.smarts = self.createReactions(self.smartsStrings)
        self.protections = self.createSubstructs(self.protectionStrings)

    def createN(self, smarts: List[str], method) -> List:
        """ Converts a list of smarts to a list of reactions. """
        reactions = []

        for S in smarts:
            S = mask(S)

            reaction = method(S)

            if reaction is None:
                raise ReactionInvalidSmartsException(S)

            reactions.append(reaction)

        return reactions

    def createReactions(self, smarts: List[str]) -> List["rdkit.Chem.rdChemReactions.ChemicalReaction"]:
        """ Converts a list of smarts to a list of reactions. """
        return self.createN(smarts, AllChem.ReactionFromSmarts)

    def createSubstructs(self, smarts: List[str]) -> List:
        """ Converts a list of smarts to a list of structures for substruct match. """
        return self.createN(smarts, Chem.MolFromSmarts)

    def react(self, molecule: Molecule):
        mol = molecule._mol
        productSet = None

        # Add protections
        for p in self.protections:
            matches = mol.GetSubstructMatches(p)

            for match in matches:
                mol.GetAtomWithIdx(match[0]).SetProp("_protected", "1")

        # Run all given reactions
        for rxn in self.smarts:
            products = rxn.RunReactants((mol, ))

            if len(products) > 0:
                productSet = products[0]

        if productSet is None:
            raise ReactionNoProductException(molecule._smiles, self.smartsStrings)

        # Replace old with new molecule
        molecule._mol = productSet[0]





predefinedReactions = {
    "AmideCoupling_Amine": {
        "smarts": [
            "[CX4:2][NX3;H3:1]>>[C:2][N:1][R]",
            "[CX4:2][NX3;H2:1]>>[C:2][N:1][R]",
            "[CX4:2][NX3;H1:1]>>[C:2][N:1][R]",
        ],
        "protections": [
            "[N;$(NC=[O,S])]",
        ],
    },

    "AmideCoupling_Acid": {
        "smarts": [
            "[C:1](=[O:2])-[OH1]>>[C:1](=[O:2])"
        ],
    },

    "CuAAC_Alkyne": {
        "smarts": [
            "[C:1]#[CH1:2]>>[c:1]1[c:2]n([R])nn1"
        ],
    },

    "CuAAC_Azide": {
        "smarts": [
            "[#6:1]-[N]=[N+]=[N-]>>[*:1]-[R]"
        ]
    },
}