

class DECLException(Exception):
    exitcode = 1


class LibraryNotInitializedError(DECLException):
    exitcode = 2


class LibraryExistsError(DECLException):
    exitcode = 3


class LibraryPropertiesNotPreCalculated(DECLException):
    exitcode = 4


class LibraryTemplateException(DECLException):
    exitcode = 10


class LibraryNoDNATemplateException(LibraryTemplateException):
    exitcode = 11


class LibraryCategoryException(DECLException):
    exitcode = 20


class LibraryCategoryNotFoundException(LibraryCategoryException):
    exitcode = 22


class LibraryCategoryExistsException(LibraryCategoryException):
    exitcode = 23


class LibraryCategoryFullException(LibraryCategoryException):
    exitcode = 24


class LibraryCategoryEmptyException(LibraryCategoryException):
    exitcode = 25


class LibraryCategoryAnchorException(LibraryCategoryException):
    exitcode = 28


class LibraryElementException(DECLException):
    exitcode = 30


class LibraryElementNotFoundException(LibraryCategoryException):
    exitcode = 32


class LibraryElementExistsException(LibraryCategoryException):
    exitcode = 33


class EvaluationException(DECLException):
    exitcode = 100


class EvaluationFileDoesNotExist(DECLException):
    exitcode = 101

class EvaluationInvalidFileFormat(DECLException):
    exitcode = 102


class MoleculeInvalidSmilesException(DECLException):
    smiles: str

    def __init__(self, smiles):
        super(MoleculeInvalidSmilesException, self).__init__()

        self.smiles = smiles

class ReactionInvalidSmartsException(DECLException):
    smarts: str

    def __init__(self, smarts):
        super(ReactionInvalidSmartsException, self).__init__()

        self.smarts = smarts

class ReactionNoProductException(DECLException):
    smiles: str
    smarts: str

    def __init__(self, smiles, smarts):
        super(ReactionNoProductException, self).__init__()

        self.smiles = smiles
        self.smarts = smarts