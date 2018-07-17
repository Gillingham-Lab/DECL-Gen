

class DECLException(Exception):
    exitcode = 1


class LibraryNotInitializedError(DECLException):
    exitcode = 2


class LibraryExistsError(DECLException):
    exitcode = 3


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