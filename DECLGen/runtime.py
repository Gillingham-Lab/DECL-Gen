import os
import sys
from pickle import dump, load
from .storage import Storage
from .codon import CodonConfig
from .exceptions import DECLException, LibraryExistsError, LibraryNotInitializedError
from .terminal import Terminal


class Runtime:
    @staticmethod
    def create(template: "Template string used for adding in diversity elements"):
        """ creates a runtime in cwd and returns the instance. """
        if os.path.exists("decl_gen.data"):
            raise LibraryExistsError()

        storage = Storage()
        storage.set_template(template)

        with open("decl_gen.data", "wb") as datafile:
            dump(storage, datafile, 4)

        return Runtime()

    @staticmethod
    def remove_init_file() -> None:
        """ removes runtime file """
        if not os.path.exists("decl_gen.data"):
            raise LibraryNotInitializedError()

        os.unlink("decl_gen.data")

    def __init__(self, anew=False):
        """ loads a runtime from cwd """

        if not os.path.exists("decl_gen.data"):
            raise LibraryNotInitializedError()

        with open("decl_gen.data", "rb") as datafile:
            storage = load(datafile)

        self.storage = storage
        CodonConfig.bases = storage.bases
        self.t = Terminal()

    def save(self):
        with open("decl_gen.data", "wb") as datafile:
            dump(self.storage, datafile, 4)

    def error_exit(self, e: DECLException):
        print("{t.red}{e}{t.normal}".format(t=self.t, e=e), file=sys.stderr)
        exit(e.exitcode)
