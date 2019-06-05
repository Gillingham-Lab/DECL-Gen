import os
import sys
from pickle import dump, load
from .storage import Storage
from .codon import CodonConfig
from .exceptions import DECLException, LibraryExistsError, LibraryNotInitializedError
from .terminal import Terminal


class Runtime:
    @staticmethod
    def create(
        template: "Template string used for adding in diversity elements",
        advanced_anchors: bool = False,
        superset_categories: bool = False,
    ):
        """ creates a runtime in cwd and returns the instance. """
        if os.path.exists("decl_gen.data"):
            raise LibraryExistsError()

        storage = Storage(
            advanced_anchors=advanced_anchors,
            superset_categories=superset_categories
        )
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
        """
        Constructor.

        Loads the library from decl_gen.data using pickle.load() and initializes the Terminal.
        :param anew:
        """

        if not os.path.exists("decl_gen.data"):
            raise LibraryNotInitializedError()

        with open("decl_gen.data", "rb") as datafile:
            storage = load(datafile)

        self.storage = storage
        CodonConfig.bases = storage.bases
        self.t = Terminal()

    def save(self):
        """
        Saves the library to decl_gen.data using pickle.dump().
        :return:
        """
        with open("decl_gen.data", "wb") as datafile:
            dump(self.storage, datafile, 4)

    def error_exit(self, e: DECLException):
        """
        Exits the Runtime for a DECLException and prints the errorcode as well as the error message.
        :param e:
        :return:
        """
        msg = "{e} (Code: {e.exitcode}".format(e)
        msg = self.t.val_bad(msg)
        print(msg, file=sys.stderr)

        exit(e.exitcode)

    def warning(self, msg: str):
        """
        Prints a warning and exists.
        :param msg:
        :return:
        """
        print(self.t.val_warning(msg))
        exit(0)