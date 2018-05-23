from DECLGen import Runtime
from DECLGen.exceptions import \
    DECLException, \
    LibraryElementExistsException
from DECLGen.cli.helpers import ProgressBar


def elm_list(id: "Category identifier"):
    r = Runtime()

    try:
        cat = r.storage.library.get_category(id)

        print("{t.bold}{a}\t{b}\t{c}{t.normal}".format(a="index", b="codon", c="smiles", t=r.t))
        for element in cat:
            print("{}\t{}\t{}".format(element.index, element.codon, element.raw_smiles))
    except DECLException as e:
        r.error_exit(e)


def elm_show():
    pass


def elm_add(
    id: "Category identifier",
    elm_smiles: "Element SMILES with the required R-Groups",
    index: "Element index. Automatically generated if not given. Can be either Number or DNA tag." = None,
    force_yes: "Force yes" = False,
):
    """ Adds a new element to an existing category """
    r = Runtime()

    try:
        index = int(index)
    except TypeError:
        pass

    try:
        cat = r.storage.library.get_category(id)

        try:
            cat.add_element(elm_smiles, index)
        except LibraryElementExistsException:
            if not force_yes:
                print("{t.red}Warning:{t.normal} This will overwrite the existing element.".format(t=r.t))

                while True:
                    answer = input("Proceed anyway (Y/n)? ")
                    answer = answer[0]
                    if answer == "Y" or answer == "n":
                        break
            else:
                answer = "Y"

            if answer == "Y":
                cat.del_element(index)
                cat.add_element(elm_smiles, index)
            else:
                print("Category was not added.")
    except DECLException as e:
        r.error_exit(e)

    r.save()


def elm_del(
    id: "Category identifier",
    index: "Element index. Automatically generated if not given. Can be either Number or DNA tag.",
):
    r = Runtime()

    try:
        cat = r.storage.library.get_category(id)
        cat.del_element(index)
    except DECLException as e:
        r.error_exit(e)

    r.save()


def elm_replace():
    pass


def elm_import(
    id: "Category identifier",
    filename: "Tab separated file to import"
):
    """Imports all diversity elements from a tab separated file and adds them to the given category"""
    r = Runtime()

    try:
        cat = r.storage.library.get_category(id)
        print("Importing...")

        progressBar = ProgressBar(r.t)
        progressBar.start()
        imported = cat.import_elements(filename, updateable=progressBar)
        progressBar.finish()

        print("Added {} compounds".format(imported))
    except DECLException as e:
        print()
        r.error_exit(e)

    r.save()
