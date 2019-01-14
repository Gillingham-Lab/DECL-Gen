import argh
from DECLGen import Runtime
from DECLGen.exceptions import \
    DECLException, \
    LibraryElementExistsException
from DECLGen.cli.helpers import ProgressBar

@argh.arg("--for-export", default=False)
def elm_list(
    id: "Category identifier",
    for_export: "Format the output for exporting the list." = False,
):
    """ Lists all elements of a given diversity element category. """
    r = Runtime()

    try:
        cat = r.storage.library.get_category(id)

        if for_export is False:
            print("{t.bold}{a}\t{b}\t{c}{t.normal}".format(a="index", b="codon", c="smiles", t=r.t))
            for element in cat:
                print("{}\t{}\t{}".format(element.index, element.codon, element.raw_smiles))
        else:
            for element in cat:
                print("{}\t{}".format(element.codon, element.raw_smiles))
    except DECLException as e:
        r.error_exit(e)


def elm_show(
        id: "Category identifier",
        index: "Element index."
):
    """ Shows a given element of a given diversity element category. """
    r = Runtime()

    try:
        cat = r.storage.library.get_category(id)
        elm = cat.get_element(index)

        print("{t.bold}{a}\t{b}\t{c}{t.normal}".format(a="index", b="codon", c="smiles", t=r.t))
        print("{}\t{}\t{}".format(elm.index, elm.codon, elm.raw_smiles))
    except DECLException as e:
        r.error_exit(e)


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
    """ Removes a diversity element with a given index."""
    r = Runtime()

    try:
        cat = r.storage.library.get_category(id)
        cat.del_element(index)
    except DECLException as e:
        r.error_exit(e)

    r.save()


def elm_replace():
    pass

@argh.arg("anchorTranslations", nargs="+")
def elm_copy(
    idFrom: "Category identifier to copy into.",
    idInto: "Category identifier to copy from.",
    anchorTranslations: "A list of R1:R3 pairs to translate the R-groups (In this case, R1 will be replaced with R3)",
):
    """ Copies diversity elements from a category into the current one. Overwrites all target cat must be empty."""
    r = Runtime()

    try:
        catInto = r.storage.library.get_category(idInto)
        catFrom = r.storage.library.get_category(idFrom)

        anchors = {}
        for anchor in anchorTranslations:
            a = anchor.split(":")
            aFrom, aInto = a[0], a[1]
            anchors[aFrom] = aInto

        print("Copying...")
        progressBar = ProgressBar(r.t)
        progressBar.start()
        imported = catInto.copy_elements_from(catFrom, anchors, updateable=progressBar)
        progressBar.finish()
        print("Added {} compounds".format(imported))

    except DECLException as e:
        print()
        r.error_exit(e)

    r.save()


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
