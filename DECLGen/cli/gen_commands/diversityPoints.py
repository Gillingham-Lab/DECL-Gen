import argh
import sys
from DECLGen import Runtime
from DECLGen.exceptions import \
    DECLException, \
    LibraryCategoryException, \
    LibraryCategoryExistsException, \
    LibraryCategoryNotFoundException


def pnt_list():
    """ Lists all diversity points. """
    r = Runtime()

    try:
        cats = r.storage.library.get_categories()

        if len(cats) == 0:
            r.warning("No diversity points defined.")
            return

        table = r.t.table((10, 8, 40, 10, 20), first_column=True, first_row=True)
        table.add_row("id", "Length", "Name", "Anchors", "Extra")

        for cat in cats:
            cat_desc = cat.describe()

            if cat.is_subset():
                cat_desc["id"] = "({})".format(cat_desc["id"])
                cat_desc["anchors"] = "-"

            table.add_row(cat_desc["id"], len(cat), cat_desc["name"], cat_desc["anchors"], cat_desc["extra"])

            if cat.is_superset():
                for subcat in cat:
                    subcat_desc = subcat.describe()
                    subcat_desc["id"] = "{}.{}".format(cat_desc["id"], subcat.codon(cat))
                    table.add_row(subcat_desc["id"], len(subcat), subcat_desc["name"], subcat_desc["anchors"], subcat_desc["extra"])

        table.display()
    except LibraryCategoryException as e:
        r.error_exit(e)


def pnt_show(id: "Category identifier"):
    """ Shows details about a specific diversity point. """
    r = Runtime()

    try:
        cat = r.storage.library.get_category(id)
        cat_desc = cat.describe()

        print(r.t.highlight("Diversity point {id}: {name}".format(**cat_desc)))

        dl = r.t.dl(20, 20, highlight_key=True, list_item="", values_justify_right=True)
        dl.add_row("R-Groups", cat_desc["anchors"])
        dl.add_row("Codon length", cat_desc["codon_length"])
        dl.add_row("Number of elements", cat_desc["elements"])

        if type(cat_desc["max_elements"]) is int:
            if int(cat_desc["elements"]) < int(cat_desc["max_elements"]):
                format = r.t.val_good
            elif int(cat_desc["elements"]) == int(cat_desc["max_elements"]):
                format = r.t.val_good
            else:
                format = r.t.val_bad
        else:
            format = r.t.val_good

        dl.add_row("Max. of elements", cat_desc["max_elements"], format)

        if "min_codon_length" in cat_desc:
            dl.add_row("Codon length needed", cat_desc["min_codon_length"])
        if "reverse_complement" in cat_desc:
            dl.add_row("Reverse complement", cat_desc["reverse_complement"])

        if "extra" in cat_desc:
            dl.add_row("Extra", cat_desc["extra"])


        dl.display()

    except DECLException as e:
        r.error_exit(e)


@argh.arg("anchors", nargs="+")
def pnt_add(
    id: "Desired category identifier",
    name: "A human-readable name for this category",
    anchors: "A list of R-Groups that MUST elements for this category must have",
    codon_length: "The desired length of the codon. Codon length is flexible if not given or if set to 0." = 0,
    reverse_complement: "Set to 1 if the codon should be put as its reverse complement into the DNA" = 0,
    force_yes: "Set to true if answers should be assumed to be yes" = False,
) -> None:
    """ Adds a new diversity point. """
    r = Runtime()

    try:
        r.storage.library.add_category(id, name, anchors, codon_length)
        r.storage.library.get_category(id).set_reverse_complement(True if reverse_complement is 1 else False)
    except LibraryCategoryExistsException:
        if not force_yes:
            print("{t.red}Warning:{t.normal} This will overwrite the existing diversity point and delete all saved elements within that point.".format(t=r.t))
            print("If you want to edit some parts, use {exec} pnt-edit instead.".format(exec=sys.argv[0]))

            while True:
                answer = input("Proceed anyway (Y/n)? ")
                answer = answer[0]
                if answer == "Y" or answer == "n":
                    break
        else:
            answer = "Y"

        if answer == "Y":
            r.storage.library.del_category(id)
            r.storage.library.add_category(id, name, anchors, codon_length)
        else:
            print("Diversity point was not added.")
    except DECLException as e:
        r.error_exit(e)

    r.save()


def pnt_del(id: "Category identifier"):
    """ Removes an existing diversity point. """
    r = Runtime()

    try:
        r.storage.library.del_category(id)
    except DECLException as e:
        r.error_exit(e)

    r.save()


def pnt_edit(
    id: "Category identifier",
    name: "A human-readable name for this category" = None,
    codon_length: "The desired length of the codon. Codon length is flexible if not given or if set to 0." = None,
    reverse_complement: "Set to 1 if the codon should be put as its reverse complement into the DNA and 0 if not." = -1,
) -> None:
    """ Changes some available parameters for a diversity point. """
    r = Runtime()

    try:
        cat = r.storage.library.get_category(id)

        if name is not None:
            cat.set_name(name)
        if codon_length is not None:
            cat.set_codon_length(int(codon_length))

        if reverse_complement > 0:
            cat.set_reverse_complement(True if reverse_complement >= 1 else False)
    except DECLException as e:
        r.error_exit(e)

    r.save()


def pnt_clear(id: "Category identifier"):
    """ Removes all diversity elements from a given diversity point. """
    r = Runtime()

    try:
        cat = r.storage.library.get_resolved_category(id)
        cat.clear()
    except DECLException as e:
        r.error_exit(e)

    r.save()