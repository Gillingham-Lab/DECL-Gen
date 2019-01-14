import argh
import sys
from DECLGen import Runtime
from DECLGen.exceptions import \
    DECLException, \
    LibraryCategoryException, \
    LibraryCategoryExistsException, \
    LibraryCategoryNotFoundException


def cat_list():
    """ Lists all registered diversity element categories """
    r = Runtime()

    try:
        cats = r.storage.library.get_categories()

        if len(cats) == 0:
            print("{t.red}No categories defined.{t.normal}".format(t=r.t))
            return
        table_head = "{t.underline}{t.bold}{id:<10} {name:<50} {anchors:<10}{t.normal}"
        table_entry = "{t.bold}{id:<10}{t.normal} {name:<50} {anchors:<10}"
        print(table_head.format(id="id", name="Name", anchors="R-Groups", t=r.t))
        for cat in cats:
            print(table_entry.format(**cat.describe(), t=r.t))
    except LibraryCategoryException as e:
        r.error_exit(e)


def cat_show(id: "Category identifier"):
    """ Shows details about a specific diversity element category """
    r = Runtime()

    try:
        cat = r.storage.library.get_category(id)
        cat_desc = cat.describe()
        print("{t.bold}Diversity Element Category <{id}>{t.normal}".format(t=r.t, **cat_desc))
        print("  {t.bold}{title:<20}{t.normal}{value:>20}".format(
            t=r.t, title="R-Groups", value=cat_desc["anchors"]))
        print("  {t.bold}{title:<20}{t.normal}{value:>20}".format(
            t=r.t, title="Codon length", value=cat_desc["codon_length"]))
        print("  {t.bold}{title:<20}{t.normal}{value:>20}".format(
            t=r.t, title="Number of elements", value=cat_desc["elements"]))

        if type(cat_desc["max_elements"]) is int:
            if int(cat_desc["elements"]) < int(cat_desc["max_elements"]):
                print("  {t.bold}{title:<20}{t.normal}{t.green}{value:>20}{t.normal}".format(
                    t=r.t, title="Max. of elements", value=cat_desc["max_elements"]))
            elif int(cat_desc["elements"]) == int(cat_desc["max_elements"]):
                print("  {t.bold}{title:<20}{t.normal}{t.green}{value:>20}{t.normal}".format(
                    t=r.t, title="Max. of elements", value=cat_desc["max_elements"]))
            else:
                print("  {t.bold}{title:<20}{t.normal}{t.green}{value:>20}{t.normal}".format(
                    t=r.t, title="Max. of elements", value=cat_desc["max_elements"]))
        else:
            print("  {t.bold}{title:<20}{t.normal}{t.green}{value:>20}{t.normal}".format(
                t=r.t, title="Max. of elements", value=cat_desc["max_elements"]))

        if "min_codon_length" in cat_desc:
            print("  {t.bold}{title:<20}{t.normal}{value:>20}".format(
                t=r.t, title="Codon length needed", value=cat_desc["min_codon_length"]))

        if "reverse_complement" in cat_desc:
            print("  {t.bold}{title:<20}{t.normal}{value:>20}".format(
                t=r.t, title="Reverse complement", value=cat_desc["reverse_complement"]))
    except DECLException as e:
        r.error_exit(e)


@argh.arg("anchors", nargs="+")
def cat_add(
    id: "Desired category identifier",
    name: "A human-readable name for this category",
    anchors: "A list of R-Groups that MUST elements for this category must have",
    codon_length: "The desired length of the codon. Codon length is flexible if not given or if set to 0." = 0,
    reverse_complement: "Set to 1 if the codon should be put as its reverse complement into the DNA" = 0,
    force_yes: "Set to true if answers should be assumed to be yes" = False,
) -> None:
    """ Adds a new diversity element category """
    r = Runtime()

    try:
        r.storage.library.add_category(id, name, anchors, codon_length)
        r.storage.library.get_category(id).set_reverse_complement(True if reverse_complement is 1 else False)
    except LibraryCategoryExistsException:
        if not force_yes:
            print("{t.red}Warning:{t.normal} This will overwrite the existing DEC and delete all saved DE.".format(t=r.t))
            print("If you want to edit some parts, use {exec} DEC-edit instead.".format(exec=sys.argv[0]))

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
            print("Category was not added.")
    except DECLException as e:
        r.error_exit(e)

    r.save()


def cat_del(id: "Category identifier"):
    """ Removes an existing diversity element category """
    r = Runtime()

    try:
        r.storage.library.del_category(id)
    except DECLException as e:
        r.error_exit(e)

    r.save()


def cat_edit(
    id: "Category identifier",
    name: "A human-readable name for this category" = None,
    codon_length: "The desired length of the codon. Codon length is flexible if not given or if set to 0." = None,
    reverse_complement: "Set to 1 if the codon should be put as its reverse complement into the DNA and 0 if not." = -1,
) -> None:
    """ Changes some available parameters for a diversity element category """
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


def cat_clear(id: "Category identifier"):
    """Removes all diversity elements from a given diversity element category"""
    r = Runtime()

    try:
        cat = r.storage.library.get_category(id)
        cat.clear()
    except DECLException as e:
        r.error_exit(e)

    r.save()