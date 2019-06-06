import argh
import sys
from DECLGen import Runtime
from DECLGen.exceptions import \
    DECLException, \
    LibraryCategoryException, \
    LibraryCategoryExistsException, \
    LibraryCategoryNotFoundException

@argh.arg("anchors", nargs="+")
def ssp_add(
    id1: "Desired primary category identifier. Used to identify the superset in the DNA strand.",
    id2: "Desired secondary category identifier. Used to identify the compound in the DNA strand, depending on id1.",
    name: "A human-readable name for this superset category",
    anchors: "A list of R-Groups that MUST elements for this superset-category must have",
    codon1_length: "The desired length of the codon for identifying the subset. Codon length is flexible if not given or if set to 0." = 0,
    codon2_length: "The desired length of the codon for identifying individual compounds of a subset." = 0,
    codon1_reverse: "Set to 1 if codon 1 should be read as its reverse complement in the DNA" = 0,
    codon2_reverse: "Set to 1 if codon 2 should be read as its reverse complement in the DNA" = 0,
    force_yes: "Set to true if answers should be assumed to be yes" = False,
) -> None:
    """ Adds a new superset diversity point (SSP).

    SSPs can be used to distinguish different subsets of a library.
    """
    r = Runtime()

    try:
        r.storage.library.add_superset_category(id1, id2, name, anchors,
                                                codon1_length=codon1_length, codon2_length=codon2_length)
        r.storage.library.get_category(id1).set_reverse_complement(True if codon1_reverse is 1 else False)
        r.storage.library.get_category(id2).set_reverse_complement(True if codon2_reverse is 1 else False)
    except LibraryCategoryExistsException:
        if not force_yes:
            answer = r.t.decide(
                "One or more of the ids you've selected are already occupied. "
                "This will overwrite the existing categories and delete stored diversity elements." + \
                "If you want to edit some parts, use {exec} DEC-edit instead".format(exec=sys.argv[0]),
                "Proceed anyway?",
                ["Y", "n"]
            )

            if answer == "Y":
                if r.storage.library.has_category(id1):
                    r.storage.library.del_category(id1)
                if r.storage.library.has_category(id2):
                    r.storage.library.del_category(id2)

                r.storage.library.add_superset_category(id1, id2, name, anchors,
                                                        codon1_length=codon1_length, codon2_length=codon2_length)
                r.storage.library.get_category(id1).set_reverse_complement(True if codon1_reverse is 1 else False)
                r.storage.library.get_category(id1).set_reverse_complement(True if codon2_reverse is 1 else False)
            else:
                print("Superset category was not added.")
    except DECLException as e:
        r.error_exit(e)

    r.save()


def ssp_pnt_list(
        ssid: "Superset category identifier.",
) -> None:
    """ Lists all registered categories of a superset diversity point."""
    r = Runtime()

    try:
        ssc = r.storage.library.get_category(ssid)
        if ssc.is_superset() is False:
            raise DECLException("Given identifier must belong to a superset category.")

        table = r.t.table((10, 10, 10, 40), first_column=True, first_row=True)
        table.add_row("Index", "Codon", "Lengh", "Name")

        for subcat in ssc:
            table.add_row(subcat.id, subcat.codon(ssc), len(subcat), subcat.name)

        table.display()
    except LibraryCategoryException as e:
        r.error_exit(e)


def ssp_pnt_add(
        ssid: "Superset category identifier.",
        name: "A human-readable name for the category.",
        index: "The desired codon or index to encode this sub category. Can also get automatically generated." = None,
) -> None:
    """ Adds a diversity point with index codon to an existing SSP. """
    r = Runtime()

    try:
        ssc = r.storage.library.get_category(ssid)
        if ssc.is_superset() is False:
            raise DECLException("Given identifier must belong to a superset category.")

        try:
            ssc.add_category(name, index)
        except LibraryCategoryExistsException:
            answer = r.t.decide(
                "This category index already exists. This will overwrite the category and remove all elements.",
                "Proceed anyway?",
                ["Y", "n"]
            )

            if answer == "Y":
                ssc.del_category(index)
                ssc.add_category(name, index)
            else:
                print("No sub category added.")
    except DECLException as e:
        r.error_exit(e)

    r.save()

def ssp_pnt_reindex(
        ssid: "Superset category identifier.",
        current_index: "The current codon or index of a sub category.",
        new_index: "The new codon or index of the sub category. Must not exist yet.",
) -> None:
    """ Moves a diversity point from a SSP to a newly given index. """
    r = Runtime()

    try:
        ssc = r.storage.library.get_category(ssid)
        if ssc.is_superset() is False:
            raise DECLException("Given identifier must belong to a superset category.")

        ssc.reindex_category(current_index, new_index)
    except DECLException as e:
        r.error_exit(e)

    ssp_pnt_list(ssid)

    r.save()

def ssp_pnt_del(
        ssid: "Superset category identifier.",
        index: "The desired codon or index to encode this sub category."
) -> None:
    """ Removes a diversity point from an existing SSP. """
    r = Runtime()

    try:
        ssc = r.storage.library.get_category(ssid)
        if ssc.is_superset() is False:
            raise DECLException("Given identifier must belong to a superset category.")

        ssc.del_category(index)
    except DECLException as e:
        r.error_exit(e)

    r.save()