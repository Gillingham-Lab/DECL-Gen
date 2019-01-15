import argh
from DECLGen import Runtime
from DECLGen.exceptions import \
    DECLException, \
    LibraryCategoryException, \
    LibraryCategoryExistsException, \
    LibraryCategoryNotFoundException

@argh.arg("anchors", nargs="+")
def ssc_add(
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
    """ Adds a new superset diversity element category.

    SSCs can be used to distinguish different subsets of a library.
    """
    r = Runtime()

    try:
        r.storage.library.add_superset_category(
            id1, id2, name, anchors,
            codon1_length=codon1_length, codon1_reverse=codon1_reverse,
            codon2_length=codon2_length, codon2_reverse=codon2_reverse
        )
    except DECLException as e:
        r.error_exit(e)

    """try:
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

    #r.save()"""