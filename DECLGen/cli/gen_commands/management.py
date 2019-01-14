import sys
import argh
from DECLGen import Runtime
from DECLGen.exceptions import LibraryExistsError, LibraryNotInitializedError


@argh.arg("--advanced-anchors", default=False)
def init(
    template: "A Smiles-String describing the chemical structure of the template. Anchor points for elements must be added as [Rn], where n is a number between 1 and 9",
    advanced_anchors: "Initializes the library in the advanced anchor mode." = False,
) -> None:
    """Initializes a decl-gen datafile in the current directory."""

    try:
        r = Runtime.create(template=template, advanced_anchors=advanced_anchors)
    except LibraryExistsError as e:
        print(
            "Cannot initialize a directory that already has a library init file." + \
            "Use {exec} clear to remove data only, or {exec} remove to remove everything".format(
                exec=sys.argv[0]
            ),
            file=sys.stderr
        )
        exit(e.exitcode)
        return

    try:
        r.save()
    except Exception as e:
        print(e)


def remove() -> None:
    """ Removes the decl-gen datafile from the current directory. """
    try:
        Runtime.remove_init_file()
    except LibraryNotInitializedError as e:
        print("Cannot remove datafile since it does not exist.", file=sys.stderr)
        exit(e.exitcode)
        return

