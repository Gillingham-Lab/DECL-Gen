import sys
import argh
from DECLGen import Runtime
from DECLGen.exceptions import LibraryExistsError, LibraryNotInitializedError


@argh.arg("--enable-advanced-anchors", default=False)
@argh.arg("--enable-superset-categories", default=False)
def init(
    template: "A Smiles-String describing the chemical structure of the template. Anchor points for elements must be added as [Rn], where n is a number between 1 and 9",
    enable_advanced_anchors: "Initializes the library in the advanced anchor mode." = False,
    enable_superset_categories: "EXPERIMENTAL: Allows to use superset categories (ssc)" = False,
) -> None:
    """Initializes a decl-gen datafile in the current directory."""

    try:
        r = Runtime.create(
            template=template,
            advanced_anchors=enable_advanced_anchors,
            superset_categories=enable_superset_categories,
        )
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

