import multiprocessing as mp
from typing import List
from csv import writer
from DECLGen import Runtime
from DECLGen.molecule import Molecule
from DECLGen.exceptions import DECLException


def lib_info():
    """ Shows information about the library """
    r = Runtime()

    description = r.storage.library.describe()
    for key in description:
        print("{t.bold}{key:<15}{t.normal} {value}".format(t=r.t, key=key + ":", value=description[key]))


def lib_edit(
    template: "New Template" = None,
):
    """ Allows editing of common library features, such as template. """
    r = Runtime()

    try:
        if template is not None:
            r.storage.library.change_template(template)
    except DECLException as e:
        r.error_exit(e)

    r.save()


def iterate_queue(library, queue):
    for item in queue:
        yield [library, item]


def process_molecules(args) -> List[Molecule]:
    library, item = args
    cat1, cats = item

    molecules = []
    for element_set in yield_helper(cat1, cats):
        smiles = library.get_molecule_smiles_by_index(element_set)
        molecule = Molecule(smiles)

        molecules.append(molecule)

    return molecules


def yield_helper(cat1, cats):
    a = 1
    ids = []
    sizes = []
    for cat in cats:
        a *= cat[0]
        ids.append(cat[1])
        sizes.append(cat[0])

    result = {cat1[0]: cat1[1]}
    for i in range(a):
        parts = {}
        for j in range(len(sizes)):
            size = sizes[j]
            id = ids[j]
            k = i % size
            i = i // size
            parts = {**parts, **{id: k}}
        elements = {**result, **parts}
        yield elements


def lib_generate():
    """ Generates physicochemical properties of the library and saves them in library-properties.tsv"""
    r = Runtime()

    queue = r.storage.library.generate_molecule_queue()

    with open("library-properties.csv", "w") as fh:
        csv_file = writer(fh)

        i = 0
        j = 0
        with mp.Pool(8) as pool:
            for molecules in pool.imap_unordered(process_molecules, iterate_queue(r.storage.library, queue)):
                i += len(molecules)
                j += 1

    print("Number of jobs: ", j)
    print("Number of molecules generated: ", i)


