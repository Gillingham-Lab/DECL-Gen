import multiprocessing as mp
from typing import List, Dict, Iterable
from csv import writer
from DECLGen import Runtime
from DECLGen.molecule import Molecule
from DECLGen.exceptions import DECLException
from DECLGen.library import Library


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


def iterate_queue(
    library: Library,
    data_fields: Dict[str, bool],
    queue: Iterable
) -> List:
    for item in queue:
        yield [library, data_fields, item]


def process_molecules(args: List) -> List:
    library, data_fields, item = args
    cat1, cats = item

    molecules = []
    for element_set in yield_helper(cat1, cats):
        smiles, codons = library.get_molecule_data_by_index(element_set)
        molecule = Molecule(smiles)

        molecules.append(codons + molecule.get_data(data_fields))

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

    queue, elements = r.storage.library.generate_molecule_queue()

    data_fields = {
        "canonical_smiles": True
    }

    with open("library-properties.csv", "w") as fh:
        csv_file = writer(fh)
        csv_file.writerow(elements + Molecule.get_data_headers(data_fields))

        i = 0
        j = 0
        with mp.Pool(8) as pool:
            for molecules in pool.imap_unordered(process_molecules, iterate_queue(r.storage.library, data_fields, queue)):
                i += len(molecules)
                j += 1
                for molecule in molecules:
                    csv_file.writerow(molecule)

    print("Number of jobs: ", j)
    print("Number of molecules generated: ", i)


