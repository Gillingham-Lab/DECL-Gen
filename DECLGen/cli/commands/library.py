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
    dna_template: "DNA Template. Use {catId} as placeholders to specify where the codon of the category should get included." = None,
):
    """ Allows editing of common library features, such as template. """
    r = Runtime()

    try:
        if template is not None:
            r.storage.library.change_template(template)

        if dna_template is not None:
            r.storage.library.change_dna_template(dna_template)
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

    dna_template = library.get_dna_template()

    molecules = []
    for element_set in yield_helper(cat1, cats):
        keys = list(element_set.keys())
        smiles, codons = library.get_molecule_data_by_index(element_set)
        molecule = Molecule(smiles)

        dna_parts = {keys[x]: codons[x] for x in range(len(keys))}

        try:
            if dna_template is not None:
                dna = library.get_formatted_dna_template(dna_parts)
            else:
                dna = ""
        except KeyError:
            raise Exception("Generation of DNA strand was not possible; Please check if all {catId} in the dna strand exist within the library.")

        molecules.append(codons + [dna] + molecule.get_data(data_fields))

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


def lib_generate(
    threads: "Number of threads" = 1,
    all: "Include all possible data" = False,
    mw: "Include molecular weight" = False,
    qed: "Quantitative estimation of drug-like properties" = False,
    tpsa: "Topological polar surface area" = False,
    tpsapermw: "Topological polar surface area per Da" = False,
    alogp: "aLogP" = False,
    hdonors: "Number of H donors" = False,
    hacceptors: "Number of H acceptors" = False,
    nhetero: "Number of hetero atoms" = False,
    rotatable: "Number of rotatable bonds" = False,
    no: "Number of N and O" = False,
    nhoh: "Number of NH and OH" = False,
    rings: "Number of rings" = False,
    maxRingSize: "Maximum ring size" = False,
    csp3: "Fraction of C atoms that are sp3" = False,
    heavyAtoms: "Number of heavy (non-hydrogen) atoms" = False,
):
    """ Generates physicochemical properties of the library and saves them in library-properties.tsv"""
    r = Runtime()

    queue, elements = r.storage.library.generate_molecule_queue()

    data_fields = {
        "canonical_smiles": True,
        "mw": mw if all is False else True,
        "qed": qed if all is False else True,
        "tpsa": tpsa if all is False else True,
        "tpsapermw": tpsapermw if all is False else True,
        "alogp": alogp if all is False else True,
        "hdonors": hdonors if all is False else True,
        "hacceptors": hacceptors if all is False else True,
        "nhetero": nhetero if all is False else True,
        "rotatable": rotatable if all is False else True,
        "no": no if all is False else True,
        "nhoh": nhoh if all is False else True,
        "rings": rings if all is False else True,
        "maxRingSize": maxRingSize if all is False else True,
        "csp3": csp3 if all is False else True,
        "heavyAtoms": heavyAtoms if all is False else True,
    }

    with open("library-properties.csv", "w") as fh:
        csv_file = writer(fh)
        csv_file.writerow(elements + ["DNA"] + Molecule.get_data_headers(data_fields))

        i = 0
        j = 0
        with mp.Pool(threads) as pool:
            try:
                for molecules in pool.imap_unordered(process_molecules, iterate_queue(r.storage.library, data_fields, queue)):
                    i += len(molecules)
                    j += 1
                    for molecule in molecules:
                        csv_file.writerow(molecule)
            except Exception as e:
                print("{t.red}Error!\n{e}{t.normal}\n".format(e=e,t=r.t))

    print("Number of jobs: ", j)
    print("Number of molecules generated: ", i)


