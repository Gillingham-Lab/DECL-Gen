from typing import List, Tuple, Union
import re
import sys
from rdkit import Chem
from Bio.Seq import Seq

_search_pattern = r"\[R([0-9]+)\]"


def sanitize(raw_template: str) -> str:
    """ This method tries to sanitize a raw_template by moving the R-group to a better position.

    R-groups are essentially placeholders for SMILES connections commonly used to form rings. As such, they must
    always be placed after an atom, or else the connection would not work.

    This method tries to correct these cases:
        sanitize("[R1]CCC") -> "C[R1]CC"
    """
    if raw_template.find(".") >= 0:
        template_parts = raw_template.split(".")
        template_parts_sanitized = []
        for part in template_parts:
            part = sanitize(part)
            template_parts_sanitized.append(part)

        raw_template = ".".join(template_parts_sanitized)
    else:
        if raw_template.startswith("[R"):
            atom_pos = None
            atom_length = None

            for atom in raw_template:
                if atom in "CNO":
                    atom_pos = raw_template.find(atom)

                    if raw_template[atom_pos + 1] in "0123456789":
                        atom_length = 2
                    else:
                        atom_length = 1

                    break

            raw_template = "{}({}){}".format(
                raw_template[atom_pos:atom_pos + atom_length],
                raw_template[0:atom_pos],
                raw_template[atom_pos + atom_length:]
            )

        raw_template = re.sub(r"(\([A-Z0-9\(\)\[\]]+\))(\(\[R[0-9]+\]\))", lambda x: "".join([x.group(2), x.group(1)]), raw_template)
        raw_template = re.sub(r"\((\[R([0-9]+)\])\)", lambda x: x.group(1), raw_template)

    return raw_template


def count_anchors(raw_template: str) -> int:
    return len(re.findall(_search_pattern, raw_template))


def get_anchors(raw_template: str) -> List[str]:
    """ Returns all anchors defined in the scaffold. """
    return sorted(["R{}".format(x) for x in re.findall(_search_pattern, raw_template)])


def _parse_callback(n) -> str:
    n = int(n.group(1))
    if n < 90:
        return "%{}".format(99 - n)
    else:
        raise ValueError("{name} only supports R labels from R0 to R89 due to technical reasons.".format(name=sys.argv[0]))


def parse(raw_template: str) -> str:
    return re.sub(_search_pattern, _parse_callback, raw_template)


def is_valid(parsed: str, anchors: List[str]) -> bool:
    test_smiles_fragments = [parse("C[{}]".format(x)) for x in anchors]
    fragments = [parsed] + test_smiles_fragments
    smiles = ".".join(fragments)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    else:
        return True


def get_codon_coordinates(template: Union[Seq, str]) -> List[Tuple[int, int]]:
    """
    Retrieves codon positions from a given template as a list of tuples with from-to coordinates. Codons must be
    indicated as N.
    :param template: The template to derive the coordinates from. A Bio.Seq.Seq object.
    :return: List of coordinate tuples (from, to), 0-indexed.
    """
    start = 0
    positions = []
    while True:
        temp_start = template.find("N", start)
        temp_end = temp_start
        for char in template[temp_start:]:
            if char != "N":
                break

            temp_end+=1

        start = temp_end
        positions.append((temp_start, temp_end))

        if template.find("N", start) < 0:
            break

    return positions

