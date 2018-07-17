from typing import List, Set, Dict, Generator, Tuple, Union, Optional
from operator import mul, itemgetter
from functools import reduce
from Bio import pairwise2, SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import multiprocessing as mp

from DECLGen.exceptions import \
    LibraryCategoryException, \
    LibraryCategoryExistsException, \
    LibraryCategoryNotFoundException, \
    LibraryTemplateException, \
    EvaluationException
from .categories import Category
from DECLGen import template, codon
from DECLGen.evaluation import \
    read_loader, \
    read_processor, \
    qc
from DECLGen.evaluation.metadata import \
    ReadfileWorkerMetadata, \
    ReadfileType, \
    ReadfileMetadata

def _check_anchor(anchor) -> bool:
    if not anchor.startswith("R"):
        raise LibraryCategoryException(
            "R-Group <{}> must start with R.".format(anchor)
        )

    if len(anchor) < 2:
        raise LibraryCategoryException(
            "R-Group <{}> must be Rn or Rnn, whereas n is a number from 0-9.".format(anchor)
        )

    if int(anchor[1:]) < 0:
        raise LibraryCategoryException(
            "R-Group <{}> must be a positive number.".format(anchor)
        )

    return True


class Library:
    storage = None
    anchors = None
    anchors_in_use = None
    categories = None
    dna_template = None

    def __init__(self, storage):
        self.storage = storage
        self.anchors = []
        self.anchors_in_use = {}
        self.categories = {}

    def set_anchors(self, anchors: List[str]):
        self.anchors = sorted(anchors)

    def change_template(self, new_template: str) -> None:
        # We must sanitize the template first in case the R group is at the beginning
        template_string = template.sanitize(new_template)

        raw_template = template_string
        smiles_template = template.parse(template_string)
        anchors = template.get_anchors(template_string)

        if len(self.anchors) > 0 and sorted(anchors) != sorted(self.anchors):
            raise LibraryTemplateException("Cannot change template with different anchors. You must remove the library first.")

        self.storage.raw_template = raw_template
        self.storage.smiles_template = smiles_template

    def change_dna_template(self, dna_template: str) -> None:
        """ Changes the DNA template of the library. {catId} can be used as a placeholders."""
        self.dna_template = dna_template

    def get_dna_template(self) -> Union[str, None]:
        """ Returns the DNA template if set. """
        return self.dna_template

    def get_formatted_stub_dna_template(self):
        """ Returns a formatted DNA template where the codons are replaced by N. """
        formatting = {}
        for catId in self.categories:
            formatting[catId] = "N" * self.categories[catId].get_codon_length()

        return self.dna_template.format(**formatting)

    def get_formatted_dna_template(self, codon_list: Dict[str, str]) -> Union[str, None]:
        for catId in codon_list.keys():
            if self.has_category(catId) is False:
                raise KeyError("Unknown category id ({})".format(catId))

            cat = self.get_category(catId)
            if cat.is_reverse_complement():
                codon_list[catId] = codon.reverse(codon_list[catId])

        return self.dna_template.format(**codon_list)

    def get_codon_summary_string(self, codon_list: Dict[str, str]) -> Union[str, None]:
        for catId in codon_list.keys():
            if self.has_category(catId) is False:
                raise KeyError("Unknown category id ({})".format(catId))

            cat = self.get_category(catId)
            if cat.is_reverse_complement():
                codon_list[catId] = codon.reverse(codon_list[catId])

        template = self.dna_template
        if template is None:
            template = "A" + "A".join(["{" + x + "}" for x in codon_list.keys()]) + "A"

        template = template.replace("}", "{").split("{")
        codon_summary = []
        for stuff in template:
            if stuff in codon_list:
                codon_summary.append(codon_list[stuff])
        return "-".join(codon_summary)


    def describe(self) -> Dict[str, Union[str, None]]:
        description = {
            "Template": self.storage.raw_template,
            "DNA-Template": self.get_dna_template(),
            "R-Groups": ", ".join(self.anchors),
            "Library Shape": ", ".join([str(x) for x in self.shape]),
            "Library Size": self.get_size(),
        }

        return description

    def get_size(self) -> int:
        return reduce(mul, self.shape, 1)

    def has_category(self, id: str):
        if id in self.categories:
            return True
        else:
            return False

    def get_categories(self) -> List[Category]:
        return list(self.categories.values())

    def get_category(self, id) -> Category:
        if not self.has_category(id):
            raise LibraryCategoryNotFoundException("Library category <{}> not found".format(id))

        return self.categories[id]

    @property
    def shape(self) -> Tuple[int]:
        shape = []
        for cat in self.categories.values():
            shape.append(len(cat))
        return tuple(shape)

    def add_category(self, id: str, name: str, anchors: List[str], codon_length: int = 0) -> bool:
        """
        Adds a new diversity element category to the library.
        :param id: The id of the category
        :param name: The name of the category
        :param anchors: A list of R-Group anchors (R1, R2 ...)
        :param codon_length: Codon length
        :return:
        """
        # Check if id is already in use
        if self.has_category(id):
            raise LibraryCategoryExistsException("A category with the id <{}> already exists".format(id))

        # Check validity of each anchor and if anchor is free
        if len(anchors) == 0:
            raise LibraryCategoryException("You must at least give 1 anchor, none were given.")

        anchors_used = []
        for anchor in anchors:
            _check_anchor(anchor)

            if anchor in self.anchors_in_use:
                anchors_used.append(anchor)

        if len(anchors_used) > 0:
            raise LibraryCategoryException(
                "Anchors <{}> are already in use.".format(", ".join(anchors_used))
            )

        # Create the new category
        self.categories[id] = Category(id, name, anchors, codon_length)
        # Marks anchors as used by cross-referencing Category
        for anchor in anchors:
            self.anchors_in_use[anchor] = self.categories[id]

        return True

    def del_category(self, id: str) -> bool:
        # Check if id is already in use
        if not self.has_category(id):
            raise LibraryCategoryNotFoundException("A category with the id <{}> does not exist exists".format(id))

        # Remove used anchors
        anchors = self.categories[id].get_anchors()
        for anchor in anchors:
            del self.anchors_in_use[anchor]

        # Remove category
        del self.categories[id]

        return True

    def get_molecule_data_by_index(self, elements: Dict[str, int]) -> str:
        fragments = []
        codons = []
        for cat_id in elements:
            fragments.append(self.categories[cat_id]._get_element_by_list_index(elements[cat_id]).parsed_smiles)
            codons.append(self.categories[cat_id]._get_element_by_list_index(elements[cat_id]).codon)

        fragments.append(self.storage.smiles_template)
        smiles = ".".join(fragments)

        return smiles, codons

    def generate_molecule_queue(self) -> Tuple[list, list]:
        categories_by_size = []
        for cat_id in self.categories:
            categories_by_size.append((len(self.categories[cat_id]), cat_id))

        categories_by_size.sort()

        queue = []
        for element1 in range(categories_by_size[0][0]):
            gen = ([categories_by_size[0][1], element1], categories_by_size[1:])
            queue.append(gen)

        return queue, [x[1] for x in categories_by_size]

    def evaluate_sequencing_results(self,
                                    r1: str,
                                    r2: Optional[str] = None,
                                    compare_n: int = 5,
                                    threads: int = 8,
                                    blocksize: int = 10,
                                    method: str = "simple",
                                    quality: Optional[Union[int, float]] = None,
                                    **kwargs
                                    ) -> None:
        """
        Evaluates given sequencing result files
        :param r1: fastq file containing the reads
        :param r2: optional fastq file containing the paired reads
        :param compare_n: A comparison factore. @ToDo: refactor!
        :param threads: Number of threads to be used
        :param blocksize: If harddrive is fast, higher values are better
        :param method: Evaluation method. A choice of qc.all.
        :return:
        """
        template_f = Seq(self.get_formatted_stub_dna_template().upper(), alphabet=IUPAC.ambiguous_dna)
        template_r = template_f.reverse_complement()

        r1_template, r2_template = self._assign_template(template_f, template_r, r1, r2, n=compare_n)

        n_positions_r1 = template.get_codon_coordinates(r1_template)
        n_positions_r2 = template.get_codon_coordinates(r2_template)

        if r1_template == template_f:
            r1_results = ReadfileType.forward
            r2_results = ReadfileType.reverse
        else:
            r1_results = ReadfileType.reverse
            r2_results = ReadfileType.forward

        r1 = ReadfileMetadata(r1, r1_template, n_positions_r1, r1_results)
        r2 = ReadfileMetadata(r2, r2_template, n_positions_r2, r2_results) if r2 is not None else None

        all_results = None

        # Sanitize checktypes
        if isinstance(method, str):
            if method in qc.Type.all:
                method = qc.Type.get(method)
            else:
                raise EvaluationException("Unknown method given ({}). Use one of {}.".format(method, ", ".join(qc.Type.all)))

        # Initialize Worker metadata
        worker_metadata = ReadfileWorkerMetadata(r1, r2,
                                                 blocksize=blocksize, method=method, compare_n=compare_n,
                                                 quality=quality, **kwargs)

        with mp.Pool(threads) as pool:
            for temp_result in pool.imap_unordered(read_processor, read_loader(worker_metadata)):
                # Do something with the result
                if all_results is None:
                    all_results = temp_result
                else:
                    all_results += temp_result

        return all_results

    def _assign_template(self, template_f: Seq, template_r: Seq, r1: str, r2: str, n: int = 5):
        reads_1 = SeqIO.parse(r1, "fastq-illumina")
        reads_2 = SeqIO.parse(r2, "fastq-illumina") if r2 is not None else None

        reads_1_template = None
        reads_2_template = None
        max_trials = 10

        for read_1 in reads_1:
            read_2 = next(reads_2) if reads_2 is not None else None

            # Try first to figure out which template fits the read better
            if read_1.seq[0:n] == template_f[0:n] and \
                    (reads_2 is None or read_2.seq[0:n] == template_r[0:n]):
                reads_1_template = template_f
                reads_2_template = template_r
                break
            elif read_1.seq[0:n] == template_r[0:n] and \
                    (reads_2 is None or read_2.seq[0:n] == template_f[0:n]):
                reads_1_template = template_r
                reads_2_template = template_f
                break
            else:
                max_trials -= 1

                if max_trials <= 0:
                    raise EvaluationException("Unable to figure out to which file is forward and reverse read. "
                                              "Are you sure the read files fit the template?")
                continue

        return reads_1_template, reads_2_template


