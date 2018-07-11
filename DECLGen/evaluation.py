from typing import List, Tuple, Union, Generator
from enum import Enum
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


class ReadfileType(Enum):
    forward = 0
    reverse = 1


class QualityControlType(Enum):
    Simple = 0
    Score = 1
    TrueAlign = 2


class ReadfileMetadata():
    filename = None
    template = None
    coordinates = None
    type = None

    def __init__(self, filename: str, template: Seq, coordinates: List[Tuple[int, int]], type=ReadfileType):
        self.filename = filename
        self.template = template
        self.coordinates = coordinates
        self.type = type

    def is_reverse(self):
        return self.type == ReadfileType.reverse

    def is_forward(self):
        return self.type == ReadfileType.forward

    def is_none(self):
        if self.filename is None:
            return True

        return False


class ReadfileWorkerMetadata():
    r1 = None
    r2 = None
    compare_n = None
    blocksize = None
    checktype = None

    def __init__(self,
                 r1: ReadfileMetadata,
                 r2: Union[ReadfileMetadata, None],
                 compare_n: int,
                 blocksize: int,
                 checktype
                 ):
        self.r1 = r1
        self.r2 = r2
        self.compare_n = compare_n
        self.blocksize = blocksize
        self.checktype = checktype

    def is_single(self):
        return True if self.r2 is None else False

    def is_paired(self):
        return False if self.r2 is None else True


class ReadBlock():
    reads = None
    metadata = None

    def __init__(self, metadata: ReadfileWorkerMetadata):
        self.reads = []
        self.metadata = metadata

    def __len__(self):
        return len(self.reads)

    def __iter__(self):
        for read in self.reads:
            yield read

    def append(self, r1, r2):
        self.reads.append((r1, r2))


class AlignmentResult():
    """ Represents the result of codon alignment. """
    _result = None
    _key_annotations = {
        "reads_processed": ("Processed Reads", "both"),
        "reads_useful": ("Useful Reads", "both"),
        "valid_pairs": ("Valid Pairs", "paired"),
        "invalid_pairs": ("Invalid Pairs", "paired"),
        "low_quality_skips": ("Low Quality skips", "single"),
        "both_low_quality_skips": ("Low Quality skips (both)", "paired"),
        "r1_low_quality_skips": ("Low Quality skips (r1)", "paired"),
        "r2_low_quality_skips": ("Low Quality skips (r2)", "paired"),
    }
    paired = False
    _codons = None

    def __init__(self, paired=True):
        self._result = {key: 0 for key in self._key_annotations}
        self._codons = {}
        self._paired = paired

    def __getitem__(self, item):
        if item not in self._result:
            raise ValueError(
                "Result only supports a set amount of keys: {}, but {} given".format(
                    ", ".join(self._result.keys()),
                    item
                )
            )

        return self._result[item]

    def __setitem__(self, item, value):
        if item not in self._result:
            raise ValueError(
                "Result only supports a set amount of keys: {}, but {} given".format(
                    ", ".join(self._result.keys()),
                    item
                )
            )

        self._result[item] = value

    def __add__(self, othr):
        r = self.__class__()
        # Add result meta
        for key in self._result:
            r[key] = self[key] + othr[key]

        # Merge codon lists
        self_codons = self.get_codons()
        othr_codons = othr.get_codons()

        for codon in self_codons:
            r.increase_codon(codon, self_codons[codon])
        for codon in othr_codons:
            r.increase_codon(codon, othr_codons[codon])

        return r

    def __str__(self):
        ret = []

        for key in self._result:
            annotation = self._key_annotations[key]
            value = self._result[key]

            if self._paired is True and annotation[1] in ["both", "paired"]:
                ret.append("{0:<30} {1}".format(annotation[0], value))
            elif self._paired is False and annotation[1] in ["both", "single"]:
                ret.append("{0:<30} {1}".format(annotation[0], value))

        return "\n".join(ret)

    def init_codon(self, codon) -> None:
        if codon not in self._codons:
            self._codons[codon] = 0

    def increase_codon(self, codon, increase=1):
        self.init_codon(codon)
        self._codons[codon] += increase

    def get_codons(self):
        return self._codons

def extract_codon(read: Seq, coordinates: List[Tuple[int, int]], invert=False) -> List[str]:
    """
    Extracts codons from a given Bio.Seq.Seq object depending on the coordinates given.

    :param read: The read to extract the codons from.
    :param coordinates: A list of coordinate-tuples for all codons.
    :param invert: If true, returns the coordinates as their reverse complement.
    :return: A list of codons as strings.
    """
    codons = []
    for start, stop in coordinates:
        codon = read[start:stop]

        if invert is True:
            codons.append(codon.reverse_complement())
        else:
            codons.append(codon)

    return codons if invert is False else codons[::-1]


def read_loader(
        data: ReadfileWorkerMetadata
) -> Generator[ReadBlock, None, None]:
    reads_1 = SeqIO.parse(data.r1.filename, "fastq-illumina")
    reads_2 = SeqIO.parse(data.r2.filename, "fastq-illumina") if data.is_paired() else None

    read_block = ReadBlock(data)

    for read_1 in reads_1:
        read_2 = next(reads_2) if reads_2 is not None else None

        read_block.append(read_1, read_2)

        # If we have enough reads, we yield the package to submit it to a thread.
        if len(read_block) == data.blocksize:
            yield read_block
            read_buffer = ReadBlock(data)

    # We yield everything now that's left over
    yield read_block


def compare_sequence(seq1: Seq, seq2: Seq, start: int = 0, stop: Union[int, None] = None):
    return seq1[start:stop] == seq2[start:stop]


def _qc_simple_per_read(read: Seq, r: ReadfileMetadata, n: int) -> bool:
    has_passed = True

    # Check beginning of read
    if compare_sequence(read, r.template, 0, n) is False:
        has_passed = False

    # Check around codon +- N for codon length L
    for positions in r.coordinates:
        L = positions[1] - positions[0]

        # Check upstream of codon
        if compare_sequence(read, r.template, positions[0] - L, positions[1] - L) is False:
            has_passed = False

        # Check downstream of codon
        if compare_sequence(read, r.template, positions[0] + L, positions[1] + L) is False:
            has_passed = False

    return has_passed


def qc_simple(
        reads: Tuple[SeqRecord, SeqRecord],
        metadata: ReadfileWorkerMetadata
) -> Tuple[bool, bool]:
    read_1, read_2 = reads
    r1 = metadata.r1
    r2 = metadata.r2

    r1_pass = _qc_simple_per_read(read_1.seq, r1, metadata.compare_n)
    r2_pass = True

    if metadata.is_paired():
        r2_pass = _qc_simple_per_read(read_2.seq, r2, metadata.compare_n)

    return r1_pass, r2_pass


def read_processor(block: ReadBlock):
    #read_buffer, r1, r2, compare_n, blocksize, checktype = args
    #r1, r1_template, n_positions_r1, r1_results = r1
    #r2, r2_template, n_positions_r2, r2_results = r2

    paired = block.metadata.is_paired()
    metadata = block.metadata
    r1 = block.metadata.r1
    r2 = block.metadata.r2

    result = AlignmentResult(paired=paired)

    for read_1, read_2 in block:
        result["reads_processed"] += 1

        # Extract codons from first read
        codons_1 = extract_codon(read_1.seq, r1.coordinates, r1.is_reverse())
        codons_2 = None

        if paired:
            codons_2 = extract_codon(read_2.seq, r2.coordinates, r2.is_reverse())

            # Check validity of the pairs
            if codons_1 == codons_2:
                result["valid_pairs"] += 1
            else:
                result["invalid_pairs"] += 1
                continue

        if block.metadata.checktype == "simple":
            r1_pass, r2_pass = qc_simple((read_1, read_2), metadata)

            if paired:
                if not r1_pass and not r2_pass:
                    result["both_low_quality_skips"] += 1
                    continue
                elif not r1_pass:
                    result["r1_low_quality_skips"] += 1
                    continue
                elif not r2_pass:
                    result["r2_low_quality_skips"] += 1
                    continue
            else:
                if not r1_pass:
                    result["low_quality_skips"] += 1
                    continue

        elif block.metadata.checktype == "pairwise_score":
            # Alignment (potentially slower down)
            score_1 = pairwise2.align.localms(read_1.seq, r1_template, 5, -3, -4, -4, score_only=True)
            max_score_1 = min(len(read_1.seq), len(r1_template)) * 5 * 0.5

            if read_2 is not None:
                score_2 = pairwise2.align.localms(read_2.seq, r2_template, 5, -3, -4, -4,
                                                  score_only=True)
                max_score_2 = min(len(read_2.seq), len(r2_template)) * 5 * 0.5

            if read_2 is not None:
                if score_1 < max_score_1 and score_2 < max_score_2:
                    result["both_low_quality_skips"] += 1
                    continue
                elif score_1 < max_score_1:
                    result["r1_low_quality_skips"] += 1
                    continue
                elif score_2 < max_score_2:
                    result["r2_low_quality_skips"] += 1
                    continue
            else:
                if score_1 < max_score_1:
                    result["low_quality_skips"] += 1
                    continue
        elif block.metadata.checktype == "pairwise_complex":
            n1 = min(len(read_1.seq), len(r1_template))
            alignment_1 = pairwise2.align.localms(read_1.seq[:n1], r1_template[:n1], 5, -3, -4, -4, one_alignment_only=True)

            if read_2 is not None:
                n2 = min(len(read_1.seq), len(r2_template))
                alignment_2 = pairwise2.align.localms(read_2.seq[:n2], r2_template[:n2], 5, -3, -4, -4, one_alignment_only=True)

            if read_2 is not None:
                if alignment_1[0].count("-") > 0 and alignment_2[0].count("-"):
                    result["both_low_quality_skips"] += 1
                    continue
                elif alignment_1[0].count("-") > 0:
                    result["r1_low_quality_skips"] += 1
                    continue
                elif alignment_2[0].count("-"):
                    result["r2_low_quality_skips"] += 1
                    continue
            elif alignment_1[0].count("-") > 0:
                result["low_quality_skips"] += 1
                continue


        result["reads_useful"] += 1
        codons_1 = tuple([str(x) for x in codons_1])
        result.increase_codon(codons_1)

    return result