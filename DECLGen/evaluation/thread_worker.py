from typing import Generator, Optional
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from .metadata import ReadfileWorkerMetadata
from .result import AlignmentResult

class ReadBlock():
    """
    Associates a block of reads and the corresponding meta data
    """
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

    def append(self, r1: SeqRecord, r2: Optional[SeqRecord]):
        """
        Appends a read pair to the block. r2 can be none.
        :param r1:
        :param r2:
        :return:
        """
        self.reads.append((r1, r2))


def read_loader(
        data: ReadfileWorkerMetadata
) -> Generator[ReadBlock, None, None]:
    """
    Iterates over the reads in the files given in data and yields the ReadBlock as soon as the blocksize is reached.
    :param data:
    :return:
    """
    # Open files with gzip if they are zipped
    if data.r1.filename.endswith(".gz"):
        handle_1 = gzip.open(data.r1.filename, "rt", encoding="ASCII")
    else:
        handle_1 = open(data.r1.filename, "r")

    reads_1 = SeqIO.parse(handle_1, "fastq")

    if data.is_paired():
        if data.r2.filename.endswith(".gz"):
            handle_2 = gzip.open(data.r2.filename, "rt", encoding="ASCII")
        else:
            handle_2 = open(data.r2.filename)

        reads_2 = SeqIO.parse(handle_2, "fastq") if data.is_paired() else None
    else:
        reads_2 = None

    read_block = ReadBlock(data)
    read_total = 0

    for read_1 in reads_1:
        read_2 = next(reads_2) if reads_2 is not None else None

        read_block.append(read_1, read_2)

        # If we have enough reads, we yield the package to submit it to a thread.
        if len(read_block) == data.blocksize:
            yield read_block
            read_block = ReadBlock(data)

        max_reads = data.get_kwarg("max_reads", None)
        read_total += 1
        if max_reads is not None and read_total >= max_reads:
            break

    # We yield everything now that's left over
    yield read_block

def read_processor(block: ReadBlock):
    """
    Processes a block of reads in a thread.
    :param block:
    :return:
    """
    paired = block.metadata.is_paired()
    metadata = block.metadata

    skip_codon_matching = False if "skip_codon_matching" not in metadata.kwargs else metadata.kwargs["skip_codon_matching"]

    save_failed = False
    if "save_failed" in metadata.kwargs:
        save_failed = metadata.kwargs["save_failed"]

    result = AlignmentResult(paired=paired)

    for read_1, read_2 in block:
        result["reads_processed"] += 1

        # Call extraction algorithm
        r1_pass, r2_pass, (codons_1, codons_2) = metadata.method((read_1, read_2), metadata)

        if paired:
            # Check and report which read was wrong
            if not r1_pass and not r2_pass:
                result["both_low_quality_skips"] += 1
                result.add_failed_read(read_1.seq, read_2.seq) if save_failed else None
                continue
            elif not r1_pass:
                result["r1_low_quality_skips"] += 1
                result.add_failed_read(read_1.seq, read_2.seq) if save_failed else None
                continue
            elif not r2_pass:
                result["r2_low_quality_skips"] += 1
                result.add_failed_read(read_1.seq, read_2.seq) if save_failed else None
                continue

            if skip_codon_matching is False:
                # If paired, check if both codons are the same.
                if codons_1 == codons_2 and codons_1 != None:
                    result["valid_pairs"] += 1
                else:
                    result["invalid_pairs"] += 1
                    result.add_failed_read(read_1.seq, read_2.seq) if save_failed else None
                    continue
            # If skip_codon_matching is used, we assume that something is not okay - like codons not covered on one strand.
            else:
                true_codons = []
                for i in range(len(codons_1)):
                    c1 = codons_1[i]
                    c2 = codons_2[i]

                    if len(c1) != len(c2):
                        result["invalid_pairs"] += 1
                        result.add_failed_read(read_1.seq, read_2.seq) if save_failed else None
                        continue

                    nc = ""

                    for y in range(len(c1)):
                        nt1 = c1[y]
                        nt2 = c2[y]

                        if nt1 == "-":
                            nc += nt2
                        elif nt2 == "-":
                            nc += nt1
                        elif nt1 == nt2:
                            nc += nt1
                        else:
                            result["invalid_pairs"] += 1
                            result.add_failed_read(read_1.seq, read_2.seq) if save_failed else None
                            continue

                    true_codons.append(nc)

                codons_1 = true_codons
                result["valid_pairs"] += 1
        else:
            # Report if read was wrong
            if not r1_pass:
                result["low_quality_skips"] += 1
                result.add_failed_read(read_1.seq, read_2.seq) if save_failed else None
                continue

        # Convert
        codons_1 = tuple([str(x) for x in codons_1])
        result.increase_codon(codons_1)

    return result