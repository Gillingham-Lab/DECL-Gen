import os
import argh
from timeit import default_timer as timer

from DECLGen import Runtime
from DECLGen.exceptions import LibraryNoDNATemplateException
from DECLGen.evaluation.qc import Type

@argh.arg("r2", nargs="?")
@argh.arg("--method", choices=Type.all)
@argh.arg("--blocksize", type=int)
@argh.arg("--quality", type=float)
@argh.arg("--max-reads", type=int)
@argh.arg("--timing", default=False)
def extract(
        r1: "Fastq file containing forward reads",
        r2: "Fastq file containing reverse reads from paired read, if available",
        threads: "Number of threads" = 8,
        method = Type.default,
        quality: "Quality number whose purpose changes depending on the method used." = None,
        blocksize: "Size of reads to load before sending them to a thread" = 10000,
        result_file: "Filename of the results. Generated automatically if not given" = None,
        max_reads = None,
        timing = False,
        save_failed: "Save failed reads" = False,
):
    r = Runtime()

    s, e = (0, 0)
    if timing:
        s = timer()

    if r.storage.library.get_dna_template() == None:
        a = LibraryNoDNATemplateException("You have not defined a DNA template.")
        r.error_exit(a)

    result = r.storage.library.evaluate_sequencing_results(
        r1,
        r2,
        threads=threads,
        blocksize=blocksize,
        method=method,
        quality=quality,
        max_reads=max_reads
    )

    # Save result
    if result_file is None:
        result_file = "".join([os.path.basename(r1).split(".")[0], ".result.csv"])

    with open(result_file, "w") as fh:
        fh.write("{0}\t{1}\n".format("Codon-Combination", "Count"))

        codon_list = result.get_codons()
        for codon in codon_list:
            fh.write("{0}\t{1}\n".format("-".join([str(x) for x in codon]), codon_list[codon]))

    with open(result_file[:-4] + ".log", "w") as fh:
        fh.write(str(result))

    if save_failed:
        with open("".join([os.path.basename(r1).split(".")[0], ".1.failed"]), "w") as fr1, \
            open("".join([os.path.basename(r1).split(".")[0], ".2.failed"]), "w") as fr2:
            for read1, read2 in result._failed_reads:
                fr1.write(read1 + "\n")
                fr2.write(read2 + "\n")

    if timing:
        e = timer()

    # Show log information
    print(result)
    print("\nAll results have been saved into {}".format(result_file))

    if timing:
        print("Time required: {:.2f}".format(e-s))