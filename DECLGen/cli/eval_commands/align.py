import os
import argh

from DECLGen import Runtime
from DECLGen.exceptions import LibraryNoDNATemplateException

@argh.arg("r2", nargs="?")
def align(
    r1: "Fastq file containing forward reads",
    r2: "Fastq file containing reverse reads from paired read, if available",
    threads: "Number of threads" = 8,
    blocksize: "Size of reads to load before sending them to a thread" = 10000,
    checktype: "Check type" = "simple",
    result_file: "Filename of the results. Generated automatically if not given" = None,
):
    r = Runtime()

    if r.storage.library.get_dna_template() == None:
        a = LibraryNoDNATemplateException("You have not defined a DNA template.")
        r.error_exit(a)

    result = r.storage.library.align(r1, r2, threads=threads, blocksize=blocksize, checktype=checktype)

    # Save result
    if result_file is None:
        result_file = os.path.basename(r1).split(".")[0] + ".result.csv"

    with open(result_file, "w") as fh:
        fh.write("{0}\t{1}\n".format("Codon-Combination", "Count"))

        codon_list = result.get_codons()
        for codon in codon_list:
            fh.write("{0}\t{1}\n".format("-".join(codon), codon_list[codon]))

    # Show log information
    print(result)

    print("\nAll results have been saved into {}".format(result_file))

    #print(r.storage.library.get_formatted_stub_dna_template())