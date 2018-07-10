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

        for codon in result["codons"]:
            fh.write("{0}\t{1}\n".format("-".join(codon), result["codons"][codon]))

    # Show log information
    print("{0:<30} {1}".format("Reads processed", result["reads_processed"]))
    print("{0:<30} {1}".format("Reads useful", result["reads_useful"]))

    if r2 is None:
        print("{0:<30} {1}".format("Low quality skips", result["low_quality_skips"]))
    else:
        print("{0:<30} {1}".format("Valid pairs", result["valid_pairs"]))
        print("{0:<30} {1}".format("Invalid pairs", result["invalid_pairs"]))
        print("{0:<30} {1}".format("Low quality skips (both)", result["both_low_quality_skips"]))
        print("{0:<30} {1}".format("Low quality skips (r1)", result["r1_low_quality_skips"]))
        print("{0:<30} {1}".format("Low quality skips (r2)", result["r2_low_quality_skips"]))

    print("\nAll results have been saved into {}".format(result_file))

    #print(r.storage.library.get_formatted_stub_dna_template())