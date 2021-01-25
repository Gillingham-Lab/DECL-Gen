import os
import argh
from timeit import default_timer as timer

from DECLGen import Runtime
from DECLGen.exceptions import LibraryNoDNATemplateException
from DECLGen.evaluation.qc import Type
from DECLGen.cli.helpers import ProgressBar


@argh.arg("r2", nargs="?")
@argh.arg("--method", choices=Type.all)
@argh.arg("--blocksize", type=int)
@argh.arg("--quality", type=float)
@argh.arg("--max-reads", type=int)
@argh.arg("--timing", default=False)
@argh.arg("--no-auto-detection", default=False)
def extract(
        r1: "Fastq file containing forward reads",
        r2: "Fastq file containing reverse reads from paired read, if available",
        threads: "Number of threads" = 8,
        method = Type.default,
        quality: "Quality number whose purpose changes depending on the method used." = None,
        blocksize: "Size of reads to load before sending them to a thread" = None,
        result_file: "Filename of the results. Generated automatically if not given" = None,
        max_reads: "Maximum amount of reads to read in. Good for quick data inspection." = None,
        timing: "Measure the timing" = False,
        save_failed: "Save failed reads" = False,
        no_auto_detection: "Deactivates the autodetection of which file is forward and which is the reverse read." = False,
        skip_codon_matching: "Skips the automated invalidation of mismatching codons." = False,
):
    r = Runtime()

    s, e = (0, 0)
    if timing:
        s = timer()

    if r.storage.library.get_dna_template() == None:
        a = LibraryNoDNATemplateException("You have not defined a DNA template.")
        r.error_exit(a)

    # Create a progress bar that can be updated.
    with ProgressBar(r.t, desc="Reading reads") as progressBar:
        # Evaluates the sequencing results using multiple threads.
        result = r.storage.library.evaluate_sequencing_results(
            r1,
            r2,
            threads=threads,
            blocksize=blocksize,
            method=method,
            quality=quality,
            max_reads=max_reads,
            progressBar=progressBar,
            no_auto_detection=no_auto_detection,
            save_failed=save_failed,
            skip_codon_matching=skip_codon_matching,
        )

    # Save result
    if result_file is None:
        result_file = "".join([os.path.basename(r1).split(".")[0], ".result.csv"])

    # Write the results to a .csv file
    with open(result_file, "w") as fh:
        fh.write("{0}\t{1}\n".format("Codon-Combination", "Count"))

        codon_list = result.get_codons()
        for codon in codon_list:
            fh.write("{0}\t{1}\n".format("-".join([str(x) for x in codon]), codon_list[codon]))

    # Write the log file, including the time required if requested.
    with open(result_file[:-4] + ".log", "w") as fh:
        fh.write(str(result))

        if timing:
            e = timer()
            print("\nTime required: {:.2f}\n".format(e-s))

    # Write the failed-to-use reads in separate files as well. Might be useful for analysis later.
    if save_failed:
        with open("".join([os.path.basename(r1).split(".")[0], "-forward.failed"]), "w") as fr1, \
            open("".join([os.path.basename(r1).split(".")[0], "-reverse.failed"]), "w") as fr2:
            for read1, read2 in result._failed_reads:
                fr1.write(str(read1) + "\n")
                fr2.write(str(read2) + "\n")

    # Take the full timing.
    if timing:
        e = timer()

    # Show log information
    print(result)
    print("\nAll results have been saved into {}".format(result_file))

    if timing:
        print("Time required: {:.2f}".format(e-s))