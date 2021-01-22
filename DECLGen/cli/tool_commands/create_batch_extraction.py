import os
import argh
from colorama import init, Fore, Style

from DECLGen.evaluation.qc import Type

@argh.arg("--method", choices=Type.all)
@argh.arg("--blocksize", type=int)
@argh.arg("--quality", type=float)
@argh.arg("--max-reads", type=int)
@argh.arg("--timing", default=False)
def create_batch_extraction(
    path: "Path to generate files from",
    batch_name: "Filename for the batch file",
    no_mates: "Set to true if you don't have mate reads'" = False,
    invert_mates: "Set to true if mates should be inverted" = False,
    threads: "Number of threads" = 4,
    method: "Extraction algorithm. Simple is quickest." = "simple",
    quality: "Quality number whose purpose changes depending on the method used." = 2,
    blocksize: "Size of reads to load before sending them to a thread" = None,
    result_file: "Filename of the results. Generated automatically if not given" = None,
    max_reads: "Maximum amount of reads to read in. Good for quick data inspection." = None,
    timing: "Measure the timing" = False,
    save_failed: "Save failed reads" = False,
    skip_codon_matching: "Skips the automated invalidation of mismatching codons." = False,
):
    """
    Generates a batch file to extract reads from all files within a given path.

    The generated batch file can be used as a seed file for a slurm array job.

    Should be called from the library working directory to get the path right.
    """
    init()

    # Error out if the path does not exist
    if not os.path.exists(path):
        print(f"{Fore.RED}The given path was not found ({path}).{Fore.RESET}")
        return 1

    # Error out of the path is not a directory
    if not os.path.isdir(path):
        print(f"{Fore.RED}The given path is not a directory ({path}).{Fore.RESET}")
        return 1

    # Error out if not in a library directory
    if not os.path.exists("decl_gen.data"):
        print(f"{Fore.YELLOW}This command should be called from a library folder (no decl_gen.data found).{Fore.RESET}")

    # List everything from current directory
    files = os.listdir(path)
    # Prepare the two mate lists
    mates1 = []
    mates2 = []

    for file in files:
        file = os.path.join(path, file)

        if os.path.isdir(file):
            continue

        if not file.endswith(".gz") and not file.endswith(".fq"):
            continue

        if no_mates:
            mates1.append(file)
        else:
            if file.endswith("1.fq.gz") or file.endswith("1.fq"):
                if invert_mates:
                    mates2.append(file)
                else:
                    mates1.append(file)
            elif file.endswith("2.fq.gz") or file.endswith("2.fq"):
                if invert_mates:
                    mates1.append(file)
                else:
                    mates2.append(file)

    mates1 = sorted(mates1)
    mates2 = sorted(mates2)

    if len(mates1) != len(mates2):
        print(f"{Fore.RED}The length of detected mates1 and mates2 don't match.{Fore.RESET}")
        return 1

    if no_mates:
        pairs = mates1
    else:
        pairs = zip(mates1, mates2)

    with open(batch_name, "w") as fw:
        for pair in pairs:
            if no_mates:
                row = [f"declEval extract {pair}"]
            else:
                row = [f"declEval extract {pair[0]} {pair[1]}"]

            if result_file:
                row.append(f"--result-file {result_file}")

            if threads:
                row.append(f"--threads 4")

            row.append(f"--method {method}")
            row.append(f"--quality {quality}")

            if blocksize:
                row.append(f"--blocksize {blocksize}")

            if timing:
                row.append("--timing")

            row.append("--no-auto-detection")

            row = " ".join(row) + "\n"
            print(row, end="")
            fw.write(row)
