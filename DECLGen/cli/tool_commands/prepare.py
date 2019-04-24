from colorama import init, Fore
import gzip
import os


def prepare(
    source: "Sequencing main folder containing each sample in a separate directory (as .tar.gz files)",
    target: "Target folder to put unzipped files into" = None,
    glue: "If true, unzipped sequencing files are glues together" = False,
):
    init()

    # Check paths
    if not os.path.exists(source):
        print("{Fore.RED}The given source path was not found ({}).{Fore.RESET}".format(source, Fore=Fore))
        return

    if target and not os.path.exists(target):
        print("{Fore.RED}The given target path was not found ({}).{Fore.RESET}".format(source, Fore=Fore))
        return

    elif not target:
        print("No target given, files will be extracted in {}.\n".format(source))
        target = source

    directories = None
    directory_count = 0

    for dirpath, dirnames, filenames in os.walk(source):
        # Root folder: Save dirnames
        if directories is None:
            directories = dirnames
            continue

        # Increase directory count by 1
        directory_count += 1

        # Sort filenames alphabetically
        filenames = sorted(filenames)
        filenames = [x for x in filenames if x.endswith("fq.gz") or x.endswith("fastq.gz")]

        print("{Fore.CYAN} - {path} ({i}/{len}){Fore.RESET}".format(path=dirpath, i=directory_count, len=len(directories), Fore=Fore))

        if len(filenames) == 0:
            print("No sequencing files found.\n")
            continue


        samplename = os.path.basename(dirpath)
        todo = {}

        for filename in filenames:
            # Lets figure out how we should put the sample together
            if not filename.startswith(samplename):
                print("{Fore.RED}File {} does not start with the name of the directory ({}). Auto-preparation failed.{Fore.RESET}".format(filename, samplename, Fore=Fore))
                break

            runname, cell_id, lane, mate = filename[len(samplename)+1:].split("_")
            mate = mate.split(".")[0]

            # If glue is used, we only distinguish the samples by using the samplename and mate identifier because the rest is assumed to be not interesting.
            # If not, we use the full filename for it (without extension)
            if glue:
                key = "_".join([samplename, mate])
            else:
                key = "_".join([samplename, runname, cell_id, lane, mate])

            if key not in todo:
                todo[key] = []

            todo[key].append(os.path.join(dirpath, filename))

        # We work through the to-do list. For each key, we create a new file based on that key.
        # The corresponding item is always a list. All files (1 or more) are read and saved into
        # the aforementioned file.
        for key in todo:
            target_filename = os.path.join(target, key + ".fq")

            with open(target_filename, "w") as fw:
                print("   - Write {}".format(os.path.basename(target_filename)))

                for subfile in todo[key]:
                    with gzip.open(subfile, "rt", encoding="ASCII") as fr:
                        print("     - Read {}".format(os.path.basename(subfile)))
                        for line in fr:
                            fw.write(line)

        print()

    print("{Fore.GREEN}Run complete.{Fore.RESET}")

