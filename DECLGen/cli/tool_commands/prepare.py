from colorama import init, Fore
import gzip
import os
import shutil


def prepare(
    source: "Sequencing main folder containing each sample in a separate directory (as .tar.gz files)",
    target: "Target folder to put unzipped files into" = None,
    glue: "If true, unzipped sequencing files are glues together" = False,
    no_decompress: "Set to true if you only want to copy the file." = False,
    simplify_name: "Simplifies the filename by using the directory name." = False,
):
    """
    Prepares a sequencing run by copying './*/sample.fq.gz' files to './sample.fq'.
    """
    init()

    # Check paths
    if not os.path.exists(source):
        print(f"{Fore.RED}The given source path was not found ({source}).{Fore.RESET}")
        return

    if target and not os.path.exists(target):
        print(f"{Fore.RED}The given target path was not found ({source}).{Fore.RESET}")
        return

    elif not target:
        print(f"No target given, files will be extracted in {source}.\n")
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

        print(f"{Fore.CYAN} - {dirpath} ({directory_count}/{len(directories)}){Fore.RESET}")

        if len(filenames) == 0:
            print("No sequencing files found.\n")
            continue


        samplename = os.path.basename(dirpath)
        todo = {}

        for filename in filenames:
            # Lets figure out how we should put the sample together
            if not filename.startswith(samplename):
                print(f"{Fore.RED}File {filename} does not start with the name of the directory ({samplename}). "
                      f"Auto-preparation failed.{Fore.RESET}")
                break

            runname, cell_id, lane, mate = filename[len(samplename)+1:].split("_")
            mate = mate.split(".")[0]

            # If glue is used, or we call for simplfied filenames, we only distinguish the samples by using the
            # samplename (same as subdirectory name) and mate identifier, as the rest is assumed to be not interesting.
            # If not, we use the full filename for it.
            # Key does not contain a file extension.
            if glue or simplify_name:
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

            # If we don't want glue, and decompression is not wanted, we can just copy the file
            if no_decompress is True and glue is False:
                shutil.copyfile(todo[key][0], target_filename + ".gz")
            # If not, we can gzip-read the file and save it into a new (gzipped-) file.
            else:
                # Determine the way to open the target file
                # If decompression is wanted (default, no_decompress==true), we use the normal open command
                # If not, we use gzip.open instead to write a compressed file.
                if no_decompress is False:
                    opener = (open, [target_filename, "w"], {})
                else:
                    target_filename += ".gz"
                    opener = (gzip.open, [target_filename, "wt"], {"encoding": "ASCII"})

                # Write into the new file
                with opener[0](*opener[1], **opener[2]) as fw:
                    print(f"   - Write {os.path.basename(target_filename)}")

                    for subfile in todo[key]:
                        with gzip.open(subfile, "rt", encoding="ASCII") as fr:
                            print(f"     - Read {os.path.basename(subfile)}")
                            for line in fr:
                                fw.write(line)

        print()

    print(f"{Fore.GREEN}Run complete.{Fore.RESET}")

