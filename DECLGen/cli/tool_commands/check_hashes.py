import hashlib
import os
from colorama import init, Fore, Style


def md5_file(file):
    """ Memory efficient calculation of md5 hashes of files. """
    md5 = hashlib.md5()

    with open(file, "rb") as fh:
        for chunk in iter(lambda: fh.read(4096), b""):
            md5.update(chunk)

    return md5.hexdigest()


def check_hashes(
    path: "Path to start searchning."
):
    """ Checks recursively the md5 hashes against a list in a given direction. """
    init()

    if not os.path.exists(path):
        print("{Fore.RED}The given path was not found ({}){Fore.RESET}".format(path, Fore=Fore))

    errors = 0
    files = 0

    # Walk through the directory
    for (dirpath, dirnames, filenames) in os.walk(path):
        hashes_should = {}
        hashes_is = {}

        # In each directory, walk through all files.
        for filename in filenames:
            # With the Illumina platform, there should be a md5.txt file contained in each sample directory.
            # This file contains a list of md5 hashes and filenames, separated by two spaces.
            # We load in this file and save the result in hashes_should (key; filename).
            if filename.lower() == "md5.txt":
                with open(os.path.join(dirpath, filename), "r") as md5File:
                    for row in md5File:
                        # Strip whitespaces from the row (gets rid of \n at the end)
                        row = row.strip()

                        # Split exactly once.
                        md5, filename = row.split("  ", 1)

                        hashes_should[filename] = md5
                continue

            # Calculate the md5 hash of the given file.
            md5 = md5_file(os.path.join(dirpath, filename))

            # Put it into hashes_is (key: filename)
            hashes_is[filename] = md5

        # If hashes_should is empty, there was no md5 file in the directory and we can skip it.
        if len(hashes_should) == 0:
            continue

        # Output results of the directory.
        print("Checking hashes in {}".format(dirpath))
        # Iterate over the keys in hashes-should
        for filename in hashes_should:
            files += 1

            # If the hashes are teh same, we print success; if not, we print the failure and increase errors by 1.
            if hashes_is[filename] == hashes_should[filename]:
                print("[{Fore.GREEN}X{Fore.RESET}] {}".format(filename, Fore=Fore))
            else:
                print("[{Fore.RED}X{Fore.RESET}] {}".format(filename, Fore=Fore))
                errors += 1
        print("\n")

    # Summarize the results.
    print("\n{} files have been checked".format(files))

    if errors == 0:
        print("{Fore.GREEN}All files are good.{Fore.RESET}".format(Fore=Fore))
    else:
        print("{Fore.RED}There have been errors in {} files.\nCheck which one failed. Try to re-unpack the tarball and check again." \
            + "If the error persists, consider to download the data again. If it still persists, contact the sequencing company. " \
            + "{Fore.RESET}".format(errors, Fore=Fore))
