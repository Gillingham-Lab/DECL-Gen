from colorama import init, Fore
import gzip
import os
from DECLGen import MoleculeContainer
from DECLGen.exceptions import MoleculeInvalidSmilesException


def react(
    inFile: "A list of smiles strings that should get modified.",
    outFile: "Filename to store the product of the reaction",
    reaction: "A predefined reaction" = None,
    rGroup: "Tne anchor (R1, R2, ..., R10, ..., Rn) to use with a predefined reaction" = None,
    smarts: "A smarts query to use instead" = None,
    protect: "Smart query to use for protection (for example, to protect amids or carbamate from modification)" = None,
):
    init()

    if not os.path.exists(inFile):
        print(f"{Fore.RED}Input file not found ({inFile}).{Fore.RESET}")
        return

    if os.path.exists(outFile):
        answer = input(f"Output file exists ({outFile}). Overwrite? (Y/n)").strip()

        if answer != "Y":
            print(f"{Fore.RED}Reaction cancelled.{Fore.RESET}")
            return

    try:
        molecules = MoleculeContainer.fromSmilesList(inFile)
    except MoleculeInvalidSmilesException as e:
        print(f"{Fore.RED}Error: {e.smiles} is not a valid smiles.{Fore.RESET}")
        exit(e.exitcode)

    print(f"{molecules.count()} molecules found in List.")






    print(f"{Fore.GREEN}Reaction complete.{Fore.RESET}")

