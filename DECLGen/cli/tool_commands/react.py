import argh
from colorama import init, Fore
import gzip
import os

from DECLGen import MoleculeContainer
from DECLGen.exceptions import MoleculeInvalidSmilesException, ReactionNoProductException
from DECLGen.reactions import Reaction, predefinedReactions


@argh.arg("--reaction", choices=list(predefinedReactions.keys()))
@argh.arg("--smarts", nargs="+")
@argh.arg("--protect", nargs="+")
def react(
    inFile: "A list of smiles strings that should get modified.",
    outFile: "Filename to store the product of the reaction",
    reaction: "A predefined reaction" = None,
    rGroup: "Tne anchor (R1, R2, ..., R10, ..., Rn) to use with a predefined reaction" = None,
    smarts: "A smarts query to use to install the R-Group. You must specify the R-group and always use the same." = None,
    protect: "Smart query to use for protection (for example, to protect amids or carbamate from modification). All are applied at the same time." = None,
    force_yes: "Set to true if answers should be assumed to be yes" = False,
):
    init()

    # Error if inFile does not exist.
    if not os.path.exists(inFile):
        print(f"{Fore.RED}Input file not found ({inFile}).{Fore.RESET}")
        return

    # Ask for overwrite permission if outFile exists.
    if os.path.exists(outFile):
        if force_yes is False:
            answer = input(f"Output file exists ({outFile}). Overwrite? (Y/n)").strip()

            if answer != "Y":
                print(f"{Fore.RED}Reaction cancelled.{Fore.RESET}")
                return
        else:
            print(f"Output file ({outFile}) exists, but will be overwritten.")

    # Check if either a reaction is given, or a smarts string.
    if reaction is None and smarts is None:
        print(f"{Fore.RED}You must give either a predefined reaction or a smarts string.{Fore.RESET}")
        return

    if reaction is not None and smarts is not None:
        print(f"{Fore.RED}You specified a predefined reaction. The smarts and protect values are ignored.{Fore.RESET}")

    # Check if rgroup is given
    if rGroup is None and reaction is not None:
        print(f"{Fore.RED}If you use a predefined reaction you MUST specify the rGroup parameter.")
        return

    # Check if rgroup is proper
    if rGroup is not None:
        try:
            r, val = rGroup[0], int(rGroup[1:])
        except ValueError:
            print(f"{Fore.RED}A R-Group must consist of the letter R followed by 1 or 2 digits.{Fore.RESET}")
            return

        if r != "R":
            print(f"{Fore.RED}A R-Group MUST start with an R.{Fore.RESET}")
            return

        if not 0 < val < 90:
            print(f"{Fore.RED}The R-Group number MUST be > 0 and < 90.{Fore.RESET}")
            return

    # Create list of molecules
    try:
        molecules = MoleculeContainer.fromSmilesList(inFile)
    except MoleculeInvalidSmilesException as e:
        print(f"{Fore.RED}Error: {e.smiles} is not a valid smiles.{Fore.RESET}")
        exit(e.exitcode)

    print(f"{molecules.count()} molecules found in List.")

    # If a predefined reaction, get if from the database
    if reaction is not None:
        reactionClass = Reaction(**predefinedReactions[reaction], rGroup=rGroup)
    elif smarts is not None:
        reactionClass = Reaction(smarts, protect)
    else:
        # For completion, should not arrive here.
        return

    # Run the reactions
    for m in molecules:
        try:
            reactionClass.react(m)
        except ReactionNoProductException as e:
            print(f"No product detected with {e.smiles}.")

    # Save the results
    molecules.export(outFile)

    # Print OK
    print(f"{Fore.GREEN}Reaction complete and saved.{Fore.RESET}")

