from colorama import init, Fore
import gzip
import os


def react(
    molecules: "A list of smiles strings that should get modified.",
    product: "Filename to store the product of the reaction",
    reaction: "A predefined reaction" = None,
    rGroup: "Tne anchor (R1, R2, ..., R10, ..., Rn) to use with a predefined reaction" = None,
    smarts: "A smarts query to use instead" = None,
    protect: "Smart query to use for protection (for example, to protect amids or carbamate from modification)" = None,
):
    init()



    print("{Fore.GREEN}Run complete.{Fore.RESET}")

