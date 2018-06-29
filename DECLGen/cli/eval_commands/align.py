from DECLGen import Runtime
from DECLGen.exceptions import LibraryNoDNATemplateException

def align(
    fastq = "FastQ Read file"
):
    r = Runtime()

    if r.storage.library.get_dna_template() == None:
        a = LibraryNoDNATemplateException("You have not defined a DNA template.")
        r.error_exit(a)

