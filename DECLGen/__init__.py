from .runtime import Runtime
from .storage import Storage
from . import template
from .exceptions import *
from .terminal import Terminal
from . import codon

DECLGen_meta = dict(
    name="DECL-Gen",
    version="v1.0-alpha",
    description="A tool to generate DNA encoded compound library from fragments and calculate their properties.",
    url="https://github.com/Gillingham-Lab/DECL-Gen",

    author="Basilius Sauter",
    author_email="basilius.sauter@unibas.ch",
    year="2018",

    install_requires=[
        "argh",
        "rdkit"
    ],

    entry_points={
        "console_scripts": [
            'declGen=DECLGen.cli.main:main',
        ]
    }
)