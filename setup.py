#!/usr/bin/env python
from setuptools import setup, find_packages

build_exe_options = {
    "includes": [
        "sys",
        "os"
    ],
    "excludes": [

    ],
}

setup(
    name="DECL-Gen",
    version="v1.0-alpha",
    description="A tool to generate DNA encoded libraries from fragments",
    url="https://github.com/Gillingham-Lab/DECL-Gen",

    author="Basilius Sauter",
    author_email="basilius.sauter@unibas.ch",

    python_requires='>3.6.0',
    packages=find_packages(),
    install_requires=[
        "argh",
        "colorama",
        "numpy",
        "scipy",
        "pandas",
        "seaborn",
        "BioPython",
        "PIL",
        "rdkit",
        "statsmodels",
    ],

    entry_points={
        "console_scripts": [
            'declGen=DECLGen.cli.generate:main',
            'declEval=DECLGen.cli.evaluate:main',
        ]
    }
)