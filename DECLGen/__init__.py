

DECLGen_meta = dict(
    name="DECL-Gen",
    version="v1.0-alpha",
    description="A tool to generate DNA encoded libraries from fragments",
    url="https://github.com/Gillingham-Lab/DECL-Gen",

    author="Basilius Sauter",
    author_email="basilius.sauter@unibas.ch",

    packages=find_packages(),
    install_requires = [
        "argh",
        "rdkit"
    ],

    entry_points = {
        "console_scripts": [
            'declGen=DECLGen.cli.main:main',
        ]
    }
)