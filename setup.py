from setuptools import setup, find_packages

from .DECLGen import DECLGen_meta

build_exe_options = {
    "includes": [
        "sys",
        "os"
    ],
    "excludes": [

    ],
}

setup(**DECLGen_meta)