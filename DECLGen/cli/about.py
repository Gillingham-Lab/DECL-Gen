import os
import sys
from colorama import init, Fore, Style
import DECLGen


def about():
    """ Shows information about the module """
    init()

    print("""
    {highlight}{name} {version} {description}
    Copyright (C) {year} {author}

    {lowlight}This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.{reset}
    """.format(**{
        "name": DECLGen.__name__,
        "version": DECLGen.__version__,
        "description": DECLGen.__description__,
        "year": DECLGen.__year__,
        "author": DECLGen.__author__,
        "highlight": Fore.LIGHTWHITE_EX,
        "lowlight": Fore.LIGHTBLACK_EX,
        "help": Fore.GREEN,
        "reset": Fore.RESET
    }))



    print(Fore.GREEN + "Use {} --help for help about how to use this application.".format(sys.argv[0]) + Fore.RESET)

    print(Fore.LIGHTBLACK_EX, end="")
    print("\nExecutables: ")
    print("  Python: {}".format(sys.executable))
    print("  Path to this app: {}".format(os.path.join(os.path.dirname( sys.argv[0]), sys.argv[0])))
    print(Fore.RESET, end="")
