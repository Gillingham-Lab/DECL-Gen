import os
import sys
from DECLGen import DECLGen_meta


def about():
    """ Shows information about the module """
    print("{name}: {version}".format(
        name=DECLGen_meta["name"],
        version=DECLGen_meta["version"]
    ))

    print("{description}".format(
        description=DECLGen_meta["description"]
    ))

    print("Use {} --help for help about how to use this application.".format(sys.argv[0]))

    print("\nExecutables: ")
    print("  Python: {}".format(sys.executable))
    print("  Path to this app: {}".format(os.path.join(os.path.dirname( sys.argv[0]), sys.argv[0])))
