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

    print("Use declGen --help for help about how to use this application.")
