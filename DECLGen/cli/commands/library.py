from DECLGen import Runtime


def lib_info():
    """ Shows information about the library """
    r = Runtime()

    description = r.storage.library.describe()
    for key in description:
        print("{key:<15} {value}".format(key=key + ":", value=description[key]))
