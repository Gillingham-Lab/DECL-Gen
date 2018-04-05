from DECLGen import Runtime


def lib_info():
    """ Shows information about the library """
    r = Runtime()

    description = r.storage.library.describe()
    for key in description:
        print("{t.bold}{key:<15}{t.normal} {value}".format(t=r.t, key=key + ":", value=description[key]))
