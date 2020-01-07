import re

rSearch = re.compile(r"\[R([1-9][0-9]?)\]")
rReplace = r"[R\1]"
AuSearch = re.compile(r"\[([1-9][0-9]?)Au\]")
AuReplace = r"[\1Au]"


class RGroup():
    R: str
    Au: str

    def __init__(self, number: int):
        self.R = "[R{}]".format(number)
        self.Au = "[{}Au]".format(number)

    @classmethod
    def fromR(cls, R):
        return cls(int(R[1:]))

    @classmethod
    def fromAu(cls, Au):
        return cls(int(Au.split("Au")[0]))


def mask(string: str) -> str:
    result = rSearch.sub(AuReplace, string)

    return result

def unmask(string: str) -> str:
    result = AuSearch.sub(rReplace, string)

    return result

