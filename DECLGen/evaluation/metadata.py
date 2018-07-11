from enum import Enum
from typing import List, Tuple, Optional
from Bio.Seq import Seq

class ReadfileType(Enum):
    forward = 0
    reverse = 1

class ReadfileMetadata():
    filename = None
    template = None
    coordinates = None
    type = None

    def __init__(self, filename: str, template: Seq, coordinates: List[Tuple[int, int]], type=ReadfileType):
        self.filename = filename
        self.template = template
        self.coordinates = coordinates
        self.type = type

    def is_reverse(self):
        return self.type == ReadfileType.reverse

    def is_forward(self):
        return self.type == ReadfileType.forward

    def is_none(self):
        if self.filename is None:
            return True

        return False


class ReadfileWorkerMetadata():
    r1 = None
    r2 = None
    compare_n = None
    blocksize = None
    checktype = None

    def __init__(self,
                 r1: ReadfileMetadata,
                 r2: Optional[ReadfileMetadata],
                 compare_n: int,
                 blocksize: int,
                 checktype,
                 ):
        """

        :param r1:
        :param r2:
        :param compare_n:
        :param blocksize:
        :param checktype: DECLGen.evaluation.qc.Type
        """
        self.r1 = r1
        self.r2 = r2
        self.compare_n = compare_n
        self.blocksize = blocksize
        self.checktype = checktype

    def is_single(self):
        return True if self.r2 is None else False

    def is_paired(self):
        return False if self.r2 is None else True