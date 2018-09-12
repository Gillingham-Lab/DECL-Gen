from .qc_simple import qc as qc_simple
from .qc_align import qc as qc_align
from .qc_score import qc as qc_score
from .qc_advanced import qc as qc_advanced

_all = [qc_simple, qc_score, qc_align, qc_advanced]

class Type():
    Simple = qc_simple
    Score = qc_score
    Align = qc_align
    Advanced = qc_advanced

    default = qc_simple.__name__
    all = [x.__name__ for x in _all]
    _all = {x.__name__: x for x in _all}

    @staticmethod
    def get(v):
        return Type._all[v]