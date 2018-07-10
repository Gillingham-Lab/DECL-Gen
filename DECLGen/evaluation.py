
class AlignmentResult():
    """ Represents the result of codon alignment. """
    _result = None
    _key_annotations = {
        "reads_processed": ("Processed Reads", "both"),
        "reads_useful": ("Useful Reads", "both"),
        "valid_pairs": ("Valid Pairs", "paired"),
        "invalid_pairs": ("Invalid Pairs", "paired"),
        "low_quality_skips": ("Low Quality skips", "single"),
        "both_low_quality_skips": ("Low Quality skips (both)", "paired"),
        "r1_low_quality_skips": ("Low Quality skips (r1)", "paired"),
        "r2_low_quality_skips": ("Low Quality skips (r2)", "paired"),
    }
    paired = False
    _codons = None

    def __init__(self, paired=True):
        self._result = {key: 0 for key in self._key_annotations}
        self._codons = {}
        self._paired = paired

    def __getitem__(self, item):
        if item not in self._result:
            raise ValueError(
                "Result only supports a set amount of keys: {}, but {} given".format(
                    ", ".join(self._result.keys()),
                    item
                )
            )

        return self._result[item]

    def __setitem__(self, item, value):
        if item not in self._result:
            raise ValueError(
                "Result only supports a set amount of keys: {}, but {} given".format(
                    ", ".join(self._result.keys()),
                    item
                )
            )

        self._result[item] = value

    def __add__(self, othr):
        r = Result()
        # Add result meta
        for key in self._result:
            r[key] = self[key] + othr[key]

        # Merge codon lists
        self_codons = self.get_codons()
        othr_codons = othr.get_codons()

        for codon in self_codons:
            r.increase_codon(codon, self_codons[codon])
        for codon in othr_codons:
            r.increase_codon(codon, othr_codons[codon])

        return r

    def __str__(self):
        ret = []

        if self._paired:
            for key in self._result:
                annotation = self._key_annotations[key]
                value = self._result[key]

                if self._paired and annotation[1] in ["both", "paired"]:
                    ret.append("{0:<30} {1}".format(annotation[0], value))
                elif self._paired and annotation[1] in ["both", "single"]:
                    ret.append("{0:<30} {1}".format(annotation[0], value))

        return "\n".join(ret)

    def init_codon(self, codon) -> None:
        if codon not in self._codons:
            self._codons[codon] = 0

    def increase_codon(self, codon, increase=1):
        self.init_codon(codon)
        self._codons[codon] += increase

    def get_codons(self):
        return self._codons