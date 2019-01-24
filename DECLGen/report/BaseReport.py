from typing import Dict, List, TypeVar, Any

Real = TypeVar("real", int, float)
class BaseReport:
    template: Dict[str, str] = {}
    filename: str
    title: str
    stats: List[Dict[str, str]]
    entries: List[Dict[str, Any]]

    def __init__(self, filename: str, title: str = None):
        self.filename = filename
        self.stats= []
        self.entries = []
        self.title = ""

        if title is not None:
            self.title = title

    @classmethod
    def copyFromOther(cls, other: "BaseReport"):
        instance = cls(other.filename, other.title)
        instance.stats = other.stats
        instance.entries = other.entries

        return instance

    def set_title(self, title: str):
        self.title = title

    def add_stats(self, name: str, value: str):
        self.stats.append({"name": name, "value": value})

    def _get_formatted_stats(self):
        raise NotImplementedError()

    def add_entry(self, rank: int, count: Real, codons: str, smiles: str, image: bytes):
        self.entries.append({
            "rank": rank,
            "count": count,
            "codon": codons,
            "smiles": smiles,
            "imgb64": image.decode("utf8")})

    def _get_formatted_table_entries(self):
        raise NotImplementedError()

    def _get_formatted_other(self):
        raise NotImplementedError()

    def save(self) -> None:
        content = str(self)

        with open(self.filename, "w") as file:
            file.write(content)

    def __str__(self) -> str:
        content = self.template["main"].format(
            title=self.title,
            stats=self._get_formatted_stats(),
            table_entries=self._get_formatted_table_entries(),
            other=self._get_formatted_other(),
        )

        return content

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.save()