from typing import Dict

from .BaseReport import BaseReport


class TextReport(BaseReport):
    template: Dict[str, str] = {}

    def _get_formatted_stats(self) -> str:
        formatted = []
        for stat in self.stats:
            formatted.append(self.template["stat_entry"].format(**stat))

        return "\n".join(formatted)

    def _get_formatted_table_entries(self):
        formatted = []
        for entry in self.entries:
            formatted.append(self.template["table_entry"].format(**entry))

        return "\n".join(formatted)


TextReport.template["main"] = """
Report: {title}
Statistics:
{stats}
"""
TextReport.template["stat_entry"] = "{name:<40}{value}"
TextReport.template["table_entry"] = """<tr>
    <td>{rank}</td>
    <td>{count}</td>
    <td>{codon}</td>
    <td class="smiles">{smiles}</td>
    <td><img src="data:image/png;base64,{imgb64}" /></td>
</tr>
"""