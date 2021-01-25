import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from typing import Optional, Tuple, List

from DECLGen.evaluation.evaluator import Evaluator


class ReportPlotter:
    evaluator: Evaluator
    minimum_axis: str
    minimum_uniques: List[float]

    def __init__(
        self,
        evaluator: Evaluator,
    ):
        self.evaluator = evaluator

        # Determine smallest axis
        minimum_axis = None
        minimum_uniques = -1
        minimum_axis_i = None
        for i, rank_column in enumerate(evaluator.codon_rank_columns):
            if minimum_axis is None:
                minimum_axis = rank_column
                minimum_axis_i = i
                minimum_uniques = evaluator.count_table_purified[rank_column].unique()
            else:
                uniques = evaluator.count_table_purified[rank_column].unique()
                if len(uniques) < len(minimum_uniques):
                    minimum_uniques = uniques
                    minimum_axis = rank_column
                    minimum_axis_i = i

        self.minimum_axis = minimum_axis
        self.minimum_uniques = minimum_uniques
        self.minimum_axis_i = minimum_axis_i

    def cross_correlation(self) -> Tuple[str, str]:
        if self.evaluator.replicates <= 1:
            return "Replicate cross-correlations", "<p>No cross-correlation possible, as number of replicates is 1"

        figure = px.scatter_matrix(
            self.evaluator.count_table_purified,
            dimensions=self.evaluator.normalized_column_names,
            title="Scatter matrix of mean counts of all replicates",
            labels={x: f"Replicate {i+1}" for i, x in enumerate(self.evaluator.normalized_column_names)},
            height=1000,
        )

        figure.update_traces(diagonal_visible=False)

        return "Replicate cross-correlations", figure.to_html(full_html=False)

    def count_scatters(self, top_hits: Optional[int] = 200):
        assert(len(self.evaluator.codon_rank_columns) == 3)

        all_data = self.evaluator.count_table_purified

        if len(self.minimum_uniques) <= 3:
            figure = make_subplots(
                rows=1, cols=len(self.minimum_uniques),
                specs=[[{"type": "scatter3d"} for x in range(len(self.minimum_uniques))]],
                column_width=[1/len(self.minimum_uniques) for x in range(len(self.minimum_uniques))],
            )
        else:
            figure = make_subplots(
                rows=len(self.minimum_uniques)//3+1, cols=3,
                specs=[[{"type": "scatter3d"} for x in range(3)] for y in range(len(self.minimum_uniques)//3)]
            )

        codons = [x for x in self.evaluator.codon_rank_columns if x != self.minimum_axis]
        codon_x, codon_y = codons[:2]

        for i, unique in enumerate(self.minimum_uniques):
            data = all_data[all_data[self.minimum_axis] == unique].nlargest(top_hits, "mean")

            row = i//3 + 1
            col = i % 3 + 1

            figure.add_trace(
                go.Scatter3d(
                    x=data[codon_x],
                    y=data[codon_y],
                    z=data["mean"],
                    mode="markers",
                    marker={
                        "color": data[codon_x],
                        "size": 2,
                        "opacity": 1.0,
                        "line": {
                            "color": "black",
                            "width": 0.01,
                        },
                    },
                    hovertext=data["Codon-Combination"],
                ),
                row=row, col=col
            )

            figure.update_xaxes({"title": {"text": codon_x}}, row=row, col=col)
            figure.update_yaxes({"title": {"text": codon_y}}, row=row, col=col)

        return "Top counts per sub codon", figure.to_html(full_html=False)

    def hits_preview(self, top_hits: Optional[int] = 250):
        if self.evaluator.diversity_elements_count == 2:
            pass
        elif self.evaluator.diversity_elements_count == 3:
            return f"3D Scatter (top {top_hits} hits)", self.three_element_3d_scatter(top_hits=top_hits)
        else:
            pass

    def three_element_3d_scatter(
        self,
        codon_x: int = 1,
        codon_y: int = 2,
        codon_z: int = 3,
        top_hits: Optional[int] = None,
    ):
        codon_x = self.evaluator.codon_rank_columns[codon_x - 1]
        codon_y = self.evaluator.codon_rank_columns[codon_y - 1]
        codon_z = self.evaluator.codon_rank_columns[codon_z - 1]
        data = self.evaluator.count_table_purified

        if top_hits:
            data = data.nlargest(top_hits, "mean")

        figure = go.Figure(
            data=go.Scatter3d(
                x=data[codon_x],
                y=data[codon_y],
                z=data[codon_z],
                mode="markers",
                marker={
                    "size": data["mean"]/data["mean"].max()*10,
                    "color": data[codon_y],
                    "opacity": 1.0,
                    "line": {
                        "color": "black",
                        "width": 0.01,
                    },
                },
                hovertext=data["Codon-Combination"],
            ),
            layout={
                "title": f"3 Codon plot, with top {len(data)} hits shown."
            }
        )

        figure_html = figure.to_html(full_html=False)

        return figure_html
