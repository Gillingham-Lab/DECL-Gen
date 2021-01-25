import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from typing import Tuple
import numpy as np


class ComparePlotter:
    data: pd.DataFrame

    def __init__(self, data):
        self.data = data

    def enrichment_plot(self) -> Tuple[str, str]:
        figure = go.Figure(
            data=go.Scattergl(
                x=self.data["mean_signal"],
                y=self.data["mean_signal"]/self.data["mean_noise"],
                mode="markers",
                marker={
                    "color": self.data["score"],
                    "colorscale": "Viridis",
                    "size": 4,
                    "opacity": 1.0,
                    "line": {
                        "color": "black",
                        "width": 0.01,
                    },
                },
                hovertext=self.data["Codon-Combination"],
            ),
            layout={
                "title": f"Enrichment plot."
            }
        )

        figure.update_xaxes({"title": {"text": "Mean counts (Signal)"}})
        figure.update_yaxes({"title": {"text": "Mean counts (Signal/Noise)"}})

        figure_html = figure.to_html(full_html=False)

        return "Enrichment plot.", figure_html

    def vulcano_plot_p(self) -> Tuple[str, str]:
        figure = go.Figure(
            data=go.Scattergl(
                x=self.data["log-Fold"],
                y=-np.log10(self.data["pValue"]),
                mode="markers",
                marker={
                    "color": self.data["score"],
                    "colorscale": "Viridis",
                    "size": 4,
                    "opacity": 1.0,
                    "line": {
                        "color": "black",
                        "width": 0.1,
                    },
                },
                hovertext=self.data["Codon-Combination"],
            ),
            layout={
                "title": f"Vulcano-Plot."
            }
        )

        figure.update_xaxes({"title": {"text": "log-fold enrichment"}})
        figure.update_yaxes({"title": {"text": "p-Value"}})

        figure_html = figure.to_html(full_html=False)

        return "Vulcano plot (p-Value)", figure_html

    def vulcano_plot_q(self) -> Tuple[str, str]:
        figure = go.Figure(
            data=go.Scattergl(
                x=self.data["log-Fold"],
                y=-np.log10(self.data["qValue"]),
                mode="markers",
                marker={
                    "color": self.data["score"],
                    "colorscale": "Viridis",
                    "size": 4,
                    "opacity": 1.0,
                    "line": {
                        "color": "black",
                        "width": 0.1,
                    },
                },
                hovertext=self.data["Codon-Combination"],
            ),
            layout={
                "title": f"Vulcano-Plot with Benjamini-Hochberg correction."
            }
        )

        figure.update_xaxes({"title": {"text": "log-fold enrichment"}})
        figure.update_yaxes({"title": {"text": "q-Value"}})

        figure_html = figure.to_html(full_html=False)

        return "Vulcano plot (q-Value)", figure_html