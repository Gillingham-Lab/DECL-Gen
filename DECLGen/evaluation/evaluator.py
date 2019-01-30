import base64
import io
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors
from matplotlib.markers import MarkerStyle
from matplotlib import path as mpath
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle
import numpy as np
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit import RDLogger
import seaborn as sns
from typing import Sequence, List, Tuple, Optional

from DECLGen.exceptions import EvaluationFileDoesNotExist, EvaluationInvalidFileFormat
from DECLGen.codon import decode

def calculate_enrichment_factor(k: int, N: int, n: int) -> float:
    """
    Normalizes the count by taking sequencing depth into account.

    This method normalizes according to F. Samain et al., J. Med. Chem. 2015, 58, 12, 5143-5149.
    :param k: copy count
    :param N: library diversity
    :param n: read count
    :return: The normalized enrichment factor e.
    """
    return k*N/n

def estimate_enrichment_confidence_boundaries(k: int, N: int, n: int) -> Tuple[float, float]:
    """
    Estimates the lower and upper confidence boundaries for the enrichment factor for 95% confidence interval.

    Bases on a Method described by L. Kuai et al, 10.1177/2472555218757718
    :param k:
    :param N:
    :param n:
    :return:
    """
    def Ř(x):
        return x**(1/2)

    lower = (Ř(k+1)-1)**2 * N/n
    upper = (Ř(k+1)+1)**2 * N/n

    return lower, upper


class Evaluator:
    count_table: pd.DataFrame
    count_table_purified: pd.DataFrame
    count_table_reduced: pd.DataFrame
    replicates: int
    count_columns: List[str]
    codon_columns: List[str]
    codon_rank_columns: List[str]
    column_numbers: List[int]
    diversity_elements_count: int

    unique_codons: int
    all_counts: int
    valid_codons: int
    valid_counts: int

    plot_format: str


    def __init__(self, *results, progressBar = None, plot_format = "png", properties: pd.DataFrame):
        self.library_size = len(properties)
        self.count_table = None
        self.progressBar = progressBar
        self.load_counts(results)
        self.plot_format = plot_format

        self.merge_properties(properties)

    @property
    def expected_missing(self):
        N = self.library_size
        p = []

        for c in self.count_columns:
            x = 1 - (1 - 1/N)**self.count_table_purified[c].sum()
            p.append(x)

        P = N
        for i in p:
            P *= i
        return self.library_size - P

    def _get_encoding(self):
        if self.plot_format == "svg":
            return "image/svg+xml;base64"
        return "image/png;base64"

    def merge_properties(self, properties: pd.DataFrame):
        self.library_size = len(properties)

        self.count_table_purified = properties.merge(self.count_table, on="Codon-Combination", how="left")
        self.count_table_purified.sort_values(by="Codon-Combination", inplace=True)
        self.count_table_purified.reset_index(drop=True)

        self.diversity_elements_count = len(self.count_table_purified["Codon-Combination"].iloc[0].split("-"))

        # We do not know the codon column names. Thats why we create our own, and rank them for plotting
        de_columns = []
        fill = False
        for col in properties.columns:
            if fill is False:
                if col == "Codon-Combination":
                    fill = True
                continue
            if fill is True:
                if col == "DNA":
                    fill = False
                else:
                    de_columns.append(col)

        self.codon_columns = de_columns
        self.codon_rank_columns = []
        for de_col in de_columns:
            new_col = "CodonRank_{}".format(de_col)
            self.count_table_purified[new_col] = self.count_table_purified[de_col].apply(decode) + 1
            self.codon_rank_columns.append(new_col)


        # Create a copy of the table where non-overlapping are removed.
        self.count_table_reduced = self.count_table_purified.dropna()
        self.count_table_reduced = self.count_table_reduced.sort_values(by=["MeanCounts", "StdDevCounts"], ascending=[False, True])
        self.count_table_reduced = self.count_table_reduced.reset_index(drop=True)

        # Create a copy of the table where non-overlapping are kept, but set to 0.
        self.count_table_purified = self.count_table_purified.fillna(0)
        self.count_table_purified = self.count_table_purified.sort_values(by=["MeanCounts", "StdDevCounts"], ascending=[False, True])
        self.count_table_purified = self.count_table_purified.reset_index(drop=True)

        # Statistics
        self.valid_codons = self.count_table_purified[self.count_table_purified["MeanCounts"] > 0]["Codon-Combination"].nunique()
        self.valid_counts = int(self.count_table_purified[self.count_columns].sum().sum())

        if self.progressBar is not None:
            self.progressBar.update(0.99)

    def save_df(self, filename):
        self.count_table_purified.to_csv(filename, sep="\t", index=False)

    def load_counts(self, results: Sequence):
        """
        Loads given result count files and merges then into one dataframe, storing the mean and standard deviation of the counts.

        The original count columns are kept intact, but they are enumerated.
        :param results:
        :return:
        """
        result = None
        filename = None
        column = None
        columns = None
        count_columns = []
        count_column_names = []
        column_numbers = []
        replicates = len(results)

        try:
            i = 0
            for filename in results:
                i+=1
                if self.progressBar is not None:
                    self.progressBar.update(i / (replicates + 1))

                df = pd.read_csv(filename, sep="\t")
                count_column_names.append(os.path.basename(filename))

                if column is None:
                    if "CountRel" in df.columns:
                        column = "CountRel"
                    elif "Count" in df.columns:
                        column = "Count"
                    else:
                        raise EvaluationInvalidFileFormat(
                            ("Cannot determine count column of result file {}. Make sure the file has been generated by " +
                            "declEval extract or declEval compare.").format(filename))

                    columns = [x for x in df.columns]
                else:
                    if [x for x in df.columns] != columns:
                        raise EvaluationInvalidFileFormat(
                            ("The columns of the result file {} does not match those loaded before. Make sure they " +
                            "were properly created.").format(filename))

                # If first result, just load in.
                if result is None:
                    result = df
                else:
                    result = result.merge(df, on="Codon-Combination", how="outer")

                new_column = column + "{}".format(i)
                column_number = format(i)

                result.rename(columns={column: new_column}, inplace=True)
                result["Enrichment{}".format(i)] = calculate_enrichment_factor(result[new_column], self.library_size, result[new_column].sum())
                lower, upper = estimate_enrichment_confidence_boundaries(result[new_column], self.library_size, result[new_column].sum())
                result["Enrichment_Lower{}".format(i)] = lower
                result["Enrichment_Upper{}".format(i)] = upper

                count_columns.append(new_column)
                column_numbers.append(column_number)

            # Fill empty fields with 0
            result.fillna(0, inplace=True)

            # Save mean and standard deviation
            result["MeanCounts"] = result[count_columns].mean(axis=1)
            result["StdDevCounts"] = result[count_columns].std(axis=1)

            result["MeanEnrichment"] = result[["Enrichment{}".format(i) for i in column_numbers]].mean(axis=1)
            result["StdEnrichment"] = result[["Enrichment{}".format(i) for i in column_numbers]].std(axis=1)

            #result["MeanEnrichment"] = calculate_enrichment_factor(result["MeanCounts"], self.library_size, result["MeanCounts"].sum() * replicates)
            #result["StdEnrichment"] = calculate_enrichment_factor(result["StdDevCounts"], self.library_size, result["MeanCounts"].sum() * replicates)
            lower, upper = estimate_enrichment_confidence_boundaries(result["MeanEnrichment"], 1, 1)

            result["MeanEnrichment_Lower"] = lower
            result["MeanEnrichment_Upper"] = upper
            result["VarEnrichment"] = result["StdEnrichment"]**2
            result.sort_values(by=["MeanEnrichment", "StdDevCounts"], ascending=False, inplace=True)
            #result.drop(columns=["Enrichment{}".format(i) for i in column_numbers], inplace=True)
            result.drop(columns=["Enrichment_Lower{}".format(i) for i in column_numbers], inplace=True)
            result.drop(columns=["Enrichment_Upper{}".format(i) for i in column_numbers], inplace=True)

        except FileNotFoundError:
            raise EvaluationFileDoesNotExist("Result file <{}> does not exist.".format(filename))

        # Save
        self.count_table = result
        self.replicates = replicates
        self.count_columns = count_columns
        self.count_column_names = count_column_names
        self.column_numbers = column_numbers

        self.count_table.rename(dict(zip(self.count_columns, self.count_column_names)), axis="columns", inplace=True)
        self.count_columns, self.count_column_names = self.count_column_names, self.count_columns

        # Statistics
        self.unique_codons = self.count_table["Codon-Combination"].nunique()
        self.all_counts = int(self.count_table[self.count_columns].sum().sum())

    def replicates_scatter(self, scale = None):
        if self.replicates <= 1:
            return

        with io.BytesIO() as image_stream:
            ax = sns.pairplot(
                self.count_table_purified[["Enrichment{}".format(i) for i in self.column_numbers]], diag_kind="kde", kind="reg", markers="+",
                plot_kws=dict(),
                diag_kws=dict(shade=True),
            )

            plt.savefig(image_stream, format=self.plot_format)
            img_base64 = base64.b64encode(image_stream.getvalue())
            plt.close()

        return img_base64

    def variance_vs_mean(self):
        if self.replicates <= 1:
            return

        with io.BytesIO() as image_stream:
            ax = sns.regplot(x="MeanEnrichment", y="VarEnrichment", data=self.count_table_purified)
            ax = sns.regplot(x="MeanEnrichment", y="MeanEnrichment", data=self.count_table_purified, scatter=False, ax=ax, label=False)

            plt.savefig(image_stream, format=self.plot_format)
            img_base64 = base64.b64encode(image_stream.getvalue())
            plt.close()

        return img_base64

    def count_histogram(self):
        with io.BytesIO() as image_stream:
            ax = sns.kdeplot(self.count_table_reduced["MeanEnrichment"], shade=True, legend="mean counts")

            plt.savefig(image_stream, format=self.plot_format)
            img_base64 = base64.b64encode(image_stream.getvalue())
            plt.close()

        return img_base64

    def two_element_3d_scatter(self, top_hits = None, codon_x = 1, codon_y = 2, anchored=False, project_on_z_plane=False):
        codon_x = self.codon_rank_columns[codon_x - 1]
        codon_y = self.codon_rank_columns[codon_y - 1]

        df = self.count_table_reduced.sort_values(by=[codon_x, codon_y], ascending=[False, False])

        if top_hits is not None:
            df_reduced = df.nlargest(top_hits, "MeanEnrichment").sort_values(by=[codon_x, codon_y], ascending=[False, False])
        else:
            df_reduced = df

        with io.BytesIO() as image_stream:
            fig, ax = Scatter3D(
                df_reduced, codon_x, codon_y, "MeanEnrichment", c="MeanEnrichment", s="MeanEnrichment", complete_data=df,
                project_on_z_plane=project_on_z_plane,
                anchored=anchored,
            )

            fig.savefig(image_stream, format=self.plot_format)
            img_base64 = base64.b64encode(image_stream.getvalue())
            plt.close()

        return img_base64

    def three_element_3d_scatter(self, top_hits = None, codon_x = 1, codon_y = 2, codon_z = 3, anchored=False, project_on_z_plane=False):
        codon_x = self.codon_rank_columns[codon_x - 1]
        codon_y = self.codon_rank_columns[codon_y - 1]
        codon_z = self.codon_rank_columns[codon_z - 1]

        df = self.count_table_reduced.sort_values(by=[codon_x, codon_y, codon_z], ascending=[False, False, False])

        if top_hits is not None:
            df_reduced = df.nlargest(top_hits, "MeanEnrichment").sort_values(by=[codon_x, codon_y, codon_z],
                                                                             ascending=[False, False, False])
        else:
            df_reduced = df

        with io.BytesIO() as image_stream:
            fig, ax = Scatter3D(
                df_reduced, codon_x, codon_y, codon_z, c="MeanEnrichment", s="MeanEnrichment",
                complete_data=df,
                project_on_z_plane=project_on_z_plane,
                anchored=anchored,
            )

            fig.savefig(image_stream, format=self.plot_format)
            img_base64 = base64.b64encode(image_stream.getvalue())
            plt.close()

        return img_base64


    def iter_hits(self, top=None):
        lg = RDLogger.logger().setLevel(RDLogger.ERROR)

        i = 0
        for index, row in self.count_table_reduced.iterrows():
            i += 1

            if i > top:
                break

            yield Hit(i, row)

        lg = RDLogger.logger().setLevel(RDLogger.INFO)

class Hit:
    rank: int
    row: pd.Series

    def __init__(self, rank, row):
        self.rank = rank
        self.row = row

    @property
    def counts(self):
        return self.row["MeanEnrichment"]

    @property
    def codons(self):
        return self.row["Codon-Combination"]

    @property
    def smiles(self):
        return self.row["Canonical smiles"]

    def draw(self) -> str:
        m = Chem.MolFromSmiles(self.row["Canonical smiles"])
        AllChem.Compute2DCoords(m)

        with io.BytesIO() as image_stream:
            img = Draw.MolToImage(m, size=(300, 300), kekulize=True, wedgeBonds=True)
            img.save(image_stream, format="png")

            img_base64 = base64.b64encode(image_stream.getvalue())

        return img_base64



def Scatter3D(
        data: pd.DataFrame,
        x: str, y: str, z: str,
        ax=None,
        anchored: bool=False,
        c: str=None,
        s: str=None,
        project_on_x_plane: bool=False,
        project_on_y_plane: bool=False,
        project_on_z_plane: bool=False,
        complete_data: pd.DataFrame=None,
        azim = -135,
        elev = 20,
):
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d', azim=azim, elev=elev)

    lim_data = complete_data if complete_data is not None else data
    x_lim = lim_data[x].min(), lim_data[x].max()
    y_lim = lim_data[y].min(), lim_data[y].max()
    z_lim = lim_data[z].min(), lim_data[z].max()

    color_kwargs = {}
    if c is not None:
        color_kwargs = {"c": data[c], "cmap": "plasma", "norm": colors.PowerNorm(gamma=1),}

    size_kwargs = {}
    if s is not None:
        size_kwargs = {"s": (data[s] / data[s].max()) ** 2 * 100}

    if project_on_z_plane is True:
        z_projection = ax.scatter(
            data[x], data[y], z_lim[0], c="lightgray", depthshade=False,
            **size_kwargs,
            zorder=-100,
            marker="."
        )

    if project_on_y_plane is True:
        y_projection = ax.scatter(
            data[x], y_lim[1], data[z], c="lightgray", depthshade=False,
            **size_kwargs,
            zorder=-100,
        )

    if project_on_x_plane is True:
        y_projection = ax.scatter(
            x_lim[1], data[y], data[z], c="lightgray", depthshade=False,
            **size_kwargs,
            zorder=-100,
        )

    if anchored is True:
        for _, point in data.iterrows():
            ax.plot([point[x], point[x]], [point[y], point[y]], [0, point[z]],
                    color='grey', linestyle='dashed', zorder=-10, linewidth=1)

    scatter = ax.scatter(
        data[x],
        data[y],
        data[z],
        depthshade=False,
        **color_kwargs,
        **size_kwargs,
        zdir="z",
        zorder=10,
    )

    if c is not None:
        cbar = fig.colorbar(scatter)

    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_zlabel(z)

    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.set_zlim(z_lim)

    ax.set_xticks([x for x in np.arange(lim_data[x].min(), lim_data[x].max() + 1, lim_data[x].max() // 4)])
    ax.set_yticks([x for x in np.arange(lim_data[y].min(), lim_data[y].max() + 1, lim_data[y].max() // 4)])

    return fig, ax



"""bottom = np.zeros_like(df_rnk["Codon1"])

sc = ax.scatter(df_rnk["Codon1"], df_rnk["Codon2"], df_val[column], c=df_val[column], cmap="YlGnBu",
                norm=colors.PowerNorm(gamma=1.5))
cbar = fig.colorbar(sc)
ax.set_xlabel("Codon1")
ax.set_ylabel("Codon2")
ax.set_zlabel("Counts")
plt.savefig(image_stream, format="png")
plt.close()
img_base64 = base64.b64encode(image_stream.getvalue())"""