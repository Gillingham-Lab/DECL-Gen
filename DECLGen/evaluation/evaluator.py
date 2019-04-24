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
from DECLGen.cli.helpers import ProgressBar

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
    """
    Class to evaluate NGS data
    """
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
    progress_bar: ProgressBar


    def __init__(self, *results, progress_bar: ProgressBar = None, plot_format: str = "png", properties: pd.DataFrame):
        """
        :param results: A list of 1 or more csv files generated from declEval extract. Multiple files equal replicates and the counts will get averaged.
        :param progress_bar: progress bar for cli interaction.
        :param plot_format: Plot format, either png or svg.
        :param properties: A pandas dataframe containing all library members and their properties. Should be a file generated by declGen lib-gen.
        """
        if (plot_format not in ["png", "svg"]):
            raise Exception("plot_format must be png or svg.")

        self.library_size = len(properties)
        self.count_table = None
        self.progress_bar = progress_bar
        self.load_counts(results)
        self.plot_format = plot_format

        self.merge_properties(properties)

    @property
    def expected_missing(self) -> int:
        """
        Returns the amount of codons that are expected to be missing from the result assuming a equal distribution of all codons.
        :return: Number of codons that are expected to be missing from the dataset.
        """
        N = self.library_size
        p = []

        for c in self.count_columns:
            x = 1 - (1 - 1/N)**self.count_table_purified[c].sum()
            p.append(x)

        P = N
        for i in p:
            P *= i
        return self.library_size - P

    def _get_encoding(self) -> str:
        """
        Returns the encoding string depending on the type for inline-representation.
        :return: The encoding for use in the HTML img tag.
        """
        if self.plot_format == "svg":
            return "image/svg+xml;base64"

        return "image/png;base64"

    def merge_properties(self, properties: pd.DataFrame):
        """
        Merges the properties generated by declGen lib-gen with the count files loaded in.
        :param properties:
        :return:
        """
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

        # We decode here the codon to generate a number used as the "codon rank".
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

        if self.progress_bar is not None:
            self.progress_bar.update(0.99)

    def save_df(self, filename):
        """
        Saves the dataframe to a file with standard parameters (sep=tab, index=False)
        :param filename:
        :return:
        """
        self.count_table_purified.to_csv(filename, sep="\t", index=False)

    def load_counts(self, results: Sequence[str]):
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

            # For every datafile given
            for filename in results:
                i+=1

                # Update the progressBar if given to show loading progress
                if self.progress_bar is not None:
                    self.progress_bar.update(i / (replicates + 1))

                # read in the data file
                df = pd.read_csv(filename, sep="\t")

                # Add a column name depending on the filename.
                count_column_names.append(os.path.basename(filename))

                # Determine the colum name to use (Count or CountRel) depending on which source has been used.
                # only detected once. Raises an exception if the column name was not found, or if two different file sources are mixed.
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
                # Of not, we must merge.
                else:
                    result = result.merge(df, on="Codon-Combination", how="outer")

                # Because there are multiple Count columns, we must enumerate them
                new_column = column + "{}".format(i)
                column_number = format(i)

                # Rename Count to CountX, where X is a integer.
                result.rename(columns={column: new_column}, inplace=True)

                # Calculate the enrichment factor and put it into a new column
                result["Enrichment{}".format(i)] = calculate_enrichment_factor(result[new_column], self.library_size, result[new_column].sum())

                # Append the new colum name to a list of count_columns, and append the number fo a list of column numbers.
                count_columns.append(new_column)
                column_numbers.append(column_number)

            # Fill empty fields with 0
            result.fillna(0, inplace=True)

            # Save mean and standard deviation
            result["MeanCounts"] = result[count_columns].mean(axis=1)
            result["StdDevCounts"] = result[count_columns].std(axis=1)

            result["MeanEnrichment"] = result[["Enrichment{}".format(i) for i in column_numbers]].mean(axis=1)
            result["StdEnrichment"] = result[["Enrichment{}".format(i) for i in column_numbers]].std(axis=1)

            # Recalculate the enrichment confidence boundaries of the mean enrichment. I'm not sure if this is valid, but its currently the best approach.
            lower, upper = estimate_enrichment_confidence_boundaries(result["MeanEnrichment"], 1, 1)
            result["MeanEnrichment_Lower"] = lower
            result["MeanEnrichment_Upper"] = upper

            # Calculate the Variance of the enrichment (var=std**2)
            result["VarEnrichment"] = result["StdEnrichment"]**2

            # Sort the table by MeanEnrichment, StdEnrichment
            result.sort_values(by=["MeanEnrichment", "StdEnrichment"], ascending=False, inplace=True)

        except FileNotFoundError:
            raise EvaluationFileDoesNotExist("Result file <{}> does not exist.".format(filename))

        # Save the results in this object.
        self.count_table = result
        self.replicates = replicates
        self.count_columns = count_columns
        self.count_column_names = count_column_names
        self.column_numbers = column_numbers

        # This here essentially tries to rename the enumarted count_columns and replaces it with the filename. Makes identification easier.
        self.count_table.rename(dict(zip(self.count_columns, self.count_column_names)), axis="columns", inplace=True)
        self.count_columns, self.count_column_names = self.count_column_names, self.count_columns

        # Statistics
        self.unique_codons = self.count_table["Codon-Combination"].nunique()
        self.all_counts = int(self.count_table[self.count_columns].sum().sum())

    def replicates_scatter(self, scale = None):
        """
        Generates a scatter plot of each replicate with each other and returns the base64 encoded image.
        :param scale:
        :return:
        """
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
        """
        Plots the variance of the sample versus the mean. Ideally, this is a linear relationship if the original sample was equally distributed.
        :return:
        """
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
        """
        Plots a histogram of the counts and returns the base64 encoded plot.
        :return:
        """
        with io.BytesIO() as image_stream:
            ax = sns.kdeplot(self.count_table_reduced["MeanEnrichment"], shade=True, legend="mean counts")

            plt.savefig(image_stream, format=self.plot_format)
            img_base64 = base64.b64encode(image_stream.getvalue())
            plt.close()

        return img_base64

    def two_element_3d_scatter(self,
        top_hits: Optional[int] = None,
        codon_x: int = 1, codon_y: int = 2,
        anchored: bool=False,
        project_on_z_plane: bool =False
    ) -> str:
        """
        Plots a 3D scatter diagram of codon_x and codon_y.
        :param top_hits: Number if hits to plot. Plots all if not given.
        :param codon_x: Codon for the x axis (1 is the first DE).
        :param codon_y: Codon for the y axis (1 is the first DE).
        :param anchored: If true, the plot will be a pin-needle plot to "anchor" the points down to the lowest xy plane.
        :param project_on_z_plane: If true, all points will be projected into the xy plane.
        :return: The base64 encoded plot.
        """
        codon_x = self.codon_rank_columns[codon_x - 1]
        codon_y = self.codon_rank_columns[codon_y - 1]

        # We must sort the dataframe by codon so that they overlay properly.
        df = self.count_table_reduced.sort_values(by=[codon_x, codon_y], ascending=[False, False])

        # Reduces the dataframe to the top N hits if given.
        if top_hits is not None:
            df_reduced = df.nlargest(top_hits, "MeanEnrichment")
        else:
            df_reduced = df

        # Create the plot
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

    def three_element_3d_scatter(self,
        top_hits: Optional[int] = None,
        codon_x: int = 1, codon_y: int = 2, codon_z: int = 3,
        anchored: bool=False,
        project_on_z_plane: bool=False
    ) -> str:
        """
        Plots a 3D scatter plot with 3 codons on the x, y and z axis. colour and size scaling is according to MeanEnrichment.
        :param top_hits:
        :param codon_x:
        :param codon_y:
        :param codon_z:
        :param anchored:
        :param project_on_z_plane:
        :return:
        """
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