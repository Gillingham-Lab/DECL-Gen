import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statsmodels as sm
from statsmodels.graphics.gofplots import qqplot_2samples
import os
import argh
from colorama import init, Fore, Style
from scipy.stats import f_oneway
from statsmodels.stats.multitest import multipletests
from DECLGen.report import HTMLReport, TextReport
from DECLGen.evaluation.plotter_compare import ComparePlotter
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit import RDLogger
import io
import base64


lg = RDLogger.logger().setLevel(RDLogger.ERROR)


def draw(row) -> bytes:
    m = Chem.MolFromSmiles(row["Canonical smiles"])
    AllChem.Compute2DCoords(m)

    with io.BytesIO() as image_stream:
        img = Draw.MolToImage(m, size=(300, 300), kekulize=True, wedgeBonds=True)
        img.save(image_stream, format="png")

        img_base64 = base64.b64encode(image_stream.getvalue())

    return img_base64


@argh.arg("--top")
def compare(
        signal_file: "Signal file",
        noise_file: "Noise file",
        result_file: "Result file",
        top: "Number of top hits to draw in Summary" = 10,
):
    init()

    if not os.path.exists(signal_file):
        print(f"{Fore.RED}Signal file {signal_file} was not found.{Fore.RESET}")
        exit(1)

    if not os.path.exists(noise_file):
        print(f"{Fore.RED}Noise file {noise_file} was not found.{Fore.RESET}")
        exit(1)

    signal = pd.read_csv(signal_file, sep="\t", index_col="Codon-Combination")
    noise = pd.read_csv(noise_file, sep="\t", index_col="Codon-Combination")

    if signal.shape[1] != noise.shape[1]:
        print(f"{Fore.RED}Noise and signal file have an unequal amount of columns. Are you sure they belong to the same library?{Fore.RESET}")
        exit(1)

    # Calculate the amount of codons
    number_of_codons = len(signal.index[0].split("-"))

    # Get the column names with normalized counts
    # For signal
    signal_old_count_columns = [x for x in signal.columns if x.startswith("N_") and x.endswith(".csv")]
    signal_count_columns = [f"Signal_{x+1}" for x in range(len(signal_old_count_columns))]
    signal.rename(columns={x: y for x, y in zip(signal_old_count_columns, signal_count_columns)}, inplace=True)
    # For noise
    noise_old_count_columns = [x for x in noise.columns if x.startswith("N_") and x.endswith(".csv")]
    noise_count_columns = [f"Noise_{x+1}" for x in range(len(noise_old_count_columns))]
    noise.rename(columns={x: y for x, y in zip(noise_old_count_columns, noise_count_columns)}, inplace=True)

    # Shift data
    signal[signal_count_columns] += signal[signal_count_columns][signal[signal_count_columns] > 0].min()
    noise[noise_count_columns] += noise[noise_count_columns][noise[noise_count_columns] > 0].min()

    # Get the interesting columns
    signal_good_columns = [*signal.columns[:number_of_codons], "Canonical smiles", *signal_count_columns, "mean"]

    # Join both dataframes
    data = signal[signal_good_columns].join(noise[[*noise_count_columns, "mean"]], rsuffix="_signal", lsuffix="_noise")

    # Reset index
    data.reset_index(inplace=True)

    # Calculate log-fold enrichment, pValue and qValue

    # Mean
    data["mean_signal"] = data[signal_count_columns].mean(axis=1)
    data["mean_noise"] = data[noise_count_columns].mean(axis=1)

    # log-fold enrichment (as log 2)
    data["log-Fold"] = np.log2(data["mean_signal"] / data["mean_noise"])

    # p-Value
    _, data["pValue"] = f_oneway(data[signal_count_columns], data[noise_count_columns], axis=1)

    # Adjust p-Value with Benjamini-Hochberg
    # multipletests returns reject, pvalue_corrected, alphacSidak, alphacBonf
    # Only the corrected p_value is of interest here.
    data["qValue"] = multipletests(data["pValue"], alpha=0.10, method="fdr_bh")[1]

    # Score
    data["score"] = -np.log10(data["qValue"]) * data["log-Fold"]
    data.sort_values(by="score", ascending=False, inplace=True)

    # Reporting time

    with HTMLReport(result_file + ".html", os.path.basename(result_file)) as html_report:
        plotter = ComparePlotter(data)

        html_report.append_plotly(*plotter.enrichment_plot())
        html_report.append_plotly(*plotter.vulcano_plot_p())
        html_report.append_plotly(*plotter.vulcano_plot_q())

        i = 0
        for rowkey, hit in data.head(top).iterrows():
            html_report.add_entry(i, hit["score"], hit["Codon-Combination"], hit["Canonical smiles"], draw(hit))
            i += 1

    html_report.save()

    # Save CSV
    data.to_csv(result_file + ".csv", sep="\t")

    print(f"{Fore.GREEN}Done{Fore.RESET}")
