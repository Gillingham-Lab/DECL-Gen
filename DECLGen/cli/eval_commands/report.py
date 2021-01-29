import os
import pandas as pd
import argh
import io
import base64
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
import matplotlib.colors as colors
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from DECLGen import Runtime
from DECLGen.report import HTMLReport, TextReport
from DECLGen.exceptions import DECLException, EvaluationFileDoesNotExist
from DECLGen.evaluation.evaluator import Evaluator
from DECLGen.cli.helpers import ProgressBar
from DECLGen.evaluation.plotter_report import ReportPlotter


def _parse_target_filename(save_as: str, result_filename: str):
    if save_as is None:
        save_as = os.path.join(os.path.dirname(result_filename), os.path.basename(result_filename) + ".report.html")
    else:
        if save_as[-5:] != ".html":
            save_as += ".html"

    save_as_radical = ".".join(save_as.split(".")[0:-1])

    return os.path.dirname(save_as), os.path.basename(save_as_radical), save_as, save_as_radical


@argh.arg("result", nargs="+")
@argh.arg("--save-as", type=str)
@argh.arg("--top")
@argh.arg("--plot-format", choices=["png", "svg", "pdf"])
def report(
        result: "Result file (result file from either from declEval extract or declEval compare)",
        save_as: "Report file name. Will be derived from result file name if not given." = None,
        top: "Number of top hits to include in the report" = 10,
        plot_format: "File formats of the plots" = "png",
):
    """
    Generates a summary of the given replicates.


    *Must be run from a declGen library folder (requires a valid decl_gen.data and library-properties.csv)*

    This command summarizes over the given replicates .csv files and reports on some parameters, such as coverage,
    known and unknown codons, and amount of counts.

    The result files are required for declGen compare.
    """
    r = Runtime()

    # Get library properties
    try:
        library_properties = r.storage.library.get_calculated_properties()
    except DECLException as e:
        r.error_exit(e)

    with ProgressBar(r.t, desc="Loading files") as progressBar:
        try:
            rr = Evaluator(*result, progress_bar=progressBar, plot_format=plot_format, properties=library_properties)
        except EvaluationFileDoesNotExist as e:
            print(r.t.val_bad(f"{e}"))

    filename_path, filename_base, filename_report, filename_radical = _parse_target_filename(save_as, result[0])

    # Save csv file
    try:
        rr.save_df(filename_radical + "-valid.csv")
    except FileNotFoundError:
        print(r.t.val_bad("It was not possible to save the valid files in {filename_radical}"))

    # Plot files
    with HTMLReport(filename_report, filename_base) as html_report:
        html_report.add_stats([
            (
                "Replicates",
                f"{rr.replicates:>10d}"
            ),
            (
                "Replicate files",
                ", ".join([x[2:] for x in rr.normalized_column_names]),
            ),
            (
                "Total counts",
                f"{rr.all_counts:>10.0f}"
            ),
            (
                "All codon combinations found",
                f"{rr.unique_codons:>10.0f}"
            ),
            (
                "Existing codons",
                f"{rr.valid_codons:>10d} "
                f"({rr.valid_codons / rr.unique_codons * 100:.1f}% of all combinations, "
                f"{rr.valid_counts / rr.all_counts * 100:.1f}% of all counts)"
            ),
            (
                "Mutated codons",
                f"{rr.unique_codons - rr.valid_codons:>10d} "
                f"({100 - rr.valid_codons / rr.unique_codons * 100:.1f}% of all combinations, "
                f"{100 - rr.valid_counts / rr.all_counts * 100:.1f}% of all counts)"
            ),
            (
                "Library size",
                f"{rr.library_size:>10d}"
            ),
            (
                "Structures not found",
                f"{rr.library_size - rr.valid_codons:>10d} "
                f"({(rr.library_size - rr.valid_codons) / rr.library_size * 100:.2f}% of library)"
            ),
            # (
            #    "Expected missing",
            #    f"{rr.expected_missing:>10.0f}"
            # ),
            (
                "Expected mean coverage",
                f"{rr.all_counts / rr.library_size / rr.replicates:>10.1f}X"
            ),
            (
                "Actual mean coverage",
                f"{rr.valid_counts / rr.library_size / rr.replicates:>10.1f}X"
            ),
        ])

        plotter = ReportPlotter(rr)

        with ProgressBar(r.t, desc="Drawing") as pb:
            # Plotting some meta information
            #html_report.append_plot("Replication scatter", rr.replicates_scatter(), rr._get_encoding())
            html_report.append_plotly(*plotter.cross_correlation())
            pb.update(0.05)

            html_report.append_plotly(*plotter.count_scatters()) if len(rr.codon_rank_columns) == 3 else None
            pb.update(0.1)

            #html_report.append_plot("Averaged histogram", rr.count_histogram(), rr._get_encoding())
            pb.update(0.2)

            html_report.append_plotly(*plotter.hits_preview()) if 3 >= len(rr.codon_rank_columns) >= 2 else None
            pb.update(0.5)

            """if rr.diversity_elements_count == 2:
                html_report.append_plot("3D Scatter", rr.two_element_3d_scatter(), rr._get_encoding())
                pb.update(0.3)

                html_report.append_plot(
                    "3D Scatter (top 200 hits)",
                    rr.two_element_3d_scatter(top_hits=200, project_on_z_plane=True, anchored=True),
                    rr._get_encoding())
                pb.update(0.4)
            elif rr.diversity_elements_count == 3:
                html_report.append_plot(
                    "3D Scatter (top 1000 hits)",
                    rr.three_element_3d_scatter(top_hits=1000, project_on_z_plane=True, anchored=True),
                    rr._get_encoding())
                pb.update(0.3)

                pb.update(0.4)"""

            i = 0
            for hit in rr.iter_hits(top):
                html_report.add_entry(hit.rank, hit.counts, hit.codons, hit.smiles, hit.draw())
                i += 1

                pb.update(0.5 + 0.5 * (i / top))

        text_report = TextReport.copyFromOther(html_report)

    html_report.save()

    print(text_report)
