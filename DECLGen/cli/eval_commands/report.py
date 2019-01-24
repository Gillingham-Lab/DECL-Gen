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
from DECLGen.exceptions import DECLException
from DECLGen.evaluation.evaluator import Evaluator
from DECLGen.cli.helpers import ProgressBar


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
def report(
        result: "Result file (result file from either from declEval extract or declEval compare)",
        save_as: "Report file name. Will be derived from result file name if not given." = None,
        top: "Number of top hits to include in the report" = 10,
):
    r = Runtime()

    # Get library properties
    try:
        library_properties = r.storage.library.get_calculated_properties()
    except DECLException as e:
        r.error_exit(e)

    with ProgressBar(r.t, desc="Loading files") as progressBar:
        rr = Evaluator(*result, progressBar=progressBar)

    rr.merge_properties(library_properties)

    filename_path, filename_base, filename_report, filename_radical = _parse_target_filename(save_as, result[0])

    with HTMLReport(filename_report, filename_base) as html_report:
        html_report.add_stats("Replicates", "{:>10d}".format(rr.replicates))
        html_report.add_stats("Total counts", "{:>10.0f}".format(rr.all_counts))
        html_report.add_stats("All codon combinations found", "{:>10d}".format(rr.unique_codons))
        html_report.add_stats("Valid combinations",
                              "{:>10d} ({:.1f}% of all combinations, {:.1f}% of all counts)".format(
                                  rr.valid_codons,
                                  rr.valid_codons / rr.unique_codons * 100,
                                  rr.valid_counts / rr.all_counts * 100
                                )
                              )
        html_report.add_stats("Invalid combinations",
                              "{:>10d} ({:.1f}% of all combinations, {:.1f}% of all counts)".format(
                                  rr.unique_codons - rr.valid_codons,
                                  100 - rr.valid_codons / rr.unique_codons * 100,
                                  100 - rr.valid_counts / rr.all_counts * 100
                                )
                              )

        html_report.add_stats("Library size", "{:>10d}".format(rr.library_size))
        html_report.add_stats("Structures not found",
                              "{:>10d} ({:.2f}% of library)".format(
                                  rr.library_size - rr.valid_codons,
                                  (rr.library_size - rr.valid_codons) / rr.library_size * 100
                                )
                              )
        html_report.add_stats("Expected average coverage", "{:>10.1f}X".format(rr.all_counts / rr.library_size / rr.replicates))
        html_report.add_stats("Actual average coverage", "{:>10.1f}X".format(rr.valid_counts / rr.library_size / rr.replicates))

        pb = ProgressBar(r.t, desc="Plotting values")
        pb.start()

        # Plotting some meta information
        html_report.append_plot("Replication scatter", rr.replicates_scatter())
        pb.update(0.1)

        html_report.append_plot("Averaged histogram", rr.count_histogram())
        pb.update(0.2)

        if rr.diversity_elements_count == 2:
            html_report.append_plot("3D Scatter", rr.two_element_3d_scatter())
            pb.update(0.3)

            html_report.append_plot("3D Scatter (top 200 hits)", rr.two_element_3d_scatter(top_hits=200, project_on_z_plane=True, anchored=True))
            pb.update(0.4)

        i = 0
        for hit in rr.iter_hits(top):
            html_report.add_entry(hit.rank, hit.counts, hit.codons, hit.smiles, hit.draw())
            i += 1

            pb.update(0.5 + 0.5*(i/top))

        pb.finish()


        text_report = TextReport.copyFromOther(html_report)

    print(text_report)

    #
    """res = pd.read_csv(result, sep="\t")
    if "CountRel" in res.columns:
        column = "CountRel"
    elif "Count" in res.columns:
        column = "Count"
    else:
        print("Cannot detect format of {}. Are you sure it was created by declEval extract or devlEval compare?".format(res))
        exit(1)
    if result2 is not None:
        # Load second result
        res2 = pd.read_csv(result2, sep="\t")
        res = res.merge(res2, on="Codon-Combination", how="outer", suffixes=["L", "R"])

        res = res.fillna(0)
        res[column] = (res[column + "L"] + res[column + "R"]) / 2

    # Merge with library
    res = res.sort_values(by=column, ascending=False, kind="mergesort")
    library = pd.read_csv("library-properties.csv", sep=",")

    merged = res.merge(library, on="Codon-Combination", how="outer")
    merged = merged.sort_values(by=column, ascending=False, kind="mergesort")
    merged_reduced = merged.dropna().reset_index(drop=True)

    merged2 = library.merge(res, on="Codon-Combination", how="left")
    merged2 = merged2.fillna(0)
    merged2 = merged2.sort_values(by=column, ascending=False, kind="mergesort")

    # Normalize reduced again
    the_mean = merged_reduced[column].mean()
    merged_reduced[column] = merged_reduced[column] / the_mean * 100
    merged2[column] = merged2[column] / the_mean * 100

    # If second sample is given, spread scatter
    if result2 is not None:
        lin_reggr = stats.linregress(merged_reduced[column + "L"], merged_reduced[column + "R"])

        with io.BytesIO() as image_stream:
            ax = merged_reduced.plot.scatter(column + "L", column + "R")
            ax.set_xlabel("Counts of " + result)
            ax.set_ylabel("Counts of " + result2)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(0.1, 10e4)
            ax.set_ylim(0.1, 10e4)
            x_regress = np.logspace(0, 4, 50)
            ax.plot(x_regress, lin_reggr[0] * x_regress + lin_reggr[1], "g-")
            plt.savefig(image_stream, format="png")
            img_base64 = base64.b64encode(image_stream.getvalue())
            plt.close()

        with io.BytesIO() as image_stream:
            ax = merged_reduced.plot.scatter(column + "L", column + "R")
            ax.set_xlabel("Counts of " + result)
            ax.set_ylabel("Counts of " + result2)
            plt.savefig(image_stream, format="png")
            img_base64_2 = base64.b64encode(image_stream.getvalue())
            plt.close()

        html_report.append_plot("Two sample scatter", img_base64)
        html_report.append_plot("Two sample scatter (no log)", img_base64_2)
        html_report.add_stats("Two-sample linear regression", "{:>10.5f}".format(lin_reggr[2]))
        html_report.add_stats("Two-sample linear p-value", "{:>10.5f}".format(lin_reggr[3]))

    # Histogram
    with io.BytesIO() as image_stream:
        ax = merged_reduced[column].plot.hist(bins=100)
        ax.set_yscale("log")
        plt.savefig(image_stream, format="svg")
        img_base64 = base64.b64encode(image_stream.getvalue())
        plt.close()
    html_report.append_plot("Count distribution", img_base64,  encoding="image/svg+xml;base64")

    # 3D-Scatter of top-1000 hits if 3-codon
    a = merged2["Codon-Combination"].iloc[0].split("-")
    if len(a) == 3:
        # Overview-Plot
        df_tops_val = merged_reduced.nlargest(1000, column)
        df_tops_val_codons = pd.DataFrame(df_tops_val["Codon-Combination"].str.split('-',2).tolist(), columns = ['Codon1','Codon2','Codon3'])
        df_tops_val = df_tops_val.join(df_tops_val_codons)
        df_tops_val = df_tops_val.sort_values(by="Codon-Combination")
        df_tops_rnk = df_tops_val.rank(method="dense")
        c = df_tops_val[column]

        with io.BytesIO() as image_stream:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            sc = ax.scatter(df_tops_rnk["Codon1"], df_tops_rnk["Codon2"], df_tops_rnk["Codon3"], c=c,
                       cmap="YlGnBu", s=(df_tops_val[column]/df_tops_val[column].max())**2*100, alpha=.8,
                       #norm=colors.LogNorm(vmin=c.min(), vmax=c.max()))
                       norm=colors.PowerNorm(gamma=0.2))
                       #     )
            cbar = fig.colorbar(sc)
            ax.set_xlabel("Codon1")
            ax.set_ylabel("Codon2")
            ax.set_zlabel("Codon3")
            plt.savefig(image_stream, format="svg")
            plt.close()
            img_base64 = base64.b64encode(image_stream.getvalue())

        html_report.append_plot("3D-Scatter of Codon/Codon/Codon", img_base64, encoding="image/svg+xml;base64")

        # Codon-1 split
        with io.BytesIO() as image_stream:

            merged_reduced_annot = pd.DataFrame(merged_reduced["Codon-Combination"].str.split('-', 2).tolist(),
                                              columns=['Codon1', 'Codon2', 'Codon3'])
            merged_reduced_annot = merged_reduced.join(merged_reduced_annot)
            merged_reduced_annot = merged_reduced_annot.sort_values(by="Codon-Combination")

            countMax = merged_reduced_annot[column].max()
            grouped = merged_reduced_annot.groupby("Codon1")
            C = 3
            R = np.ceil(len(grouped)/3)
            P = 0

            fig = plt.figure(figsize=(C*5,R*3))
            for name, group in grouped:
                P += 1
                countMax = group[column].max()

                ax = fig.add_subplot(R, C, P, projection='3d')

                c = group.rank(method="dense")
                #c[group[column] < 200] = np.nan
                c[group[column] < group.nlargest(500, column).min()[column]] = np.nan

                group = group[c.isna() == False]
                c = c[c.isna() == False]

                sc= ax.scatter(c["Codon2"], c["Codon3"], group[column], c=group[column], cmap="YlGnBu",
                               norm=colors.PowerNorm(gamma=0.2))
                cbar = fig.colorbar(sc)
                ax.set_xlabel("Codon2")
                ax.set_ylabel("Codon3")
                ax.set_zlabel("Count")
                ax.set_xlim(0, c["Codon2"].max())
                ax.set_ylim(0, c["Codon3"].max())
                #ax.set_zlim(0, countMax)
                ax.set_title(name)

            plt.savefig(image_stream, format="svg")
            plt.close()
            img_base64 = base64.b64encode(image_stream.getvalue())

        html_report.append_plot("Codon1 Split scatters", img_base64, encoding="image/svg+xml;base64")

    elif len(a) == 2:
        df_val = merged_reduced.copy()
        df_val_codons = pd.DataFrame(df_val["Codon-Combination"].str.split('-',1).tolist(), columns = ['Codon1','Codon2'])
        df_val = df_val.join(df_val_codons)
        df_val = df_val.sort_values(by="Codon-Combination")
        df_rnk = df_val.rank(method="dense")
        c = df_val[column]

        with io.BytesIO() as image_stream:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            bottom = np.zeros_like(df_rnk["Codon1"])

            sc = ax.scatter(df_rnk["Codon1"], df_rnk["Codon2"], df_val[column], c=df_val[column], cmap="YlGnBu", norm=colors.PowerNorm(gamma=0.5))
            cbar = fig.colorbar(sc)
            ax.set_xlabel("Codon1")
            ax.set_ylabel("Codon2")
            ax.set_zlabel("Counts")
            plt.savefig(image_stream, format="svg")
            plt.close()
            img_base64 = base64.b64encode(image_stream.getvalue())

        html_report.append_plot("3D-Bar of Codon/Codon/Count", img_base64,  encoding="image/svg+xml;base64")

        df_rnk[df_val[column] < df_val.nlargest(200, column).min()[column]] = np.nan

        with io.BytesIO() as image_stream:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            bottom = np.zeros_like(df_rnk["Codon1"])

            sc = ax.scatter(df_rnk["Codon1"], df_rnk["Codon2"], df_val[column], c=df_val[column], cmap="YlGnBu", norm=colors.PowerNorm(gamma=0.5))
            cbar = fig.colorbar(sc)
            ax.set_xlabel("Codon1")
            ax.set_ylabel("Codon2")
            ax.set_zlabel("Counts")
            plt.savefig(image_stream, format="svg")
            plt.close()
            img_base64 = base64.b64encode(image_stream.getvalue())

        html_report.append_plot("3D-Bar of Codon/Codon/Count (reduced, top 200)", img_base64,  encoding="image/svg+xml;base64")

    #print(len(res), len(merged), len(merged[merged["Canonical smiles"].isna()]), len(merged[merged[column].isna()]))

    unknown_structures_count = merged[merged["Canonical smiles"].isna()][column].sum()
    all_structures_count = merged[column].sum()

    total_structures = len(res)
    unknown_structures = len(merged[merged["Canonical smiles"].isna()])
    uncovered_structures = len(merged[merged[column].isna()])

    #print(merged_reduced)
    valid_structures = len(merged_reduced)

    html_report.add_stats("Total count", "{:>10.0f}".format(all_structures_count))
    html_report.add_stats("Unknown count", "{:>10.0f}".format(unknown_structures_count))

    html_report.add_stats("Total structures", "{:>10d}".format(total_structures))
    html_report.add_stats("Unknown structures", "{:>10d} ({:.1f}%, {:.1f}% of counts)".format(unknown_structures, unknown_structures / total_structures * 100, unknown_structures_count/all_structures_count*100))
    html_report.add_stats("Library structures not in dataset", "{:>10d} ({:.1f}%)".format(uncovered_structures, uncovered_structures/len(library)*100))
    html_report.add_stats("Valid structures", "{:>10d} ({:.1f}% of results, {:.1f}% of all structures)".format(valid_structures,
                                                                                          valid_structures / total_structures * 100,
                                                                                          valid_structures / len(library) * 100))

    text_report.stats = html_report.stats

    top_hits: pd.DataFrame = merged_reduced[:top].reset_index(drop=True)

    for i in range(len(top_hits)):
        smiles = top_hits["Canonical smiles"][i]
        relCount = top_hits[column][i]
        codonComb  = top_hits["Codon-Combination"][i]
        m = Chem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(m)

        with io.BytesIO() as image_stream:
            img = Draw.MolToImage(m, size=(300, 300), kekulize=True, wedgeBonds=True)
            img.save(image_stream, format="png")

            img_base64 = base64.b64encode(image_stream.getvalue())

        html_report.add_entry(i+1, relCount, codonComb, smiles, img_base64)

    html_report.save()

    print(text_report)

    merged_reduced.to_csv(".".join(save_as.split(".")[0:-1]) + "-valid.csv", sep="\t", index=False)
    """