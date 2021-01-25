import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statsmodels as sm
from statsmodels.graphics.gofplots import qqplot_2samples
import os


def compare(
        signal_file: "Signal file",
        noise_file: "Noise file",
        result_file: "Result file",
):
    """
    Compares a signal file to a noise file and reports relative enrichment.
    :param signal_file:
    :param noise_file:
    :return:
    """

    if len(signal_file.split(",")) > 0:
        signal_files = signal_file.split(",")
        signal = None
        files = len(signal_files)
        print("#Signal files", files)
        for signal_file in signal_files:
            signal_file = signal_file.strip()

            if signal is None:
                signal = pd.read_csv(signal_file, sep="\t", index_col=0)
                continue

            signal_b = pd.read_csv(signal_file, sep="\t", index_col=0)
            joined = signal.join(signal_b, how="outer", rsuffix="B")
            joined[joined.isna() == True] = 0

            joined["Count"] += joined["CountB"]
            signal = joined[["Count"]]

        signal = signal[["Count"]]/files
    else:
        signal = pd.read_csv(signal_file, sep="\t")

    if dummy_file is not None:
        if len(dummy_file.split(",")) > 0:
            dummy_files = dummy_file.split(",")
            dummy = None
            files = len(dummy_files)
            print("#Dummy files", files)
            for dummy_file in dummy_files:
                dummy_file = dummy_file.strip()

                if dummy is None:
                    dummy = pd.read_csv(dummy_file, sep="\t", index_col=0)
                    continue

                dummy_b = pd.read_csv(dummy_file, sep="\t", index_col=0)
                joined = dummy.join(dummy_b, how="outer", rsuffix="D")
                joined[joined.isna() == True] = 0

                joined["Count"] += joined["CountD"]
                dummy = joined[["Count"]]

            dummy = dummy[["Count"]] / files
        else:
            dummy = pd.read_csv(signal_file, sep="\t")

        print("Using dummy enrichment ...")
        signal = signal.join(dummy, how="outer", rsuffix="D", lsuffix="O")
        signal[joined.isna() == True] = 1
        signal["Count"] = signal["CountO"] / signal["CountD"]

    if len(noise_file.split(",")) > 0:
        noise_files = noise_file.split(",")
        noise = None
        files = len(noise_files)
        print("#Noise files", files)
        for noise_file in noise_files:
            if noise is None:
                noise = pd.read_csv(noise_file, sep="\t", index_col=0)
                continue

            noise_b = pd.read_csv(noise_file, sep="\t", index_col=0)
            joined = noise.join(noise_b, how="outer", rsuffix="B")
            joined[joined.isna() == True] = 0

            joined["Count"] += joined["CountB"]
            noise = joined[["Count"]]

        noise = noise[["Count"]] / files
    else:
        noise = pd.read_csv(noise_file, sep="\t")

    joined = signal.join(noise, how="outer", lsuffix="S", rsuffix="N").dropna()
    joined["CountRel"] = joined["CountS"] / joined["CountN"]

    x = np.zeros(100)
    y = np.zeros(100)
    for i in range(100):
        q = i/100
        x[i] = joined["CountS"].quantile(q)
        y[i] = joined["CountN"].quantile(q)

    qqplot_2samples(joined["CountS"], joined["CountN"])

    #plt.scatter(x, y)
    plt.savefig(result_file + ".png", dpi=300)
    joined.to_csv(result_file + ".csv", sep="\t", index=True)

    #"""