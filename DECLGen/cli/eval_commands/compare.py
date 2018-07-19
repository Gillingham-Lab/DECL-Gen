import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statsmodels as sm
from statsmodels.graphics.gofplots import qqplot_2samples

def compare(
        signal_file: "Signal file",
        noise_file: "Noise file",
        result_file: "Result file"
):
    """
    Compares a signal file to a noise file and reports relative enrichment.
    :param signal_file:
    :param noise_file:
    :return:
    """

    signal = pd.read_csv(signal_file, sep="\t")
    noise = pd.read_csv(noise_file, sep="\t")

    joined = signal.merge(noise, how="outer", on="Codon-Combination", suffixes=["S", "N"]).dropna()
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
    joined.to_csv(result_file + ".csv", sep="\t", index=False)