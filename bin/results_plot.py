#!/usr/bin/env python
import sys

import pandas as pd
import seaborn as sns

if __name__ == '__main__':
    results = sys.argv[1]

    df = pd.read_csv(results)

    df.rename(columns={"rho": "rho_values_evaluated"}, inplace=True)

    # plot
    sns.set_theme(style="darkgrid")
    g = sns.lineplot(data=df, x="rho_values_evaluated", y="log_likelihood_sum").set_title("rho v. log likelihood sum")

    # bbox_inches="tight", prevents axis labels from going out of frame
    g.figure.savefig("results_plot.png", bbox_inches="tight", dpi=500)
