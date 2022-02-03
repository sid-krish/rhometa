#!/usr/bin/env python
import sys
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt

if __name__ == '__main__':

    final_results = sys.argv[1]
    # final_results = "../Output/rho_20_sam_20_gen_10000/rho_20_sam_20_gen_10000_final_results.csv"

    df = pd.read_csv(final_results)

    df.drop(columns=["bootstrap_sample"], inplace=True)

    df.rename(columns={"rho": "rho_values_evaluated", "likelihood_sums": "log_likelihood_sums"}, inplace=True)


    # plot
    sns.set_theme(style="darkgrid")
    # Passing the entire dataset in long-form mode will aggregate over repeated values to show the mean and 95% confidence interval
    g = sns.lineplot(data=df,x="rho_values_evaluated", y="log_likelihood_sums").set_title("Custom Estimator Metagenomic Reads")
    # plt.xlim(0, 3)
    # bbox_inches="tight", prevents axis labels from going out of frame
    g.figure.savefig("final_results_plot.png", bbox_inches="tight", dpi=500)