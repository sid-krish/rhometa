#!/usr/bin/env python
import sys

import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt

agg_results_csv = sys.argv[1]

agg_results_df = pd.read_csv(agg_results_csv)

# Plot results
sns.set_theme(style="white")

ax = sns.boxplot(data=agg_results_df, x="mean", palette="Accent")

ax.set(xlabel="Estimated \u03C1 (mean)")

ax.figure.savefig("agg_recom_results_plot.png", dpi=500)