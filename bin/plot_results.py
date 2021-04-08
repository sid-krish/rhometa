#! /usr/bin/env python
import sys
import seaborn as sns
import pandas as pd


collected_results = sys.argv[1]

# collected_results = "../Output/Results/collected_results.csv"

df = pd.read_csv(collected_results)

df["sample_size,genome_size"] = df.apply(lambda x: str(x["sample_size"])+','+str(x["genome_size"]), axis="columns")

df = df.sort_values(by=["sample_size,genome_size"])

df["abs_max_lk"] = df["max_lk"].abs()

sns.set_theme(style="whitegrid", palette="Blues_r")
# Flare and crest palettes are nice

g = sns.FacetGrid(df, col="rho")
g.map_dataframe(sns.barplot, x="max_rho", y="sample_size,genome_size")
g.set_axis_labels("estimated_rho", "sample_size,genome_size")

g.fig.suptitle('Custom Recombination Rate Estimator : rho vs estimated_rho', y=1.05)
g.savefig("rho_comparision.png", dpi=500)


g2 = sns.FacetGrid(df, col="rho")
g2.map_dataframe(sns.barplot, x="abs_max_lk", y="sample_size,genome_size")
g2.set_axis_labels("abs(max_lk)", "sample_size,genome_size")

g2.fig.suptitle('Custom Recombination Rate Estimator : rho vs max_likelihood', y=1.05)
g2.savefig("max_lk_comparision.png", dpi=500)