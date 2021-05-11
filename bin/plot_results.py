#! /usr/bin/env python
import sys
import seaborn as sns
import pandas as pd


collected_results = sys.argv[1]

# collected_results = "../Output/Results/collected_results.csv"

df = pd.read_csv(collected_results)

df = df.sort_values(by=["depth"])

# df["abs_max_lk"] = df["max_lk"].abs()

sns.set_theme(style="whitegrid", palette="Blues_r")
# Flare and crest palettes are nice

g = sns.lineplot(x="depth", y="max_rho",data=df)
g.figure.savefig("depth_v_max_rho.png", dpi=500)