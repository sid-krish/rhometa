#! /usr/bin/env python
import seaborn as sns
import pandas as pd

# collected_results = sys.argv[1]

collected_results = "collected_aggregate_results_filtered.csv"

df = pd.read_csv(collected_results)

df.drop(columns=["run"], inplace=True)

# plot
sns.set_theme(style="darkgrid")
# Passing the entire dataset in long-form mode will aggregate over repeated values to show the mean and 95% confidence interval
g = sns.barplot(data=df,x="collection_date", y="population_recombination_rate", color='#4878d0', capsize=.2)
g.set_xticklabels(g.get_xticklabels(),rotation=90)
# bbox_inches="tight", prevents axis labels from going out of frame
g.figure.savefig("time_series_plot.png", bbox_inches="tight", dpi=500)
