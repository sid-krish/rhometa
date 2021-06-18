import pandas as pd
import seaborn as sns

simulations = "simulation_2.csv"

df = pd.read_csv(simulations)

sns.set(style="darkgrid")
g = sns.FacetGrid(df, col="type", sharey=False, hue="type", height=3.5, aspect=1.5)
g.map_dataframe(sns.lineplot, x="rho_for_estimator", y="total_of_log_likelihood")
g.set_titles('{col_name}')
g.set_axis_labels("rho_values_evaluated", "log_likelihood_sums")

g.tight_layout()
g.savefig("simulation_2.png", dpi=500)
