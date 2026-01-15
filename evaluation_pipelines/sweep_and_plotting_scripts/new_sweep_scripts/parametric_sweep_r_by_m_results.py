#!/usr/bin/env python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

theta_results_df = pd.read_csv("/home/sid/Github/rhometa/collected_theta_results_subsample.csv")
rho_results_df = pd.read_csv("/home/sid/Github/rhometa/collected_rho_results_subsample.csv")

# Merge theta and rho results
merged_df = pd.merge(
    rho_results_df,
    theta_results_df,
    left_on=["rho_sim","theta_sim","sample_size_sim","depth_sim","genome_size_sim","seed_sim"],
    right_on=["rho_sim","theta_sim","sample_size_sim","depth_sim","genome_size_sim","seed_sim"]
)

# Drop unnecessary columns
# merged_df.drop(
#     columns=[],
#     inplace=True,
# )


# Compute final rho per site by theta per site ratio
merged_df["rho_ps_by_theta_ps"] = (
    merged_df["per_site_rho"] / merged_df["theta_per_site_mean_depth"]
)

# Compute final r/m ratio
# Based on the formula: r/m = ρ (per site)/θ (per site) * tract length * substitution probability.
# As outlined in https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004041
merged_df["r_by_m"] = (
    merged_df["rho_ps_by_theta_ps"]
    * merged_df["recom_tract_len"]
    * merged_df["subst_probability"]
)


# Save final dataset
merged_df.to_csv(
    "/home/sid/Github/rhometa/collected_final_results_subsample.csv", index=False
)

ax = sns.catplot(data=merged_df, x="scaled_rho_sim", y="r_by_m", hue="sample_size_sim",
                     col="depth_sim", col_wrap=2, sharex=True, sharey=True, palette="Reds", kind="box",linewidth=0.2)

ax.set(xlabel="Simulated \u03C1", ylabel="r/m")

plt.savefig("rhometa_r_by_m_results.png", dpi=500)

#### deviation plot

merged_df["deviation"] = (merged_df['r_by_m'] - merged_df["scaled_rho_sim"]) / merged_df["scaled_rho_sim"]

ax2 = sns.catplot(data=merged_df, x="scaled_rho_sim", y="deviation", hue="sample_size_sim",
                    col="depth_sim", col_wrap=2, sharex=True, sharey=True, kind="box", palette="RdYlGn",
                    linewidth=0.2)

ax2.set(ylim=(-1, 1), xlabel="Simulated \u03C1", ylabel="Deviation")

ax2.map(plt.axhline, y=0, ls='--', color='g', linewidth=2)

ax2.savefig("rhometa_r_by_m_deviation_results.png", dpi=500)

