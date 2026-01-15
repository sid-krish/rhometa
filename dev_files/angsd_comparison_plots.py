#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# === Input files ===
rhometa_file = "/home/sid/Github/rhometa/original_theta_sweep.csv"
angsd_file   = "/home/sid/Github/rhometa/collected_theta_results_subsample.csv"

# Load data
df_rhometa = pd.read_csv(rhometa_file)
df_angsd   = pd.read_csv(angsd_file)

merge_cols = [
    "rho_sim","theta_sim","scaled_theta_sim",
    "sample_size_sim","depth_sim","genome_size_sim","seed_sim"
]

# Merge
df_merged = pd.merge(
    df_rhometa[merge_cols+["tps_mean_depth"]],
    df_angsd[merge_cols+["angsd_theta"]],
    on=merge_cols, how="inner"
)

# Reshape to long format for seaborn
df_long = pd.melt(
    df_merged,
    id_vars=["scaled_theta_sim"],
    value_vars=["tps_mean_depth","angsd_theta"],
    var_name="method",
    value_name="theta_est"
)

# Clean method labels
df_long["method"] = df_long["method"].map({
    "tps_mean_depth": "Rhometa (mean depth)",
    "angsd_theta": "ANGSD"
})

# Plot
plt.figure(figsize=(9,6))
sns.boxplot(data=df_long, x="scaled_theta_sim", y="theta_est", hue="method", width=0.6)

# Add truth lines
for t in sorted(df_long["scaled_theta_sim"].unique()):
    plt.axhline(y=t, linestyle="--", color="black")

plt.xlabel("Simulated θ")
plt.ylabel("Estimated θ")
plt.title("Rhometa vs ANGSD")
plt.legend(title="Method", loc="upper left")
plt.tight_layout()
plt.savefig("theta_comparison_seaborn.png", dpi=300)
plt.show()