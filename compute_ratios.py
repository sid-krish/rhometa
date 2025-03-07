import pandas as pd

# Load data
df = pd.read_csv("collected_merged_results.csv")

# Compute mean values for rho, per_site_rho, and subst_probability
gdf = df.groupby(["identifier", "reference"]).agg(
    mean_rho=("rho", "mean"),
    mean_per_site_rho=("per_site_rho", "mean"),
    mean_subst_probability=("subst_probability", "mean")  # Taking mean of subst_probability
).reset_index()

# Merge mean values back into the original DataFrame
merged = df.merge(gdf, on=["identifier", "reference"], how="left")

# Drop unnecessary columns (keeping recom_tract_len since it's needed)
merged.drop(columns=["seed", "rho", "per_site_rho", "subst_probability", "log_likelihood_sum"], inplace=True)

# Remove duplicate rows based on remaining columns
final = merged.drop_duplicates().copy()

# Compute final rho per site by theta per site ratio
final["rho_ps_by_theta_ps"] = final["mean_per_site_rho"] / final["theta_per_site_mean_depth"]

# Compute final r/m ratio
# Based on the formula: r/m = ρ (per site)/θ (per site) * tract length * substitution probability. As outlined in https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004041
final["r_by_m"] = final["mean_per_site_rho"] / final["theta_per_site_mean_depth"] * final["recom_tract_len"] * final["mean_subst_probability"]

# Save final dataset
final.to_csv("collected_merged_final.csv", index=False)
