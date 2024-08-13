import pandas as pd

# Logic
# df = pd.read_csv("/Users/sid/Documents/GitHub/ch2/test2.csv")

# test = df.groupby(["c","d","a"])['b'].mean().reset_index(name="mean")

# merged = df.merge(test, how='left',left_on=["c","d","a"], right_on=["c","d","a"])

# merged.drop(columns=["b"], inplace=True)

# final=merged.drop_duplicates().copy()

# final.drop(columns=["a"], inplace=True)

df = pd.read_csv("")

gdf = df.groupby(["identifier", "reference"])["rho"].mean().reset_index(name="mean_rho")

merged = df.merge(gdf, how='left',left_on=["identifier", "reference"], right_on=["identifier", "reference"])

merged["mean_rho_per_site"] = merged["mean_rho"] / (2 * merged["recom_tract_len"]) # To get c divide by 2t. 2N_e is still present unchanged

merged.drop(columns=["seed", "rho", "rho_per_site", "log_likelihood_sum", "recom_tract_len"], inplace=True)

final = merged.drop_duplicates().copy()

final['rho_ps_by_theta_ps'] = final["mean_rho_per_site"] / final["theta_per_site_mean_depth"]

final.to_csv("collected_merged_mean_rho_rps_by_tps.csv", index=False)
