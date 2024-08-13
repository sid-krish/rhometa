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

merged["mean_rho_per_site"] = merged["mean_rho"] / (2 * merged["recom_tract_len"]) 
# To get c divide by 2t. 2N_e is still present unchanged, so it would be 2Nec
# r in 2Ner is 2ct. Refer to eqn 4 in https://academic.oup.com/genetics/article/160/3/1231/6052507. The 1-e^(-d_{ij}/t) term dropping out (asymptoting to 1). so that leaves us with r_genome = 2ct.
# that factor of 2 in the per site recombination rate is because a gene conversion has two endpoints, but a crossing-over has only one. and the lookup table is for single crossover
# one way to take the intuition is to imagine a gene conversion event with length 2000bp, being analysed with 500nt long illumina fragments. 
# Because there are no read pairs that span 2000bp, the breakpoints show up separately in the data - it's a single event, but looks like two separate 
# recombination events because they are too far apart to draw a connection between them with short read data

merged.drop(columns=["seed", "rho", "rho_per_site", "log_likelihood_sum", "recom_tract_len"], inplace=True)

final = merged.drop_duplicates().copy()

final['rho_ps_by_theta_ps'] = final["mean_rho_per_site"] / final["theta_per_site_mean_depth"]

final.to_csv("collected_merged_mean_rho_rps_by_tps.csv", index=False)
