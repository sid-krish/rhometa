#!/usr/bin/env python
import sys

import pandas as pd

collected_results_csv = sys.argv[1]

# collected_results_csv = "/Users/Sid/Documents/GitHub/metagenome_pop_recom_rate/s_pnemonia/Recom_Est_Output/ERR025203_Aligned_final_results_summary.csv"

results_list = []
with open(collected_results_csv, 'r') as results:
    for line in results:
        if line.startswith("Theta_Est"):
            results_list.append(line.strip().split(","))


df_cols = "Theta_Est,mean,min,5%,25%,50%,75%,95%,max".split(",")
results_df = pd.DataFrame(results_list,columns=df_cols)

results_df.to_csv("theta_results.csv", index=False)

