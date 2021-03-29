#!/usr/bin/env python
import sys

import pandas as pd

collected_likelihoods_csv = sys.argv[1]
theta = sys.argv[2]

# collected_likelihoods_csv = "../Output/Results/collected_likelihoods.csv"
# theta = 0.01

df = pd.read_csv(collected_likelihoods_csv)

# If multiple total_of_log_likelihoods are same, sort such that the ones connected to smallest rho_for_estimator
# are first
df.sort_values(by=["total_of_log_likelihood", "rho_for_estimator"], ascending=[False, True], inplace=True)

# df.reset_index(inplace=True, drop=True)

max_rho = df.iloc[0][0]
max_lk = df.iloc[0][1]

# max_rho = round(df.iloc[0][0], 2)
# max_lk = round(df.iloc[0][1], 2) # rounding to make it look cleaner

df.sort_values(by="rho_for_estimator", inplace=True)
# df = df.round(2) # just so it looks cleaner
# with open(theta_txt, 'r') as file:
#     theta = file.readline()

with open("final_results.txt", 'w') as file:  # open in write mode (create new file)
    file.write(f"Custom Pairwise Recombination Rate Estimator\n")
    file.write(f"\n")

    file.write(f'-' * 45 + '\n')
    file.write(f"Theta estimate = {theta}\n")
    file.write(f"Max rho = {max_rho}, Max lk = {max_lk}\n")
    file.write(f'-' * 45 + '\n')

    file.write(f"\n")
    df.to_csv(file, mode = 'a', index = None)  # open in append mode (add to new file)
