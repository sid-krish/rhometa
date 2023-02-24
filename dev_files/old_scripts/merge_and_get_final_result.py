#!/usr/bin/env python
import glob
import sys

import pandas as pd

theta = sys.argv[1]

# Script looks as current working directory
extension = 'txt'
all_filenames = [i for i in glob.glob('*.{}'.format(extension))]

all_filenames.sort(key=lambda x: int(x[:-4].split('_')[-1]))

# all_filenames = all_filenames[0:2] # for testing just a few files

# combine all files in the list
combined_csv = pd.concat([pd.read_csv(f, skiprows=7, index_col="rho_for_estimator") for f in all_filenames], axis=1)

combined_csv["final_sum_for_rho"] = combined_csv.sum(axis=1)

# Drop cols that are no longer needed
combined_csv.drop(["total_of_log_likelihood"], axis=1, inplace=True)

# reset index for sorting by values
combined_csv.reset_index(inplace=True)

# If multiple total_of_log_likelihoods are same, sort such that the ones connected to smallest rho_for_estimator
# are first
combined_csv.sort_values(by=["final_sum_for_rho", "rho_for_estimator"], ascending=[False, True], inplace=True)

max_rho = combined_csv.iloc[0][0]  # try max_rho, max_lk = df.loc[0, ['name_col1', 'name_col2']]
max_lk = combined_csv.iloc[0][1]

# Now that the required values, resort for final output
combined_csv.sort_values(by="rho_for_estimator", inplace=True)

with open(f"Aggregate_results_for_full_dataset.out", 'w') as file:  # open in write mode (create new file)
    file.write(f"Custom Metagenome Pairwise Recombination Rate Estimator\n")
    file.write(f"\n")

    file.write(f"Aggregate results for complete dataset\n")
    file.write(f"\n")

    file.write(f'-' * 55 + '\n')
    file.write(f"Theta estimate = {theta}\n")
    file.write(f"Max rho = {max_rho}, Max lk = {max_lk}\n")
    file.write(f'-' * 55 + '\n')

    file.write(f"\n")
    combined_csv.to_csv(file, mode='a', index=None)  # open in append mode (add to new file)
