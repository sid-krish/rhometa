#!/usr/bin/env python
import pandas as pd
import sys

def compute_final_angsd_theta(input_file):

    # Read the input CSV file
    df = pd.read_table(input_file)

    tW_sum = df['tW'].sum()
    nSites_sum = df['nSites'].sum() #The final column is the effetive number of sites with data in the window. -> so not all references sites only the useable ones.

    # Calculate the final thetaW value
    final_theta = tW_sum / nSites_sum

    return final_theta

if __name__ == "__main__":
    angsd_tsv = sys.argv[1]
    final_theta_output = compute_final_angsd_theta(angsd_tsv)
    print(final_theta_output)