#!/usr/bin/env python
import pandas as pd
from rho_splines import compute_splines


def main(pairwise_biallelic_table, lookup_formatted_table, lookup_table_rho_vals, lookup_table):

    # Convert to config
    # Optimised for ploidy = 1, refer to _get_configs in haplotype reader for full logic

    df = pd.DataFrame(index=range(lookup_formatted_table.shape[0]), columns=['0', '1', '2', '3', '4', '5', '6', '7'])
    
    df['0'] = lookup_formatted_table["00"].values
    df['1'] = lookup_formatted_table["01"].values
    df['3'] = lookup_formatted_table["10"].values
    df['4'] = lookup_formatted_table["11"].values
    df.fillna(value=0, inplace=True)

    configs = df.to_numpy(dtype='int64')  # has to be int 64

    df_lookup_table = pd.read_csv(lookup_table, index_col="00 01 10 11")

    likelihood_table_for_configs = compute_splines(configs, df_lookup_table)

    # print(likelihood_table_for_configs)

    likelihoods_for_configs = likelihood_table_for_configs[0]

    table_ids_for_eq3 = likelihood_table_for_configs[2]

    df_likelihoods = pd.DataFrame(data=likelihoods_for_configs)

    df_likelihoods.columns = lookup_table_rho_vals

    df_likelihoods["d_ij"] = [pos_b - pos_a for ref, pos_a, pos_b in pairwise_biallelic_table.index.to_series()]

    return df_likelihoods, table_ids_for_eq3
