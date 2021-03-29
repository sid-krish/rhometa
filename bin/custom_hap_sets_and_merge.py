#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
from ldpop import rhos_from_string
from rho_splines import compute_splines


def d_ij(values):
    positions = [int(i) for i in values.strip("()").split(',')]
    pos1,pos2 = positions[0], positions[1]

    return pos2 - pos1


if __name__ == '__main__':
    lookup_table = sys.argv[1]
    pairwise_biallelic_table = sys.argv[2]
    lookup_format_csv = sys.argv[3]
    lookup_table_rho_range = sys.argv[4]

    df_lookup_format = pd.read_csv(lookup_format_csv, index_col="RefPos")

    # Convert to config
    # Optimised for ploidy = 1, refer to _get_configs in haplotype reader for full logic

    df = pd.DataFrame(index=range(df_lookup_format.shape[0]), columns=['0', '1', '2', '3', '4', '5', '6', '7'])
    
    df['0'] = df_lookup_format["00"].values
    df['1'] = df_lookup_format["01"].values
    df['3'] = df_lookup_format["10"].values
    df['4'] = df_lookup_format["11"].values
    df.fillna(value=0, inplace=True)

    configs = df.to_numpy(dtype='int64')  # has to be int 64
    # df2 = pd.DataFrame(data=configs)
    # df2.to_csv("lk_configs.csv")

    lookup_table_cols = ["Type", "#", "00", "01", "10", "11", "Rho"] + rhos_from_string(lookup_table_rho_range)
    df_lookup_table = pd.read_table(lookup_table, sep='\s+', skiprows=5, names=lookup_table_cols)

    df_lookup_table["00 01 10 11"] = df_lookup_table["00"].astype(str) + ' ' + \
                                     df_lookup_table["01"].astype(str) + ' ' + \
                                     df_lookup_table["10"].astype(str) + ' ' + \
                                     df_lookup_table["11"].astype(str)

    df_lookup_table.drop(columns=["Type", "#", "00", "01", "10", "11", "Rho"], inplace=True)

    df_lookup_table.set_index("00 01 10 11", inplace=True)

    likelihood_table_for_configs = compute_splines(configs, df_lookup_table)

    likelihoods_for_configs = likelihood_table_for_configs[0]

    df_likelihoods = pd.DataFrame(data=likelihoods_for_configs)

    df_likelihoods.columns = rhos_from_string(lookup_table_rho_range)

    df_pairwise_biallelic = pd.read_csv(pairwise_biallelic_table)

    # d_ij = df_pairwise_biallelic["d_ij"].tolist()

    df_likelihoods["d_ij"] = df_pairwise_biallelic["RefPos"].apply(lambda x: d_ij(x))

    df_likelihoods.to_csv("eq3.csv", index=None)
