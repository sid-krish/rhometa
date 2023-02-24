#!/usr/bin/env python
import sys

import pandas as pd
from ldpop import rhos_from_string
from rho_splines import compute_splines


def d_ij(values):
    positions = [int(i) for i in values.strip("()").split(',')]
    pos1, pos2 = positions[0], positions[1]

    return pos2 - pos1


if __name__ == '__main__':
    # lookup_table = sys.argv[1]
    pairwise_biallelic_table = sys.argv[1]
    lookup_format_csv = sys.argv[2]
    lookup_table_rho_range = sys.argv[3]
    depth = int(sys.argv[4])
    lookup_table = sys.argv[5]

    df_lookup_format = pd.read_csv(lookup_format_csv, index_col="RefPos_0-based")

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

    df_lookup_table = pd.read_csv(lookup_table, index_col="00 01 10 11")

    likelihood_table_for_configs = compute_splines(configs, df_lookup_table)

    likelihoods_for_configs = likelihood_table_for_configs[0]

    df_likelihoods = pd.DataFrame(data=likelihoods_for_configs)

    df_likelihoods.columns = rhos_from_string(lookup_table_rho_range)

    df_pairwise_biallelic = pd.read_csv(pairwise_biallelic_table)

    # d_ij = df_pairwise_biallelic["d_ij"].tolist()

    df_likelihoods["d_ij"] = df_pairwise_biallelic["RefPos_0-based"].apply(lambda x: d_ij(x))

    df_likelihoods.to_csv(f"eq3_depth_{depth}.csv", index=None)
