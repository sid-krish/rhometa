#!/usr/bin/env python

import sys
import pandas as pd
from multiprocessing import Pool
from ldpop import rhos_from_string
from utility import downsample


def downsample_lookup_table(arg_list):  # argument list to make things simpler for multiple processing algorithm

    lookup_table = arg_list[0]
    lookup_table_rho_range =  arg_list[1]
    downsample_value = arg_list[2]

    lookup_table_cols = ["Type", "#", "00", "01", "10", "11", "Rho"] + rhos_from_string(lookup_table_rho_range)
    df_lookup_table = pd.read_table(lookup_table, sep='\s+', skiprows=5, names=lookup_table_cols)

    df_lookup_table["00 01 10 11"] = df_lookup_table["00"].astype(str) + ' ' + \
                                     df_lookup_table["01"].astype(str) + ' ' + \
                                     df_lookup_table["10"].astype(str) + ' ' + \
                                     df_lookup_table["11"].astype(str)

    df_lookup_table.drop(columns=["Type", "#", "00", "01", "10", "11", "Rho"], inplace=True)

    df_lookup_table.set_index("00 01 10 11", inplace=True)

    downsampled_lookuptable_df = downsample(df_lookup_table, downsample_value)

    lk_downsampled_out_dir = f"lk_downsampled_{downsample_value}.csv"

    downsampled_lookuptable_df.to_csv(lk_downsampled_out_dir, index_label="00 01 10 11")

    return None


if __name__ == '__main__':
    lookup_table = sys.argv[1]
    lookup_table_rho_range = sys.argv[2]
    lk_table_max_depth = int(sys.argv[3])
    num_cores = int(sys.argv[4])

    # lookup_table = "../Lookup_tables/lookup_table.txt"
    # lookup_table_rho_range = "101,100"
    # lk_table_max_depth = 50
    # num_cores = 4

    downsample_range = [i for i in range(3, lk_table_max_depth + 1)]

    downsample_lookup_table_arg_list = []
    for i in downsample_range:
        downsample_lookup_table_arg_list.append([lookup_table, lookup_table_rho_range, i])

    # Parallel method execution
    with Pool(processes = num_cores) as p:
        p.map(downsample_lookup_table, downsample_lookup_table_arg_list)





