#!/usr/bin/env python
import pandas as pd
from ldpop import rhos_from_string
from utility import downsample

lookup_table = "lk_m_0.005_rho_0-1_s200.txt"
lookup_table_rho_range = "101,1"
# downsample_value = 100

for i in range(2,int(200 + 1)): # plus 1, so that the final val is processed
    downsample_value = i

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
