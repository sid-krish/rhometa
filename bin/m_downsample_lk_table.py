#!/usr/bin/env python

import sys
import pandas as pd

from ldpop import rhos_from_string
from utility import downsample


# argument list to make things simpler for multiple processing algorithm
def downsample_lookup_table(lookup_table, lookup_table_rho_range, downsample_value, output_fmt='csv'):

    lookup_table_cols = ["Type", "#", "00", "01", "10",
                         "11", "Rho"] + rhos_from_string(lookup_table_rho_range)
                         
    df_lookup_table = pd.read_table(
        lookup_table, sep='\s+', skiprows=5, names=lookup_table_cols)

    df_lookup_table["00 01 10 11"] = df_lookup_table["00"].astype(str) + ' ' + \
        df_lookup_table["01"].astype(str) + ' ' + \
        df_lookup_table["10"].astype(str) + ' ' + \
        df_lookup_table["11"].astype(str)

    df_lookup_table.drop(
        columns=["Type", "#", "00", "01", "10", "11", "Rho"], inplace=True)
    df_lookup_table.index.set_names("00 01 10 11", inplace=True)

    downsampled_lookuptable_df = downsample(df_lookup_table, downsample_value)

    lk_downsampled_out_dir = f"lk_downsampled_{downsample_value}.csv"

    if output_fmt == 'pickle':
        downsampled_lookuptable_df.to_pickle(lk_downsampled_out_dir)
    elif output_fmt == 'csv':
        downsampled_lookuptable_df.to_csv(lk_downsampled_out_dir)
    else:
        raise RuntimeError('unimplemented format: {}'.formt(output_fmt))

    return None


if __name__ == '__main__':
    lookup_table = sys.argv[1]
    lookup_table_rho_range = sys.argv[2]
    downsample_value = int(sys.argv[3])
    output_fmt = sys.argv[4]

    downsample_lookup_table(lookup_table, lookup_table_rho_range, downsample_value, output_fmt)
