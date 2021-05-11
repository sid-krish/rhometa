#!/usr/bin/env python
import sys

import numpy as np
import pandas as pd

from ldpop import rhos_from_string


if __name__ == '__main__':
    genome_len = int(sys.argv[1])
    recom_tract_len = int(sys.argv[2])
    ldpop_rho_range = sys.argv[3]
    merged_pairwise_and_lookup_table = sys.argv[4]
    depth = int(sys.argv[5])

    # genome_len = int(20000)
    # recom_tract_len = int(500)
    # ldpop_rho_range = "101,100"
    # merged_pairwise_and_lookup_table = "../Output/rho_30_sam_4_gen_20000/rho_30_sam_4_gen_20000_eq3.csv"

    rho_vals_to_test = np.array(rhos_from_string(ldpop_rho_range), dtype="float64")

    max_rho_val_to_test = rho_vals_to_test[-1]  # should be largest rho in lookup table, should not go beyond this

    eq3_df = pd.read_csv(merged_pairwise_and_lookup_table, usecols=["d_ij"])
    d_ij = eq3_df["d_ij"].to_numpy(dtype="int32")

    # a = np.array([0, 1, 2]) # np.tile(a, 2) # Out[10]: array([0, 1, 2, 0, 1, 2])
    d_ij_vals_to_process = np.tile(d_ij, rho_vals_to_test.shape[0])
    rho_vals_to_process = np.tile(rho_vals_to_test, d_ij.shape[0])
    rho_vals_to_process.sort()

    p_ij = 2 * rho_vals_to_process * (1 - np.exp(-1 * d_ij_vals_to_process / recom_tract_len))
    p_ij = np.minimum(p_ij, max_rho_val_to_test)  # if val over max rho, select max rho (which would be minimum)

    p_ij_reshape = p_ij.reshape(rho_vals_to_test.shape[0], d_ij.shape[0]).T

    p_ij_df = pd.DataFrame(data=p_ij_reshape, columns=rho_vals_to_test)
    p_ij_df.to_csv(f"p_ij_grid_depth_{depth}.csv", index=False)
