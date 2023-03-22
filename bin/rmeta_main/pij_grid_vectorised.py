#!/usr/bin/env python
import numpy as np
import pandas as pd


def main(recom_tract_len, lookup_table_rho_vals, merged_eq3_table):
    rho_vals_to_test = np.array(lookup_table_rho_vals, dtype="float64")

    # should be largest rho in lookup table, should not go beyond this
    max_rho_val_to_test = rho_vals_to_test[-1]

    d_ij = merged_eq3_table["d_ij"].to_numpy(dtype="int32")

    # a = np.array([0, 1, 2]) # np.tile(a, 2) # Out[10]: array([0, 1, 2, 0, 1, 2])
    d_ij_vals_to_process = np.tile(d_ij, rho_vals_to_test.shape[0])
    rho_vals_to_process = np.tile(rho_vals_to_test, d_ij.shape[0])
    rho_vals_to_process.sort()

    # LDhat implentation in pairdip.c lk_calc https://github.com/auton1/LDhat/blob/master/pairdip.c
    p_ij = (
        2
        * rho_vals_to_process
        * (1 - np.exp(-1 * d_ij_vals_to_process / recom_tract_len))
    )

    # if val over max rho, select max rho (which would be minimum)
    p_ij = np.minimum(p_ij, max_rho_val_to_test)

    p_ij_reshape = p_ij.reshape(rho_vals_to_test.shape[0], d_ij.shape[0]).T

    p_ij_df = pd.DataFrame(data=p_ij_reshape, columns=rho_vals_to_test)

    return p_ij_df
