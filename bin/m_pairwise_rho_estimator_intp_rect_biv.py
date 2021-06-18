#!/usr/bin/env python
import csv
import sys

import numpy as np
import pandas as pd
from scipy.interpolate import RectBivariateSpline

from ldpop import rhos_from_string


def eq2(x, y, z, x_new, y_new, merged_eq3_table):
    # going with kx=1, ky=1 produces expected results, same as linear interpolator
    f = RectBivariateSpline(x, y, z, kx=1, ky=1)

    intpltd_res = f.ev(xi=x_new, yi=y_new)

    intpltd_res_reshape = intpltd_res.reshape(y.shape[0], merged_eq3_table.shape[0]).T

    return intpltd_res_reshape


def export_to_csv(total_log_likelihoods, depth):
    with open(f"collected_likelihoods_depth_{depth}.csv", 'w') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["rho_for_estimator", "total_of_log_likelihood"])

        for idx, val in total_log_likelihoods.iteritems():
            csv_writer.writerow([idx, val])

    return None


def main(merged_eq3_table, table_ids_for_eq3, p_ij_grid, lookup_table, depth):

    merged_eq3_table.drop(columns=["d_ij"], inplace=True)  # was used for p_ij grid and no longer needed

    df_lookup_table = pd.read_csv(lookup_table)
    df_lookup_table.drop(columns=["00 01 10 11"], inplace=True)
    # from here the index of the df_lookup_table (default 0-based row identifier) will be used as lookup config id

    x = df_lookup_table.index.values.astype("int32")  # row index 0,1,2,3...
    y = df_lookup_table.columns.values.astype("float64")  # the rows 0.0, 1.0, 2.0...
    z = df_lookup_table.to_numpy().astype("float64")  # values in cells
    
    # a = np.array([0, 1, 2]) # np.tile(a, 2) # Out[10]: array([0, 1, 2, 0, 1, 2])
    x_new = np.tile(table_ids_for_eq3, y.shape[0])
    y_transposed = p_ij_grid.to_numpy().astype("float64").T  # must be transposed
    y_new = y_transposed.ravel()  # same as flatten() but a bit faster
    
    # Equation 2
    intpltd_res = eq2(x, y, z, x_new, y_new, merged_eq3_table)
    intpltd_res_df = pd.DataFrame(data=intpltd_res, columns=y)

    return intpltd_res_df
    
    # Equation 1
    # total_log_likelihoods = intpltd_res_df.sum(axis=0)
    
    # Output to csv
    # export_to_csv(total_log_likelihoods, depth)
