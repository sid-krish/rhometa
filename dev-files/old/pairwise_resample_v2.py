#!/usr/bin/env python
import collections
import sys

import numpy as np
import pandas as pd
from sklearn.utils import resample


def get_resampled_values(row, seed, numSamples):
    ''' Downsampling algorithm'''

    row = row.dropna()
    row_idx = row.index.values
    row_vals = np.array(row.values, dtype='int32')  # perhaps log smoothing can done here to reduce impact of outliers?

    # In [151]: np.repeat(['a','b','c'], [1,2,3], axis=0)
    # Out[151]: array(['`a', 'b', 'b', 'c', 'c', 'c'], dtype='<U1')

    samples = np.repeat(row_idx, row_vals, axis=0)

    # Stratify to retain under represented samples, seems to work better in some cases than others
    # better to leave out for this particular issue
    resampled = resample(samples, replace=False, n_samples=numSamples, random_state=seed)

    # count them up
    final = collections.Counter(resampled)

    # Counter objects have a dictionary interface except that they return a zero count for
    # missing items instead of raising a KeyError:
    return final.get("AA"), final.get("AC"), final.get("AG"), final.get("AT"), \
           final.get("CA"), final.get("CC"), final.get("CG"), final.get("CT"), \
           final.get("GA"), final.get("GC"), final.get("GG"), final.get("GT"), \
           final.get("TA"), final.get("TC"), final.get("TG"), final.get("TT")


if __name__ == '__main__':
    pairwise_table = sys.argv[1]
    seed = int(sys.argv[2])
    numSamples = int(sys.argv[3])

    # pairwise_table = "../Output/rho_45_sam_30_gen_40000/rho_45_sam_30_gen_40000_pairwise_table.csv"
    # seed = 123
    # numSamples = 30

    pt_df = pd.read_csv(pairwise_table, index_col="RefPos_0-based")

    site_pairs = pt_df.columns.values

    # pt_df = pt_df.loc[~(pt_df == 0).all(axis=1)]  # drop rows where all rows vals are 0

    # select rows where row sum is >= numSamples to down sample to other wise its up sampling
    pt_df = pt_df[pt_df.sum(axis=1) >= numSamples]

    # Remove outlier rows, this and the above line can be used to essentialy slice out entries that fall within a certain range
    pt_df = pt_df[pt_df.sum(axis=1) <= 150]

    pt_df = pt_df.replace(to_replace=0, value=np.nan)

    resampled_df = pd.DataFrame(index=pt_df.index.values, columns=site_pairs)

    resampled_df[site_pairs] = pt_df.apply(lambda x: get_resampled_values(x, seed, numSamples), axis=1,
                                           result_type='expand')

    resampled_df = resampled_df.replace(to_replace=np.nan, value=0)

    # resampled_df.to_csv("referenceDownSampled.csv")

    resampled_df.to_csv("pairwise_resampled.csv", index_label="RefPos_0-based")
