#!/usr/bin/env python
import collections
from os import replace
import sys

import numpy as np
import pandas as pd


def resample(row, seed, numSamples):
    """ Perform up/down sampling by random sampling with replacement """

    # reproducible random state
    rng = np.random.default_rng(seed)

    # the 16 site pairs
    site_pairs = np.array(row.index, dtype="str")

    # site pair counts
    # site_pair_counts = np.array(row.values, dtype="int32")
    site_pair_counts = np.log1p(np.array(row.values, dtype="int32")) # smooth values with log

    # draw up to numSamples using the observed probabilities
    random_draw = rng.choice(site_pairs, numSamples, p=site_pair_counts / site_pair_counts.sum())
    # ISSUE 1: all pairwise entries have to be greater than numSamples or else replace = false won't work
    # it shouldn't work like this, as long as all pairwise values sum to numSamples it should work.
    # will need to find another solution

    # ISSUE 2: I want 1 to stay 1 not upscaled and it can't be done with replace

    # count them up
    final = collections.Counter(random_draw)

    # Counter objects have a dictionary interface except that they return a zero count for
    # missing items instead of raising a KeyError:
    return final.get("AA"), final.get("AC"), final.get("AG"), final.get("AT"),\
           final.get("CA"), final.get("CC"), final.get("CG"), final.get("CT"),\
           final.get("GA"), final.get("GC"), final.get("GG"), final.get("GT"),\
           final.get("TA"), final.get("TC"), final.get("TG"), final.get("TT")


if __name__ == '__main__':
    pairwise_table = sys.argv[1]
    seed = int(sys.argv[2])
    numSamples = int(sys.argv[3])

    # pairwise_table = "../Output/rho_15_sam_10_gen_20000/rho_15_sam_10_gen_20000_pairwise_table.csv"
    # seed = 123
    # numSamples = 10

    pt_df = pd.read_csv(pairwise_table, index_col="RefPos")

    site_pairs = pt_df.columns.values

    pt_df = pt_df.loc[~(pt_df == 0).all(axis=1)]  # drop rows where all rows vals are 0

    # pt_df = pt_df[pt_df.sum(axis=1) > numSamples] # select rows where row sum is >= sample size other wise its upsample

    resampled_df = pd.DataFrame(index=pt_df.index.values, columns=site_pairs)

    resampled_df[site_pairs] = pt_df.apply(lambda x: resample(x, seed, numSamples), axis=1, result_type='expand')

    resampled_df = resampled_df.replace(to_replace=np.nan, value=0)

    # resampled_df.to_csv("referenceDownSampled.csv")

    resampled_df.to_csv("pairwise_resampled.csv",index_label="RefPos")
