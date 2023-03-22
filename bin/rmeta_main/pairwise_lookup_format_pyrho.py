#!/usr/bin/env python

import numpy as np


def allele_distribution(row, allele_freqs):
    # drop site pairs that don't have values then get remaining pairs as list
    pairs = row.dropna().index.to_list()

    # get alleles in site 1 then get unique
    sites_1_alleles = set([list(i)[0] for i in pairs])

    sites_2_alleles = set([list(i)[1] for i in pairs])

    # filter allele_freqs for required alleles
    s1_global_allele_counts = {i: allele_freqs[i] for i in sites_1_alleles}

    s2_global_allele_counts = {i: allele_freqs[i] for i in sites_2_alleles}

    # https://www.kite.com/python/answers/how-to-find-the-max-value-in-a-dictionary-in-python
    s1_max_allele = max(s1_global_allele_counts, key=s1_global_allele_counts.get)
    s1_min_allele = min(s1_global_allele_counts, key=s1_global_allele_counts.get)

    s2_max_allele = max(s2_global_allele_counts, key=s2_global_allele_counts.get)
    s2_min_allele = min(s2_global_allele_counts, key=s2_global_allele_counts.get)

    # 0 is max allele, 1 is min allele
    alleles_as_00 = f"{s1_max_allele}{s2_max_allele}"
    alleles_as_01 = f"{s1_max_allele}{s2_min_allele}"
    alleles_as_10 = f"{s1_min_allele}{s2_max_allele}"
    alleles_as_11 = f"{s1_min_allele}{s2_min_allele}"

    config_00 = row[alleles_as_00]
    config_01 = row[alleles_as_01]
    config_10 = row[alleles_as_10]
    config_11 = row[alleles_as_11]

    # going with the order "00,01,10,11" for this is how ldpop orders the lookup table
    lookup_config = [config_00, config_01, config_10, config_11]

    return lookup_config


def main(pairwise_biallelic_table):
    df = pairwise_biallelic_table.replace(to_replace=0, value=np.nan)

    allele_freqs = {"A": 2, "C": 3, "G": 4, "T": 5}  # pyrho allele digitisation

    df[["00", "01", "10", "11"]] = df.apply(
        lambda x: allele_distribution(x, allele_freqs), axis=1, result_type="expand"
    )

    df = df.replace(to_replace=np.nan, value=0)

    cols_to_drop = [
        "AA",
        "AC",
        "AG",
        "AT",
        "CA",
        "CC",
        "CG",
        "CT",
        "GA",
        "GC",
        "GG",
        "GT",
        "TA",
        "TC",
        "TG",
        "TT",
    ]

    df.drop(columns=cols_to_drop, inplace=True)

    return df
