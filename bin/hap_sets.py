#!/usr/bin/env python

import allel
from haplotype_reader import parse_vcf_to_genos, genos_to_configs
from rho_splines import compute_splines
from ldpop import rhos_from_string
import pandas as pd

if __name__ == '__main__':
    lookup_table = "../Output/rho_15_sam_10_gen_20000/rho_15_sam_10_gen_20000_lookupTable.txt"
    # pairwise_biallelic_table = "pairwise_biallelic_table.csv"
    vcf_file = "../Output/rho_15_sam_10_gen_20000/calls.vcf"
    lookup_table_rho_range = "101,100"
    genome_size = 20000
    ploidy = 1

    callset = allel.read_vcf(vcf_file, numbers={'GT': ploidy})  # Filter ploidy 1
    cols = sorted(callset.keys())
    ht = allel.HaplotypeArray(callset['calldata/GT'])

    # genos_to_configs = genos_to_configs(ht, genome_size, ploidy)
    # configs = genos_to_configs[0]
    #
    # lookup_table_cols = ["Type", "#", "00", "01", "10", "11", "Rho"] + rhos_from_string(lookup_table_rho_range)
    # df_lookup_table = pd.read_table(lookup_table, sep='\s+', skiprows=5, names=lookup_table_cols)
    #
    # df_lookup_table["00 01 10 11"] = df_lookup_table["00"].astype(str) + ' ' + \
    #                                  df_lookup_table["01"].astype(str) + ' ' + \
    #                                  df_lookup_table["10"].astype(str) + ' ' + \
    #                                  df_lookup_table["11"].astype(str)
    #
    # df_lookup_table.drop(columns=["Type", "#", "00", "01", "10", "11", "Rho"], inplace=True)
    #
    # df_lookup_table.set_index("00 01 10 11", inplace=True)
    #
    # likelihood_table_for_configs = compute_splines(configs, df_lookup_table)

    # likelihoods_for_configs = likelihood_table_for_configs[0]