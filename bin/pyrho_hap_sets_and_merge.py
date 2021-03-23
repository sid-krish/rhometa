#!/usr/bin/env python
import sys

import pandas as pd

# pyrho functions
from haplotype_reader import parse_vcf_to_genos, genos_to_configs
from ldpop import rhos_from_string
from rho_splines import compute_splines


if __name__ == '__main__':
    # lookup_table = sys.argv[1]
    # pairwise_biallelic_table = sys.argv[2]
    # vcf_file = sys.argv[3]
    # lookup_table_rho_range = sys.argv[4]
    # genome_size = int(sys.argv[5])

    lookup_table = "../Output/rho_15_sam_10_gen_20000/rho_15_sam_10_gen_20000_lookupTable.txt"
    # pairwise_biallelic_table = "pairwise_biallelic_table.csv"
    vcf_file = "../Output/rho_15_sam_10_gen_20000/calls.vcf"
    lookup_table_rho_range = "101,100"
    genome_size = 20000

    vcf_to_genos = parse_vcf_to_genos(vcf_file, 1)
    # genos = seqs_to_genos[0]

    # genos_to_configs = genos_to_configs(genos, genome_size, 1)
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
    #
    # df_lookup_table.drop(columns=["Type", "#", "00", "01", "10", "11", "Rho"], inplace=True)
    #
    # df_lookup_table.set_index("00 01 10 11", inplace=True)
    #
    # likelihood_table_for_configs = compute_splines(configs, df_lookup_table)
    #
    # likelihoods_for_configs = likelihood_table_for_configs[0]
    #
    # df_likelihoods = pd.DataFrame(data=likelihoods_for_configs)
    #
    # df_likelihoods.columns = rhos_from_string(lookup_table_rho_range)
    #
    # df_pairwise_biallelic = pd.read_csv(pairwise_biallelic_table)
    #
    # d_ij = df_pairwise_biallelic["d_ij"].tolist()
    #
    # df_likelihoods["d_ij"] = d_ij
    #
    # df_likelihoods.to_csv("eq3.csv", index=None)

    # in utility.py
    # have a look at
    # if n00 < n11:
    #     n00, n11 = n11, n00
    # if n01 < n10:
    #     n01, n10 = n10, n01
    # if n11 + n01 > sample_size//2:
    #     n00, n01, n10, n11 = n01, n00, n11, n10