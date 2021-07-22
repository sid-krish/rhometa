#!/usr/bin/env python

import sys
from math import log

import numpy as np
import pandas as pd
import pysam

import seaborn as sns
sns.set_theme(style="darkgrid")


def get_var_pos_from_vcf(vcf_file):
    f = pysam.VariantFile(vcf_file)

    var_pos = {i.pos for i in f.fetch()}  # set comprehension to remove duplicates

    f.close()

    return var_pos


def depth_distribution(pileup):
    pileup_fmt_cols = ["Sequence", "Position", "Reference_Base", "Depth", "Read_Results", "Quality"]
    pileup_df = pd.read_table(pileup, names=pileup_fmt_cols)  # pileup file is 1-based index

    pileup_df.drop(columns=["Sequence", "Reference_Base", "Read_Results", "Quality"], inplace=True)

    # Watterson estimate only applicable for depth > 1
    pileup_df = pileup_df[pileup_df["Depth"] > 1]

    return pileup_df


def watterson_estimate(segregating_sites, genome_len, samples):
    """
    Theta is at population level
    Might have to do theta/2 at some point for haplotypes? -> keeping original formula to maintain parity with ldhat
    """

    k = 1
    calc_sum = 0
    while k < samples:
        calc_sum += 1 / k
        k += 1

    sum_inverse = 1 / calc_sum

    calc_log = log(genome_len / (genome_len - segregating_sites))

    theta = sum_inverse * calc_log

    return np.round(theta, decimals=5)  # rounding to 5 decimals, reasonable accuracy


def plot_depth(pileup_df):

    g = sns.histplot(x="Depth", data=pileup_df)
    g.set(title='Counts of depths/read counts across genome')
    g.figure.savefig("depth_distribution.png", dpi=500)
    g.figure.clf()  # clear figure for next plot

    return None


def plot_theta(num_variant_positions, genome_size, min_depth, max_depth):

    theta_vals_for_depth_range = []
    for i in range(min_depth, max_depth + 1):
        theta = watterson_estimate(num_variant_positions,genome_size,i)
        theta_vals_for_depth_range.append(theta)

    depths = [i for i in range(min_depth, max_depth + 1)]
    f = sns.lineplot(x=depths, y=theta_vals_for_depth_range)
    f.set(xlabel='Depth', ylabel='Theta', title='Theta estimated for min-max depth range')
    f.figure.savefig("theta_estimates.png", dpi=500)
    f.figure.clf()  # clear figure for next plot

    return None


if __name__ == '__main__':

    genome_size = int(sys.argv[1])
    pileup = sys.argv[2]
    vcf = sys.argv[3]

    # genome_size = 10000
    # pileup = "../Output/rho_10_sam_10_gen_10000/rho_10_sam_10_gen_10000_Aligned_sorted.pileup"
    # vcf = "../Output/rho_10_sam_10_gen_10000/rho_10_sam_10_gen_10000_lofreqOut.vcf"

    num_variant_positions = len(get_var_pos_from_vcf(vcf))
    pileup_df = depth_distribution(pileup)
    plot_depth(pileup_df)

    # Summary stats for depth
    depth_summary_stats_df = pileup_df["Depth"].describe(percentiles = [.05,.25, .5, .75, .95])
    depth_summary_stats_df = pd.DataFrame(depth_summary_stats_df).T
    depth_summary_stats_df = depth_summary_stats_df.round(0)  # makes sense to round depth vals to nearest int
    depth_summary_stats_df = depth_summary_stats_df.drop(columns=["count", "std"]) # Drop cols that can't be used for theta

    # Plot theta for min - max depth range
    min_depth, max_depth = int(depth_summary_stats_df["min"][0]), int(depth_summary_stats_df["max"][0])
    plot_theta(num_variant_positions, genome_size, min_depth, max_depth)

    # Summary stats for theta based on depth summary stats
    theta_sum_stats_df = pd.DataFrame()
    theta_sum_stats_df = depth_summary_stats_df.applymap(lambda x: watterson_estimate(num_variant_positions,
                                                                                      genome_size,
                                                                                      x))

    # Export all summary stats
    all_summary_stats_df = depth_summary_stats_df.append(theta_sum_stats_df)
    new_idx = ["Depth", "Theta_Est"]
    all_summary_stats_df.index = new_idx
    all_summary_stats_df.to_csv("Theta_estimate_stats.csv")