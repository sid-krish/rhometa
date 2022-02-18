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
    Theta per site
    """

    # calculate k
    k = log(genome_len / (genome_len - segregating_sites))

    # calculate a_n
    a_n = 0
    for i in range(1, samples):
        a_n += 1/i

    # original theta formula
    theta = k/a_n

    return theta


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

    # genome_size = 100000
    # pileup = "../Theta_Est_Output/rho_0.001_theta_0.01_sample_size_25_depth_4_genome_size_100000_seed_0_final_Aligned_sorted.pileup"
    # vcf = "../Theta_Est_Output/rho_0.001_theta_0.01_sample_size_25_depth_4_genome_size_100000_seed_0_final_freeBayesOut.vcf"

    num_variant_positions = len(get_var_pos_from_vcf(vcf))
    pileup_df = depth_distribution(pileup)
    plot_depth(pileup_df)

    # Summary stats for depth
    depth_summary_stats_df = pileup_df["Depth"].describe(percentiles = [.5])
    depth_summary_stats_df = pd.DataFrame(depth_summary_stats_df).T
    depth_summary_stats_df = depth_summary_stats_df.round(0)  # makes sense to round depth vals to the nearest int
    depth_summary_stats_df = depth_summary_stats_df.drop(columns=["count", "std"])  # Drop cols that can't be used for theta

    # Plot theta for min - max depth range
    min_depth, max_depth = int(depth_summary_stats_df["min"][0]), int(depth_summary_stats_df["max"][0])
    plot_theta(num_variant_positions, genome_size, min_depth, max_depth)

    # Summary stats for theta based on depth summary stats
    depth_summary_stats_df = depth_summary_stats_df.drop(columns=["min", "max"])
    theta_sum_stats_df = pd.DataFrame()
    theta_sum_stats_df = depth_summary_stats_df.applymap(lambda x: watterson_estimate(num_variant_positions,
                                                                                      genome_size,
                                                                                      x))

    # Export all summary stats
    all_summary_stats_df = depth_summary_stats_df.append(theta_sum_stats_df)
    all_summary_stats_df.rename(columns={'mean': 'mean_depth', '50%': 'median_depth'}, inplace=True)
    all_summary_stats_df = pd.melt(all_summary_stats_df, value_vars=["mean_depth", "median_depth"])
    new_idx = ["mean_depth", "theta_per_site_at_mean_depth", "median_depth", "theta_per_site_at_median_depth"]
    all_summary_stats_df.index = new_idx
    all_summary_stats_df = all_summary_stats_df.drop(columns=["variable"])
    all_summary_stats_df.to_csv("Theta_estimate_stats.csv", header=False)