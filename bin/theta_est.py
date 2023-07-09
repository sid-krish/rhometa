#!/usr/bin/env python

import sys
from math import log

import numba as nb
import pandas as pd
import pysam


def get_var_pos_from_vcf(vcf_file):
    f = pysam.VariantFile(vcf_file)

    # set comprehension to remove duplicates
    var_pos = {i.pos for i in f.fetch()}

    f.close()

    return var_pos


def depth_distribution(pileup):
    pileup_fmt_cols = [
        "Sequence",
        "Position",
        "Reference_Base",
        "Depth",
        "Read_Results",
        "Quality",
    ]

    # pileup file is 1-based index
    pileup_df = pd.read_table(pileup, names=pileup_fmt_cols)

    pileup_df.drop(
        columns=["Sequence", "Reference_Base", "Read_Results", "Quality"], inplace=True
    )

    # Watterson estimate only applicable for depth > 1
    pileup_df = pileup_df[pileup_df["Depth"] > 1]

    return pileup_df


@nb.jit(nopython=True)
def watterson_estimate(segregating_sites, genome_len, samples):
    """
    Theta per site
    """
    # k and a_n are traditional notations https://en.wikipedia.org/wiki/Watterson_estimator
    # calculate k (per site)
    k = log(genome_len / (genome_len - segregating_sites))

    # calculate a_n
    a_n = 0
    for i in range(1, int(samples)):
        a_n += 1 / i

    # ldhat watterson per site
    theta = k / a_n

    return theta


if __name__ == "__main__":
    genome_size = int(sys.argv[1])
    pileup = sys.argv[2]
    vcf = sys.argv[3]

    # genome_size = 50000
    # pileup = "/Users/sid/Documents/GitHub/rhometa/Theta_Est_Output/theta_estimate/recom_0.005_tract_1000_mutation_0.005_sample_size_20_depth_20_genome_size_50000_seed_123_final_Aligned_sorted.pileup"
    # vcf = "/Users/sid/Documents/GitHub/rhometa/Theta_Est_Output/freebayes/recom_0.005_tract_1000_mutation_0.005_sample_size_20_depth_20_genome_size_50000_seed_123_final_freebayes_filt.vcf"

    num_variant_positions = len(get_var_pos_from_vcf(vcf))
    pileup_df = depth_distribution(pileup)

    # Summary stats for depth
    depth_summary_stats_df = pileup_df["Depth"].describe(percentiles=[0.5])
    depth_summary_stats_df = pd.DataFrame(depth_summary_stats_df).T
    # makes sense to round depth vals to the nearest int
    depth_summary_stats_df = depth_summary_stats_df.round(0)

    depth_summary_stats_df = depth_summary_stats_df.drop(
        columns=["count", "std"]
    )  # Drop cols that can't be used for theta

    # Summary stats for theta based on depth summary stats
    depth_summary_stats_df = depth_summary_stats_df.drop(columns=["min", "max"])
    depth_summary_stats_df.rename(columns={"mean": "pileup_mean_depth", "50%": "pileup_median_depth"}, inplace=True)

    depth_summary_stats_df["pileup_mean_depth"] = [
        int(i) for i in depth_summary_stats_df["pileup_mean_depth"]
    ]

    depth_summary_stats_df["pileup_median_depth"] = [
        int(i) for i in depth_summary_stats_df["pileup_median_depth"]
    ]

    # Summary stats for theta based on depth summary stats
    theta_sum_stats_df = pd.DataFrame()
    theta_sum_stats_df = depth_summary_stats_df.applymap(
        lambda x: watterson_estimate(num_variant_positions, genome_size, x)
    ).copy()
    theta_sum_stats_df.rename(columns={"pileup_mean_depth": "theta_per_site_mean_depth",
                                       "pileup_median_depth": "theta_per_site_median_depth"}, inplace=True)

    theta_sum_stats_df["theta_per_site_mean_depth"] = [
        float("%.3g" % i) for i in theta_sum_stats_df["theta_per_site_mean_depth"]
    ]

    theta_sum_stats_df["theta_per_site_median_depth"] = [
        float("%.3g" % i) for i in theta_sum_stats_df["theta_per_site_median_depth"]
    ]

    # Export all summary stats
    all_summary_stats_df = pd.concat([depth_summary_stats_df, theta_sum_stats_df], axis=1)

    all_summary_stats_df = all_summary_stats_df.drop(columns=["pileup_mean_depth", "theta_per_site_mean_depth"])
    all_summary_stats_df["genome_length"] = genome_size
    all_summary_stats_df.to_csv("Theta_estimate_stats.csv", index=False)
