#!/usr/bin/env python

import subprocess
import sys

import pandas as pd


def depth_distribution(pileup):
    pileup_fmt_cols = ["Sequence", "Position", "Reference_Base", "Depth", "Read_Results", "Quality"]
    pileup_df = pd.read_table(pileup, names=pileup_fmt_cols)  # pileup file is 1-based index

    pileup_df.drop(columns=["Sequence", "Reference_Base", "Read_Results", "Quality"], inplace=True)

    return pileup_df


if __name__ == '__main__':
    pileup = sys.argv[1]
    depth_range = sys.argv[2]
    bam_file = sys.argv[3]
    seed = int(sys.argv[4])

    # pileup = "../Theta_Est_Output/Aligned_sorted.pileup"
    # depth_range = "3,100"
    # bam_file = "../rho_10.0_theta_0.01_genome_size_100000.0_depth_5.0_seed_123.0_Aligned.bam"

    depth_lower_limit, depth_upper_limit = [int(i) for i in depth_range.split(',')]

    pileup_df = depth_distribution(pileup)
    pileup_depth_max = pileup_df["Depth"].max()

    depth_upper_limit_relative_proportion = round(depth_upper_limit / pileup_depth_max, 3)

    if depth_upper_limit_relative_proportion < 1.0:
        # the only time subsampling is required and possible
        # subsample seed fixed to 0 so that same results are produced each time
        frac = str(depth_upper_limit_relative_proportion).split('.')[1]
        subprocess.run(f"samtools view -b -s seed.{frac} {bam_file} > subsampled.bam", shell=True)

    else:
        # Just renaming and returning original file for nextflow
        subprocess.run(f"mv {bam_file} subsampled.bam", shell=True)
