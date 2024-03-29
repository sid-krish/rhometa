#!/usr/bin/env python

import subprocess
import sys

# import pandas as pd
import numpy as np


# def depth_distribution(pileup):
#     pileup_fmt_cols = ["Sequence", "Position", "Reference_Base", "Depth", "Read_Results", "Quality"]
#     pileup_df = pd.read_table(pileup, usecols=[3], dtype='int16')  # pileup file is 1-based index
#
#     print(pileup_df.memory_usage(deep=True))
#
#     # pileup_df.drop(columns=["Sequence", "Reference_Base", "Read_Results", "Quality"], inplace=True)
#
#     return pileup_df


def numpy_depth_max(pileup):
    pileup_np = np.loadtxt(pileup, dtype="int16", delimiter="\t", usecols=(3))

    return pileup_np.max()


if __name__ == "__main__":
    pileup = sys.argv[1]
    depth_range = sys.argv[2]
    bam_file = sys.argv[3]
    seed = int(sys.argv[4])

    # pileup = "123_paired_end_Aligned_sorted.pileup"
    # depth_range = "3,50"
    # bam_file = "123_paired_end_Aligned_sorted.bam"

    depth_lower_limit, depth_upper_limit = [int(i) for i in depth_range.split(",")]

    pileup_depth_max = numpy_depth_max(pileup)

    depth_upper_limit_relative_proportion = round(
        depth_upper_limit / pileup_depth_max, 3
    )

    if depth_upper_limit_relative_proportion < 1.0:
        # the only time subsampling is required and possible
        print(depth_upper_limit, pileup_depth_max, depth_upper_limit_relative_proportion, end='')
        frac = str(depth_upper_limit_relative_proportion).split(".")[1]
        subprocess.run(
            f"samtools view -b -s {seed}.{frac} {bam_file} > subsampled.bam", shell=True
        )

    else:
        # Just renaming and returning original file for nextflow
        print(depth_upper_limit, pileup_depth_max, "1.0", end='')
        subprocess.run(f"mv {bam_file} subsampled.bam", shell=True)
        
