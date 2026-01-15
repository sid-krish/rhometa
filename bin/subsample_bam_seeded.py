#!/usr/bin/env python

import subprocess
import sys

# import pandas as pd
import numpy as np

def numpy_depth_max(read_depth):
    pileup_np = np.loadtxt(read_depth, dtype="int16", delimiter="\t", usecols=(2))

    return pileup_np.max()


if __name__ == "__main__":
    read_depth = sys.argv[1]
    depth_range = sys.argv[2]
    bam_file = sys.argv[3]
    seed = int(sys.argv[4])

    depth_lower_limit, depth_upper_limit = [int(i) for i in depth_range.split(",")]

    read_depth_max = numpy_depth_max(read_depth)

    depth_upper_limit_relative_proportion = round(
        depth_upper_limit / read_depth_max, 3
    )

    if depth_upper_limit_relative_proportion < 1.0:
        # the only time subsampling is required and possible
        print(depth_upper_limit, read_depth_max, depth_upper_limit_relative_proportion, end='')
        frac = str(depth_upper_limit_relative_proportion).split(".")[1]
        subprocess.run(
            f"samtools view -b -s {seed}.{frac} {bam_file} > subsampled.bam", shell=True
        )

    else:
        # Just renaming and returning original file for nextflow
        print(depth_upper_limit, read_depth_max, "1.0", end='')
        subprocess.run(f"mv {bam_file} subsampled.bam", shell=True)
        
