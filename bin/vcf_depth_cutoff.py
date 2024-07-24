#!/usr/bin/env python
import sys
import numpy as np
import pysam
import seaborn as sns

def get_depth_vals(vcf_file):
    depth_vals = []
    f = pysam.VariantFile(vcf_file)

    depth_vals = [i.info["DP"] for i in f]

    return np.array(depth_vals)


if __name__ == "__main__":
    # vcf_file = sys.argv[1]
    # cutoff_percentage = 100 - int(sys.argv[2]) # top n percent to cut off

    vcf_file = "/mnt/c/Users/sid/Documents/GitHub/rhometa/test/1_ERR315858_GCA_902540595_filtered_subsampled_freebayes_raw.vcf"
    cutoff_percentage = 100 - int(5)

    depth_vals = get_depth_vals(vcf_file)
    percentile_val = np.percentile(depth_vals, cutoff_percentage)

    print(round(percentile_val))
    filt_depth_vals = depth_vals[depth_vals <= round(percentile_val)]

    g = sns.histplot(x=filt_depth_vals)
    g.set(title="vcf_depth_vals")
    g.figure.savefig("filt_depth_distribution.png", dpi=500)
    g.figure.clf()  # clear figure for next plot

