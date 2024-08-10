#!/usr/bin/env python
import sys
import numpy as np
import pysam
import subprocess


def get_depth_vals(vcf_file):
    f = pysam.VariantFile(vcf_file)

    depth_vals = [i.info["DP"] for i in f]

    return np.array(depth_vals)


if __name__ == "__main__":
    cpus = int(sys.argv[1])
    vcf_file = sys.argv[2]
    snp_qual = int(sys.argv[3])
    min_snp_depth = int(sys.argv[4])
    cutoff_percentage = 100 - int(sys.argv[5])  # top n percent to cut off

    min_RO_AO = 2 

    '''
    Rationale for top n percent depth cutoff
    anomalous depth is a symptom of some kind of phenomenon that violates the recombination model assumptions. 
    It could be that a gene has undergone one or more segmental duplication events, or it could be that it's 
    a slow-evolving gene and reads from other similar species are therefore able to map because the sequence identity is high enough (ie "conserved regions"). 
    '''
    depth_vals = get_depth_vals(vcf_file)
    percentile_val = np.percentile(depth_vals[depth_vals >= min_snp_depth], cutoff_percentage)
    # Since we only keep entries with values more than the min_snp_depth cut off, the top n% cutoff is evaluated based on  remaning values where values >= min_snp_depth

    command = f"""bcftools filter --threads {cpus} -i 'TYPE="snp" && QUAL>={snp_qual} && FORMAT/DP>={min_snp_depth} && FORMAT/DP<={round(percentile_val)} && FORMAT/RO>={min_RO_AO} && FORMAT/AO>={min_RO_AO}' {vcf_file} > freebayes_filt.vcf"""

    process = subprocess.run(command, shell=True, check=True)
