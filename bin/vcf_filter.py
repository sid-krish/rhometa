#!/usr/bin/env python
import sys
import pysam
import numpy as np


def get_depth_vals(vcf_file):

    depth_vals = [i.info["DP"] for i in vcf_file]

    return np.array(depth_vals)


if __name__ == "__main__":
    cpus = int(sys.argv[1])
    vcf_file = sys.argv[2]
    snp_qual = int(sys.argv[3])
    min_snp_depth = int(sys.argv[4])
    cutoff_percentage = 100 - int(sys.argv[5])  # top n percent to cut off

    # min_RO_AO = 2 not checking RO only AO filter applied at freebayes command

    '''
    Rationale for top n percent depth cutoff
    anomalous depth is a symptom of some kind of phenomenon that violates the recombination model assumptions. 
    It could be that a gene has undergone one or more segmental duplication events, or it could be that it's 
    a slow-evolving gene and reads from other similar species are therefore able to map because the sequence identity is high enough (ie "conserved regions"). 
    '''

    vcf_in = pysam.VariantFile(vcf_file)

    bcf_out = pysam.VariantFile('freebayes_filt.vcf', 'w', header=vcf_in.header)

    depth_vals = get_depth_vals(vcf_in)
    percentile_val = np.percentile(depth_vals[depth_vals >= min_snp_depth], cutoff_percentage)

    for i in vcf_in.fetch():
        if i.info["TYPE"] in [('snp',), ('snp','snp'), ('snp','snp','snp')]:
            if i.qual >= snp_qual and i.info["DP"] >= min_snp_depth and i.info["DP"] <= round(percentile_val):
                # print(i)
                bcf_out.write(i)

    bcf_out.close()
