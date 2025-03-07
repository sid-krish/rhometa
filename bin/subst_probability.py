#!/usr/bin/env python

import sys
import pandas as pd
import pysam


def get_AO_and_RO(vcf_file):
    # the variant allele frequencies can help compute a substitution probability
    # use pysam to simply extract and return the full list of AF values

    pos_af = []
    f = pysam.VariantFile(vcf_file)

    for i in f:
        print(i)
        pos_af.append([i.pos,i.info["AO"][0],i.info["RO"]])
        # i.info["AO"] outputs a tuple. For the rhometa pipeline there is only one value in the tuple
        # (in the case of multi allele sites there may be more) 
        # @Aaron I'm unsure why we dont see more is it because I filter by SNPs only?

    return pos_af


def get_mean_sub_probability(df_pos_af,genome_len):
    '''
     If the variant is at abundance x then the substitution probability can be thought of 
     as the probability of choosing a pair of cells for recombination where one has the 
     variant and the other doesn't. The probability for that site is 2 * x * (1-x).
     For a site with no variants anywhere in the population, you can see that the 
     probability of introducing mutation via recombination is zero, and when a variant 
     is at abundance 0.5 the probability is 0.5. To get the overall 
     genome-wide probability you'll need to average 
     across all sites: sum(2 * x_i * (1-x_i) for i in sites)) / |sites|.
     @Aaron I reworte this bit in a way that I understand better
    '''

    sub_probabilities = []

    for i in df_pos_af["AF"]:
        sp = 2 * i * (1-i)
        sub_probabilities.append(sp)

    # print(sub_probabilities)

    mean_sub_probability = sum(sub_probabilities)/genome_len

    return mean_sub_probability


if __name__ == "__main__":
    # Objective: ρ (per site)/θ (per site) * tract length * substitution probability = r/m
    # Compute substitution probability for recombination (nu) in this script

    # genome_len = int(sys.argv[1])
    # vcf_file = sys.argv[2]

    genome_len = 50000
    vcf_file = "/mnt/c/Users/sid/Documents/GitHub/rhometa/Rho_Est_Output/freebayes/123_paired_end_freebayes_filt.vcf"
    # The bams are subsampled but it's from the rho estimation pipeline, 
    # we have to use the variants actually considered for the rho estimate.

    pos_AO_RO_list = get_AO_and_RO(vcf_file)

    df_AO_RO = pd.DataFrame(pos_AO_RO_list, columns=["Pos","AO","RO"])
    # @Aaron I only seem to be observing one AO value, possibly due to fitlers applied with freebayes?
    # Logic for calculating Allele Frequency
    # Forumla for calculating AF https://help.galaxyproject.org/t/calculating-variant-allele-frequency-from-freebayes-vcf/1630

    df_AO_RO["AF"] = df_AO_RO.apply(lambda x : x["AO"]/(x["AO"]+x["RO"]), axis=1) 

    mean_sub_probability = get_mean_sub_probability(df_AO_RO, genome_len)

    # Output mean_sub_probability along side theta estimates
    print(mean_sub_probability)