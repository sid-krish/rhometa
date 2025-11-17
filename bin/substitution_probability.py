#!/usr/bin/env python

import sys
import pysam


def get_variant_allele_frequency(vcf_file):
    """
    The variant allele frequencies can help compute a substitution probability
    Compute variant allele frequency for each position in the VCF file
    Logic and formula used for calculating AF https://help.galaxyproject.org/t/calculating-variant-allele-frequency-from-freebayes-vcf/1630
    AF = AO / (AO + RO), assuming  VAF as the frequency of the most observed variant allele:
    """
    pos_af = []
    f = pysam.VariantFile(vcf_file)

    for i in f:
        # split AO on comma
        # find the maximum in the resulting list of values and use it as the numerator for your ratio
        numerator = max(i.info["AO"])

        # sum up all values in the list and add the RO value => use the result as the denominator
        denominator = sum(i.info["AO"]) + i.info["RO"]

        # divide the numerator by the denominator to get the AF
        pos_af.append(numerator/denominator)

    return pos_af


def get_mean_sub_probability(AF_list, genome_len):
    '''
     If the variant is at abundance x then the substitution probability can be thought of 
     as the probability of choosing a pair of cells for recombination where one has the 
     variant and the other doesn't. The probability for that site is 2 * x * (1-x).
     For a site with no variants anywhere in the population, you can see that the 
     probability of introducing mutation via recombination is zero, and when a variant 
     is at abundance 0.5 the probability is 0.5. To get the overall 
     genome-wide probability you'll need to average 
     across all sites: sum(2 * x_i * (1-x_i) for i in sites)) / |sites|.
    '''

    sub_probabilities = []

    for i in AF_list:
        sp = 2 * i * (1-i)
        sub_probabilities.append(sp)

    mean_sub_probability = sum(sub_probabilities)/genome_len

    return mean_sub_probability


if __name__ == "__main__":
    # Objective: ρ (per site)/θ (per site) * tract length * substitution probability = r/m
    # Compute substitution probability for recombination (nu) in this script

    genome_len = int(sys.argv[1])
    vcf_file = sys.argv[2]

    # The bams are subsampled and from the rho estimation pipeline, 
    # we have to use the variants actually considered for the rho estimate.

    AF_list = get_variant_allele_frequency(vcf_file)

    mean_sub_probability = get_mean_sub_probability(AF_list, genome_len)

    print(mean_sub_probability)