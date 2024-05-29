import pandas as pd
import pysam


def get_allelel_freqs(vcf_file):
    # the variant allele frequencies might be able to help compute a substitution probability
    # use pysam to simply extract and return the full list of AF values

    pos_af = []
    f = pysam.VariantFile(vcf_file)

    for i in f:
        pos_af.append([i.pos,i.info["AF"][0]]) #i.info["AF"] outputs a tuple with the first entry containing the value. Also tested with other info tags

    return pos_af


def get_mean_sub_probability(df_pos_af,genome_len):
    #To get the overall genome-wide probability you'll need to average across all sites: \sum_{i \in sites} (2 * x_i * (1-x_i)) / |sites|. 

    sub_probabilities = []

    for i in df_pos_af["af"]:
        sp = 2 * i * (1-i)
        sub_probabilities.append(sp)

    print(sub_probabilities)

    mean_sub_probability = sum(sub_probabilities)/genome_len

    return mean_sub_probability


# def main():
# Objective: ρ (per site)/θ (per site) * tract length * substitution probability = r/m
# Compute substitution probability for recombination (nu) - not universal, then rest is available

genome_len = 10000
vcf_file = "139794_R_GCF_002101315_filtered_subsampled_freebayes_filt.vcf"

pos_af_list = get_allelel_freqs(vcf_file)

df_pos_af = pd.DataFrame(pos_af_list, columns=["pos","af"])

mean_sub_probability = get_mean_sub_probability(df_pos_af, genome_len)

# output mean_sub_probability along side theta estimates
print(mean_sub_probability)