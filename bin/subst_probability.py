import pandas as pd
import pysam


def get_AO_and_RO(vcf_file):
    # the variant allele frequencies might be able to help compute a substitution probability
    # use pysam to simply extract and return the full list of AF values

    pos_af = []
    f = pysam.VariantFile(vcf_file)

    for i in f:
        pos_af.append([i.pos,i.info["AO"][0],i.info["RO"]]) #i.info["AO"] outputs a tuple with the first entry containing the value. A few tags are like this

    return pos_af


def get_mean_sub_probability(df_pos_af,genome_len):
    #To get the overall genome-wide probability you'll need to average across all sites: \sum_{i \in sites} (2 * x_i * (1-x_i)) / |sites|. 

    sub_probabilities = []

    for i in df_pos_af["AF"]:
        sp = 2 * i * (1-i)
        sub_probabilities.append(sp)

    # print(sub_probabilities)

    mean_sub_probability = sum(sub_probabilities)/genome_len

    return mean_sub_probability


if __name__ == "__main__":
    # Objective: ρ (per site)/θ (per site) * tract length * substitution probability = r/m
    # Compute substitution probability for recombination (nu) - not universal, then rest is available

    genome_len = 1200090
    vcf_file = "139794_R_GCF_002101315_filtered_subsampled_freebayes_filt.vcf"

    pos_AO_RO_list = get_AO_and_RO(vcf_file)

    df_AO_RO = pd.DataFrame(pos_AO_RO_list, columns=["Pos","AO","RO"])

    df_AO_RO["AF"] = df_AO_RO.apply(lambda x : x["AO"]/(x["AO"]+x["RO"]), axis=1)

    mean_sub_probability = get_mean_sub_probability(df_AO_RO, genome_len)

    # # output mean_sub_probability along side theta estimates
    print(mean_sub_probability)