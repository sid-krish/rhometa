import pysam

def load_vcf():


    return vcf


def get_allelel_freqs(vcf_file):
    #the variant allele frequencies might be able to help compute a substitution probability

    f = pysam.VariantFile(vcf_file)

    # use pysam to simply extract and return the full list of AF values
    #### Or maybe I can even extract out the column with pandas, Aaron? 

    return af_list


def get_mean_sub_probability(af_list,genome_len):
    #To get the overall genome-wide probability you'll need to average across all sites: \sum_{i \in sites} (2 * x_i * (1-x_i)) / |sites|. 

    sub_probabilities = []

    for i in af_list:
        sp = 2 * i * (1-i)
        sub_probabilities.append(sp)

    mean_sub_probability = sub_probabilities/genome_len

    return sub_probability


def main():
    # Objective: ρ (per site)/θ (per site) * tract length * substitution probability = r/m
    # Compute substitution probability, then rest is available

    genome_len = "some_val"
    vcf_file_loc = "some_location"

    vcf = load_vcf(vcf_file_loc)

    af_list = get_allelel_freqs(vcf)

    mean_sub_probability = get_mean_sub_probability(af_list,genome_len)

    # output mean_sub_probability along side theta estimates