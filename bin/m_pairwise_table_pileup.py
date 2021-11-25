#!/usr/bin/env python

import pandas as pd
import sys
import pysam

## Just for testing do not use. The ordering of bases is different to the normal method

def get_var_pos_from_vcf(vcf_file):
    f = pysam.VariantFile(vcf_file)

    var_pos = list({i.pos for i in f.fetch()})  # set comprehension to remove duplicates, then back to list
    var_pos.sort()

    return var_pos


def get_final_ref_pos_list(var_pos, windows_size):
    final_ref_pos_list = []
    while var_pos:  # while list not empty
        start = var_pos.pop(0)
        for i in var_pos:  # once last item is popped this will stop
            difference = i - start
            # this version of the program is limited to looking at positions within reads as such if the distance
            # between position pairs are beyond read length it won't find anything
            if difference <= windows_size:
                # final_ref_pos_list.append((start, i))
                # pysam uses 0-based index. Above line used earlier is wrong.
                final_ref_pos_list.append((start - 1, i - 1))

    return final_ref_pos_list

# def get_final_ref_pos_list(var_pos, windows_size):
#     final_ref_pos_list = []
#     while var_pos:  # while list not empty
#         start = var_pos.pop(0)
#         for i in var_pos:  # once last item is popped this will stop
#             final_ref_pos_list.append((start - 1, i - 1))

#     return final_ref_pos_list


def pileup_lists(bamfile):
    ref_pos = []
    coverage = []
    pileups = []

    for pileupcolumn in bamfile.pileup():
        # print("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
        # print(pileupcolumn.get_query_sequences())
        ref_pos.append(pileupcolumn.reference_pos)
        coverage.append(pileupcolumn.n)
        pileups.append(pileupcolumn.get_query_sequences())

    bamfile.close()

    return ref_pos,coverage,pileups


def get_init_df(final_ref_pos_list, baseCombinations):
    init_df = pd.DataFrame(index=final_ref_pos_list, columns=baseCombinations)
    init_df = init_df.fillna(value=0)  # all values will be NaN, convert to 0 instead

    return init_df


def process_variant_pairs(df_pileups,final_ref_pos_list, init_df):

    for pos1, pos2 in final_ref_pos_list:
        # get relevant rows from df
        # print(f"{pos1} {df_pileups.loc[pos1].to_list()}")
        pos1_list = df_pileups.loc[pos1].to_list()
        pos2_list = df_pileups.loc[pos2].to_list()

        # process elements in arrays and append to dataframe until a 'N' is hit

        for pos1_array_val, pos2_list_val in zip(pos1_list,pos2_list):

            pair = f"{pos1_array_val}{pos2_list_val}"
            if 'N' in pair:
                continue
            
            else:
                init_df.at[(pos1, pos2), pair] += 1  # get value by index and column

    return init_df


if __name__ == '__main__':
    bamfile_name = sys.argv[1]
    vcf_file = sys.argv[2]
    num_cores = int(sys.argv[3])
    window_size = int(sys.argv[4])

    # vcf_file = "../Sim_Gen_Output/lofreqOut.vcf"
    # bamfile_name = "../Sim_Gen_Output/Aligned_sorted.bam"
    # num_cores = 4
    # window_size = 300

    bamfile = pysam.AlignmentFile(bamfile_name, "rb", threads=num_cores)
    ref_pos, coverage, pileups = pileup_lists(bamfile)

    df_pileups = pd.DataFrame(pileups)
    df_pileups = df_pileups.fillna('N') # replace none with 'N
    df_pileups = df_pileups.applymap(lambda s: s.upper()) # make everything upper case
    df_pileups["ref_pos"] = ref_pos
    df_pileups.set_index("ref_pos", inplace=True)


    variant_positions = get_var_pos_from_vcf(vcf_file)

    final_ref_pos_list = get_final_ref_pos_list(variant_positions, window_size)

    base_combinations = ["AA", "AC", "AG", "AT", "CA", "CC",
                         "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]

    init_df = get_init_df(final_ref_pos_list, base_combinations)

    pairwise_table = process_variant_pairs(df_pileups,final_ref_pos_list, init_df)

    pairwise_table.to_pickle("pairwise_table.pkl")