#!/usr/bin/env python
import sys

import numpy as np
import pandas as pd
import pysam


def get_var_pos_from_vcf(vcf_file):
    f = pysam.VariantFile(vcf_file)

    var_pos = list({i.pos for i in f.fetch()})  # set comprehension to remove duplicates, then back to list
    var_pos.sort()

    return var_pos


def get_var_pos_pairs(var_pos):
    """
    No need to adjust var pos pairs to be within read len since full genomes are being used
    There is a significant increase in data that needs to be processed as a result
    """

    pos1 = []
    pos2 = []
    while var_pos:  # while list not empty
        start = var_pos.pop(0)
        for i in var_pos:  # once last item is popped this will stop
            pos1.append(start-1)  # -1 because pysam uses 0-based indexing
            pos2.append(i-1)

    return pos1, pos2


def get_init_df(var_pos_pairs, baseCombinations):
    init_df = pd.DataFrame(index=var_pos_pairs, columns=baseCombinations)
    init_df = init_df.fillna(value=0)  # all values will be NaN, convert to 0 instead

    return init_df


def pattern_match(var_pos_pairs, df, init_df):

    # order of bases is based on pileup order
    for var_pos_0, var_pos_1 in var_pos_pairs:

        # Isolate relevant rows to be processed
        var_pos_0_row = df.loc[(var_pos_0)]
        var_pos_1_row = df.loc[(var_pos_1)]

        # Identify the smallest number of bases first, that will be the number pairs that can be found for the var pos set
        smallest_num_bases = min(var_pos_0_row.get("Num_bases_aligned"), var_pos_1_row.get("Num_bases_aligned"))
        var_pos_0_bases = var_pos_0_row.get("Bases").strip('[]').replace('\'','').replace(' ','').split(',')
        var_pos_1_bases = var_pos_1_row.get("Bases").strip('[]').replace('\'','').replace(' ','').split(',')

        # Get pair combination up to smallest_num_bases
        for i in range(smallest_num_bases):
            # print(f"pair: {var_pos_0_bases[i]}, {var_pos_1_bases[i]}")
            pair = f"{var_pos_0_bases[i]}{var_pos_1_bases[i]}"
            init_df.at[(var_pos_0, var_pos_1), pair] += 1

    return init_df


def export_final_df(final_df):
    final_df.index.name = "RefPos"
    final_df.to_csv("pairwise_table.csv")

    return None


if __name__ == "__main__":
    variants_file = sys.argv[1]
    vcf_file = sys.argv[2]

    # variants_file = "variants_from_pileup.csv"
    # vcf_file = "lofreqOut.vcf"

    df = pd.read_csv(variants_file, index_col="RefPos_0-based")

    variant_positions = get_var_pos_from_vcf(vcf_file) # better to do it this way to avoid false positives from manual checking pileups

    pos1, pos2 = get_var_pos_pairs(variant_positions)

    var_pos_pairs = list(zip(pos1,pos2))

    baseCombinations = ["AA", "AC", "AG", "AT", "CA", "CC",
                        "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]

    init_df = get_init_df(var_pos_pairs, baseCombinations)

    final_df = pattern_match(var_pos_pairs, df, init_df)

    # d_ij = np.array(pos2, dtype="int32") - np.array(pos1, dtype="int32")  # will be needed later

    # final_df["d_ij"] = d_ij

    export_final_df(final_df)
