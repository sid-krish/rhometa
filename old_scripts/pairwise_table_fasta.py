#!/usr/bin/env python
import sys

import numpy as np
import pandas as pd
from Bio import SeqIO


def get_var_pos_from_variants_file(variants_file):
    df = pd.read_csv(variants_file)

    var_pos = df["RefPos_0-based"].to_list()

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
            pos1.append(start)
            pos2.append(i)

    return pos1, pos2


def get_init_df(var_pos_pairs, baseCombinations):
    init_df = pd.DataFrame(index=var_pos_pairs, columns=baseCombinations)
    init_df = init_df.fillna(value=0)  # all values will be NaN, convert to 0 instead

    return init_df


def pattern_match(input_fasta, var_pos_pairs, init_df):
    records = SeqIO.parse(input_fasta, "fasta")  # gets exhausted once accessed

    for record in records:
        genome_array = np.array(record.seq, dtype=np.str_)
        genome_pos_limit = genome_array.size - 1

        for var_pos_1, var_pos_2 in var_pos_pairs:
            # check if (size of current genome - 1) is less than largest pair in var_pos_pairs, to avoid
            # out of bound errors
            if var_pos_2 <= genome_pos_limit:
                base_at_var_pos_1 = genome_array[var_pos_1]
                base_at_var_pos_2 = genome_array[var_pos_2]

                pair = f"{base_at_var_pos_1}{base_at_var_pos_2}"
                init_df.at[(var_pos_1, var_pos_2), pair] += 1  # identify value by index and column

    return init_df


def export_final_df(final_df):
    final_df.index.name = "RefPos_0-based"
    final_df.to_csv("pairwise_table.csv")

    return None


if __name__ == "__main__":
    variants_file = sys.argv[1]
    input_fasta = sys.argv[2]

    # variants_file = "../Output/rho_30_sam_4_gen_20000/rho_30_sam_4_gen_20000_variants_in_fasta.csv"
    # input_fasta = "../Output/rho_30_sam_4_gen_20000/rho_30_sam_4_gen_20000_reformatted.fa"

    variant_positions = get_var_pos_from_variants_file(variants_file)

    pos1, pos2 = get_var_pos_pairs(variant_positions)

    var_pos_pairs = list(zip(pos1,pos2))

    baseCombinations = ["AA", "AC", "AG", "AT", "CA", "CC",
                        "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]

    init_df = get_init_df(var_pos_pairs, baseCombinations)

    final_df = pattern_match(input_fasta, var_pos_pairs, init_df)

    d_ij = np.array(pos2, dtype="int32") - np.array(pos1, dtype="int32")  # will be needed later

    final_df["d_ij"] = d_ij

    export_final_df(final_df)
