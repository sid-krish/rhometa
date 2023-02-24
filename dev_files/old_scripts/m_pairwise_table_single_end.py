#!/usr/bin/env python

import pandas as pd
import pysam


def get_var_pos_from_vcf(vcf_file):
    f = pysam.VariantFile(vcf_file)

    var_pos = list({i.pos for i in f.fetch()})  # set comprehension to remove duplicates, then back to list
    var_pos.sort()

    return var_pos


def get_final_ref_pos_list(var_pos, max_read_len):
    final_ref_pos_list = []
    while var_pos:  # while list not empty
        start = var_pos.pop(0)
        for i in var_pos:  # once last item is popped this will stop
            difference = i - start
            # this version of the program is limited to looking at positions within reads as such if the distance
            # between position pairs are beyond read length it won't find anything
            if difference <= max_read_len:
                # final_ref_pos_list.append((start, i))
                # pysam uses 0-based index. Above line used earlier is wrong.
                final_ref_pos_list.append((start - 1, i - 1))

    return final_ref_pos_list


def get_init_df(final_ref_pos_list, baseCombinations):
    init_df = pd.DataFrame(index=final_ref_pos_list, columns=baseCombinations)
    init_df = init_df.fillna(value=0)  # all values will be NaN, convert to 0 instead

    return init_df


def pattern_match(samFile, final_ref_pos_list, df):
    for read in samFile.fetch(until_eof=True):
        # until_eof=True, also includes unmapped. Perhaps these should be filtered before hand?
        # needed if bam is not indexed

        read_ref_positions = read.get_reference_positions(full_length=True)
        # By default, this method only returns positions in the reference that are within the alignment.
        # If full_length is set, None values will be included for any soft-clipped or unaligned positions
        # within the read. The returned list will thus be of the same length as the read.

        # positions that are being looked at will be non soft-clipped
        # as those ref pos will be none as mentioned above and will be ignored
        read_query_seq = read.query_sequence

        ref_positions_set = set(read_ref_positions)  # for searching
        zip_ref_pos_and_query_seq = zip(read_ref_positions, read_query_seq)
        req_pos_and_query_seq_dict = {ref_pos: query_seq for ref_pos, query_seq in zip_ref_pos_and_query_seq}

        for pos1, pos2 in final_ref_pos_list:
            if pos1 in ref_positions_set and pos2 in ref_positions_set:
                base_at_pos1 = req_pos_and_query_seq_dict.get(pos1)
                base_at_pos2 = req_pos_and_query_seq_dict.get(pos2)

                if base_at_pos1 != 'N' and base_at_pos2 != 'N':
                    pair = str(base_at_pos1 + base_at_pos2)
                    df.at[(pos1, pos2), pair] += 1  # identify value by index and column
    return df


def main(bam, vcf_file, num_cores, fragment_len):
    max_read_len = fragment_len

    bam_file = pysam.AlignmentFile(bam, "rb", threads=num_cores)

    base_combinations = ["AA", "AC", "AG", "AT", "CA", "CC",
                         "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]

    variant_positions = get_var_pos_from_vcf(vcf_file)

    final_ref_pos_list = get_final_ref_pos_list(variant_positions, max_read_len)

    init_df = get_init_df(final_ref_pos_list, base_combinations)

    pairwise_table = pattern_match(bam_file, final_ref_pos_list, init_df)

    return pairwise_table
