#!/usr/bin/env python

import sys

import pandas as pd
import pysam


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
            # difference = i - start
            # this version of the program is limited to looking at positions within reads as such if the distance
            # between position pairs are beyond read length it won't find anything
            # if difference <= windows_size:
                # final_ref_pos_list.append((start, i))
                # pysam uses 0-based index. Above line used earlier is wrong.
            final_ref_pos_list.append((start - 1, i - 1))

    return final_ref_pos_list


def get_init_df(final_ref_pos_list, baseCombinations):
    init_df = pd.DataFrame(index=final_ref_pos_list, columns=baseCombinations)
    init_df = init_df.fillna(value=0)  # all values will be NaN, convert to 0 instead

    return init_df


def pattern_match(bam_File, final_ref_pos_list, df):
    read_1 = read_2 = None  # initialise before loop

    # until_eof=True, also includes unmapped. Perhaps these should be filtered before hand?
    # needed if bam is not indexed
    for read in bam_File.fetch(until_eof=True):
        # Given that the reads are next to each other the pair are combined and processed accordingly

        if read.is_read1:
            read_1 = read
            read_1_ref_positions = read_1.positions
            read_1_query_seq = read_1.query_sequence  # using query seq for now, includes soft clipped

        elif read.is_read2:
            read_2 = read
            read_2_ref_positions = read_2.positions
            read_2_query_seq = read_2.query_sequence  # using query seq for now, includes soft clipped

        if read_1 and read_2:
            # Combine and process

            # combined read 1 and read 2 ref positions, the 2 ref positions we are searching should be in here
            combined_ref_positions = read_1_ref_positions + read_2_ref_positions
            combined_ref_positions_set = set(combined_ref_positions)  # for searching
            combined_query_sequences = read_1_query_seq + read_2_query_seq

            zip_ref_pos_and_query_seq = list(zip(combined_ref_positions, combined_query_sequences))

            req_pos_and_query_seq_dict = {ref_pos: query_seq for ref_pos, query_seq in zip_ref_pos_and_query_seq}

            for pos1, pos2 in final_ref_pos_list:
                if pos1 in combined_ref_positions_set and pos2 in combined_ref_positions_set:
                    base_at_pos1 = req_pos_and_query_seq_dict.get(pos1)
                    base_at_pos2 = req_pos_and_query_seq_dict.get(pos2)

                    if base_at_pos1 != 'N' and base_at_pos2 != 'N':
                        pair = f"{base_at_pos1}{base_at_pos2}"
                        df.at[(pos1, pos2), pair] += 1  # get value by index and column

            # Reset for next pair of reads
            read_1 = read_2 = None

    bam_File.close()

    return df


def export_final_df(final_df):
    final_df.index.name = "RefPos"
    final_df.to_csv("pairwise_table.csv")

    return None


if __name__ == "__main__":
    window_size = int(sys.argv[1]) # read_len * 2 + mean_frag_len + std_dv * 2.
    bam = sys.argv[2]
    vcf_file = sys.argv[3]

    window_size = 900  # read_len * 2 + mean_frag_len + std_dv * 2. Hard cored for testing
    # bam = "Aligned.bam"  # untouched bwa output
    # vcf_file = "lofreqOut.vcf"

    bam_file = pysam.AlignmentFile(bam, "rb", threads=4)

    # genomeSize = int(bamfile.header.get("SQ")[0].get("LN"))  # Header information can be accessed like a dictionary

    base_combinations = ["AA", "AC", "AG", "AT", "CA", "CC",
                         "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]

    variant_positions = get_var_pos_from_vcf(vcf_file)

    final_ref_pos_list = get_final_ref_pos_list(variant_positions, window_size)

    init_df = get_init_df(final_ref_pos_list, base_combinations)

    final_df = pattern_match(bam_file, final_ref_pos_list, init_df)

    export_final_df(final_df)
