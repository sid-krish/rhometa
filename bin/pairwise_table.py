#!/usr/bin/env python

import sys

import pandas as pd
import pysam


def loadAlignmentFile(bamfile):
    samfile = pysam.AlignmentFile(bamfile, "rb", threads=4)  # threads for decompression

    return samfile


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


def patternMatch(samFile, final_ref_pos_list, df):
    for read in samFile.fetch():
        # Without a contig or region all mapped reads in the file will be fetched.
        # The reads will be returned ordered by reference sequence, which will not necessarily be the order within the file.

        # Each read accessed moves down the reference sequence
        # If matches_only = True, only matched bases are returned - no None on either side.
        alnPairs = read.get_aligned_pairs(matches_only=True)
        # ^ a list of aligned read (query) and reference positions.
        # ^ Query position values will always be within the range of 0 to readLength - 1
        alnPairsReverseDict = {refPos: queryPos for queryPos, refPos in
                               alnPairs}  # reverse and convert to dictionary, key is refPos and values are queryPos
        readRefPosSet = alnPairsReverseDict.keys()  # Keys views are set-like since their entries are unique and hashable

        for pos1, pos2 in final_ref_pos_list:
            if pos1 in readRefPosSet and pos2 in readRefPosSet:
                queryPos1, queryPos2 = alnPairsReverseDict.get(pos1), alnPairsReverseDict.get(pos2)
                # baseAtPos1, baseAtPos2 = read.query_sequence[queryPos1], read.query_sequence[queryPos2]
                # This version doesn't have soft clipped bases, above does
                if queryPos1 <= (read.query_alignment_length - 1) and queryPos2 <= (read.query_alignment_length - 1):
                    # this is incorrect I mistakenly only considered the softclipped happening at 5' end
                    # maybe just keep soft clipped for now
                    baseAtPos1 = read.query_alignment_sequence[queryPos1]
                    baseAtPos2 = read.query_alignment_sequence[queryPos2]
                    if baseAtPos1 != 'N' and baseAtPos2 != 'N':
                        pair = str(baseAtPos1 + baseAtPos2)
                        df.at[(pos1, pos2), pair] += 1  # identify value by index and column
    return df


def export_final_df(final_df):
    final_df.index.name = "RefPos"
    final_df.to_csv("pairwise_table.csv")

    return None


if __name__ == "__main__":
    max_read_len = int(sys.argv[1])
    bam = sys.argv[2]
    vcf_file = sys.argv[3]

    # maxFragmentLen = 150
    # bam = "../Output/rho_30_sam_20_gen_30000/rho_30_sam_20_gen_30000_Aligned.csorted_fm_md.bam"
    # vcf_file = "../Output/rho_30_sam_20_gen_30000/rho_30_sam_20_gen_30000_lofreqOut.vcf"

    samFile = loadAlignmentFile(bam)
    genomeSize = int(samFile.header.get("SQ")[0].get("LN"))  # Header information can be accessed like a dictionary

    baseCombinations = ["AA", "AC", "AG", "AT", "CA", "CC",
                        "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]

    variant_positions = get_var_pos_from_vcf(vcf_file)

    final_ref_pos_list = get_final_ref_pos_list(variant_positions, max_read_len)

    init_df = get_init_df(final_ref_pos_list, baseCombinations)

    final_df = patternMatch(samFile, final_ref_pos_list, init_df)

    export_final_df(final_df)
