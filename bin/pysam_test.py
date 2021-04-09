#!/usr/bin/env python

import pysam

bam = "../Output/rho_10_sam_10_gen_10000/rho_10_sam_10_gen_10000_Aligned.csorted_fm_md.bam"

bamfile = pysam.AlignmentFile(bam, "rb")

count = 1
for read in bamfile.fetch():
    alnPairs = read.get_aligned_pairs(matches_only=True)
    # if (read.query_alignment_sequence != read.query_sequence):
    print(count)
    print(read.query_name)
        # print(read.query_alignment_sequence)
        # print(read.query_sequence)
        # print(read.query_alignment_length)

    count += 1