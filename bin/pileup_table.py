#!/usr/bin/env python

import csv

import pysam

samfile = pysam.AlignmentFile("Aligned.sorted.sam")

with open('pileup_table.csv', 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile, delimiter=',')
    csv_writer.writerow(["RefPos_0-based", "Num_bases_aligned", "Alleles", "Allele_counts", "Bases"])

    for pileupcolumn in samfile.pileup():
        ref_pos = pileupcolumn.reference_pos
        num_aligned_bases = pileupcolumn.get_num_aligned()
        bases = pileupcolumn.get_query_sequences()
        bases_upper = [base.upper() for base in bases]
        alleles = set(bases_upper)
        allele_counts = len(alleles)

        csv_writer.writerow([ref_pos, num_aligned_bases, alleles, allele_counts, bases_upper])

    samfile.close()
