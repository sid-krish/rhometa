#!/usr/bin/env python
import sys

from Bio import SeqIO

inputFile = sys.argv[1]

records = list(SeqIO.parse(inputFile, "fasta"))
# print(">entry_" + records[0].id)  # first record
# print(records[0].seq)

with open('reformatted.fa', 'w') as fileOut:
    for i in records:
        fileOut.write(">genome_" + i.id + '\n')
        fileOut.write(str(i.seq) + '\n')
