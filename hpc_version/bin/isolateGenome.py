#!/usr/bin/env python
import sys

from Bio import SeqIO

inputFile = sys.argv[1]

records = list(SeqIO.parse(inputFile, "fasta"))
# print(">entry_" + records[0].id)  # first record
# print(records[0].seq)

with open('firstGenome.fa', 'w') as fileOut:
    fileOut.write('>' + records[0].id + '\n')
    fileOut.write(str(records[0].seq))
