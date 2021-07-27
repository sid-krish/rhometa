#!/usr/bin/env python
import sys

from pysam import Fastafile

inputFile = sys.argv[1]

# inputFile = "rho_15_sam_15_gen_10000_seqgenOut.fa"

fasta = Fastafile(inputFile)

ref_and_len = zip(fasta.references, fasta.lengths)

with open('reformatted.fa', 'w') as fileOut:
    for identifier, length in ref_and_len:
        sequence = fasta.fetch(reference=identifier, start=0, end=length)

        fileOut.write(f">genome_{identifier}\n")
        fileOut.write(f"{sequence}\n")

fasta.close()
