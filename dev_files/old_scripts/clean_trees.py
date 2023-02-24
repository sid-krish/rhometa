#!/usr/bin/env python
import sys


def clean_trees(inputFile, outputFile):
    with open(inputFile, 'r') as fileIn, open(outputFile, 'w') as fileOut:
        for line in fileIn:
            if line.startswith('['):
                # This gets rid of SITE and other unnecessary data
                fileOut.write(line)


if __name__ == '__main__':
    inputFile = sys.argv[1]
    # inputFile = "../Output(ref)/Output_mflen_300_r_15_s_15_g_10000/rho_15_sam_15_gen_10000_trees.txt"

    outputFile = "cleanTrees.txt"

    clean_trees(inputFile, outputFile)
