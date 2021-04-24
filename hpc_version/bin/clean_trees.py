#!/usr/bin/env python
import sys

inputFile = sys.argv[1]

with open(inputFile, 'r') as fileIn, open("cleanTrees.txt", 'w') as fileOut:
    for line in fileIn:
        if line.startswith('['):
            # This gets rid of SITE and other unnecessary data
            fileOut.write(line)
