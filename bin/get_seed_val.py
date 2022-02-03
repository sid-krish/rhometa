#!/usr/bin/env python
import sys

prefix = sys.argv[1]

prefix_split = prefix.split("_")
seed = prefix_split[-3]

print(seed)
