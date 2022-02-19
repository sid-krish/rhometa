#!/usr/bin/env python
import sys
from math import log

genome_len = int(sys.argv[1])
segregating_sites = int(sys.argv[2])
num_samples = int(sys.argv[3])

# calculate k
k = log(genome_len / (genome_len - segregating_sites))
# k = segregating_sites

# calculate a_n
a_n = 0
for i in range(1, num_samples):
    a_n += 1/i

# original theta formula
theta = k/a_n
print(theta)