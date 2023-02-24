#!/usr/bin/env python
import sys

prefix = sys.argv[1]

# prefix="rho_0.0_theta_0.0_sample_size_20_depth_1_genome_size_100000_seed_1_final_"

prefix_split = prefix.split("_")
samples = prefix_split[-10]

print(samples)
