#!/usr/bin/env python
import sys

bam = sys.argv[1]
prepend_file = sys.argv[2]
seed = sys.argv[3]

# bam = "a.s_bs_ds.bam"
# prepend_file = "none"

if prepend_file == "none":
    file_name = bam.rsplit(".",1)[0]
    print(f"{seed}_{file_name}_", end = '')

else:
    print(f"{seed}_{prepend_file}", end = '')
