#!/usr/bin/env python

from pathlib import Path
import os

bams = []
entries = Path('.')
for entry in entries.iterdir():
    # print(entry.name)
    if str(entry.name).endswith(".bam"):
        bams.append(entry.name)

samtools_bam_list = ' '.join(bams)
# print(samtools_bam_list)

os.system(f"samtools merge -n merged.bam {samtools_bam_list}")