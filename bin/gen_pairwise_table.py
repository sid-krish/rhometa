#!/usr/bin/env python
import sys

import m_pairwise_table_for_paired_end

if __name__ == '__main__':
    seq_type = int(sys.argv[1]) # 0 - single end, 1 - paired end
    bam_file = sys.argv[2]
    vcf_file = sys.argv[3]

    # seq_type = 1  # 0 - single end, 1 - paired end
    # bam_file = "rho_20_sam_20_gen_10000_Aligned.bam"
    # vcf_file = "rho_20_sam_20_gen_10000_lofreqOut.vcf"

    if seq_type == 0:
        pairwise_table = None  # TODO

    elif seq_type == 1:
        pairwise_table = m_pairwise_table_for_paired_end.main(bam_file, vcf_file)

    pairwise_table.to_pickle("pairwise_table.pkl")
