#!/usr/bin/env python
import sys

import m_pairwise_table_for_paired_end
import m_pairwise_table_single_end

if __name__ == '__main__':
    seq_type = int(sys.argv[1]) # 0 - single end, 1 - paired end
    bam_file = sys.argv[2]
    vcf_file = sys.argv[3]
    num_cores = int(sys.argv[4])
    fragment_len = int(sys.argv[5])

    # seq_type = 0  # 0 - single end, 1 - paired end
    # bam_file = "../Output/rho_10_sam_10_gen_10000/rho_10_sam_10_gen_10000_Aligned.bam"
    # vcf_file = "../Output/rho_10_sam_10_gen_10000/rho_10_sam_10_gen_10000_lofreqOut.vcf"
    # num_cores = 4

    if seq_type == 0:
        pairwise_table =  m_pairwise_table_single_end.main(bam_file, vcf_file, num_cores, fragment_len)

    elif seq_type == 1:
        pairwise_table = m_pairwise_table_for_paired_end.main(bam_file, vcf_file, num_cores, fragment_len)

    pairwise_table.to_pickle("pairwise_table.pkl")
