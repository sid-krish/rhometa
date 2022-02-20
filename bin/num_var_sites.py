#!/usr/bin/env python

import sys
import pysam

def get_var_pos_from_vcf(vcf_file):
    f = pysam.VariantFile(vcf_file)

    var_pos = {i.pos for i in f.fetch()}  # set comprehension to remove duplicates

    f.close()

    return var_pos


if __name__ == '__main__':
    vcf = sys.argv[1]

    num_variant_positions = len(get_var_pos_from_vcf(vcf))
    sys.stdout.write(str(num_variant_positions)) # for next steps