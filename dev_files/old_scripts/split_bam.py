import pysam


input_bam = "35793_R_combined_filtered.bam"

sam = pysam.AlignmentFile(input_bam, 'rb')

references = sam.references
print(f"{references=}")

ref_lens = sam.lengths
print(f"{ref_lens=}")

# print(sam.header.to_dict())

# print([(i,j) for i,j in (zip(references, ref_lens))])
for ref, ref_len in (zip(references, ref_lens)):
    sam_iter = sam.fetch(contig=ref, until_eof=True)

#     new_header = {'HD': {'VN': '1.0'},
#                   'SQ': [{'SN': ref, 'LN': ref_len}]}

    with pysam.AlignmentFile(f"{ref}.bam", "wb", reference_names=[ref], reference_lengths=[ref_len]) as out_bam:
        for entry in sam_iter:
            entry.reference_id = 0
            out_bam.write(entry)

    out_bam.close()
