#!/usr/bin/env python
import multiprocessing
import os
import subprocess
import sys
from collections import defaultdict, namedtuple

import numpy as np
import pandas as pd
import pysam
import tqdm
from numba import njit
from numba.typed import List


def get_var_pos_from_vcf(vcf_file):
    """
    Extract variants from the VCF file.
    :param vcf_file: the vcf file to read
    :return: a per-reference dict of unique variant positions
    """
    f = pysam.VariantFile(vcf_file)

    var_pos = defaultdict(set)
    for snp in f.fetch():
        var_pos[snp.chrom].add(snp.pos)

    return var_pos


@njit(nogil=True)
def find_nearby(sites, window):
    """
    Create candidate variant pairs based on a window of separation. Performance
    is gained by writing a Numba compliant method.
    :param sites: the list of sites from a reference
    :param window: the window range measured in nucleotides
    :return: a list of unique pair tuples
    """
    ret = List()
    for i in range(len(sites)):
        for j in range(i + 1, len(sites)):
            dij = sites[j] - sites[i]
            if dij <= window:
                ret.append((sites[i] - 1, sites[j] - 1))
    return ret


def get_final_ref_pos_list(ref_sites, window):
    ref_pairs = {}
    for ref, sites in ref_sites.items():
        # re-wrap Numba array as a regular list as returned object is
        # 5x slower to access and cannot be pickled
        ref_pairs[ref] = list(
            find_nearby(np.sort(np.array(list(sites), dtype=np.int64)), window)
        )
    return ref_pairs


def exe_exists(exe_name: str) -> bool:
    """
    Check that a executable exists on the Path.
    :param exe_name: the base executable name
    :return: True, an executable file named exe_name exists and has executable bit set
    """
    p, f = os.path.split(exe_name)
    assert not p, "include only the base file name, no path specification"

    for pn in os.environ["PATH"].split(":"):
        full_path = os.path.join(pn, exe_name)
        if os.path.isfile(full_path) and os.access(full_path, os.X_OK):
            return True
    return False


# might need changes for single end
def count_bam_reads(
    file_name: str,
    paired: bool = False,
    mapped: bool = False,
    mapq: int = None,
    max_cpu: int = None,
) -> int:
    """
    Use samtools to quickly count the number of non-header lines in a bam file. This is assumed to equal
    the number of mapped reads.
    :param file_name: a bam file to scan (neither sorted nor an index is required)
    :param paired: when True, count only reads of mapped pairs (primary alignments only)
    :param mapped: when True, count only reads which are mapped
    :param mapq: count only reads with mapping quality greater than this value
    :param max_cpu: set the maximum number of CPUS to use in counting (otherwise all cores)
    :return: estimated number of mapped reads
    """
    assert exe_exists("samtools"), "required tool samtools was not found on path"
    assert exe_exists("wc"), "required tool wc was not found on path"
    assert not (paired and mapped), "Cannot set paired and mapped simultaneously"

    if not os.path.exists(file_name):
        raise IOError("{} does not exist".format(file_name))
    if not os.path.isfile(file_name):
        raise IOError("{} is not a file".format(file_name))

    opts = ["samtools", "view", "-c"]
    if max_cpu is None:
        max_cpu = multiprocessing.cpu_count()
    opts.append("-@{}".format(max_cpu))

    if paired:
        opts.append("-F0xC0C")
    elif mapped:
        opts.append("-F0xC00")

    if mapq is not None:
        assert 0 <= mapq <= 60, "mapq must be in the range [0,60]"
        opts.append("-q{}".format(mapq))

    proc = subprocess.Popen(
        opts + [file_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    value_txt = proc.stdout.readline().strip()
    try:
        count = int(value_txt)
        if paired and count % 2 != 0:
            print("When counting paired reads, the total was not divisible by 2.")
        return count
    except ValueError:
        raise RuntimeError(
            "Encountered a problem determining alignment count. samtools returned [{}]".format(
                value_txt
            )
        )


def remove_unaligned(in_pos, in_seq):
    """
    Conveinence method for removing elements from both position and sequence
    lists when the position is reported as "None". This indicates an unaligned
    position, which is of no use to this analysis.

    :param in_pos: input list of reference positions
    :param in_seq: input list of read-pair nucleotides
    :return: each list without Ns, if they occured.
    """
    out_pos = []
    out_seq = []
    for i in range(len(in_pos)):
        if in_pos[i] is None:
            continue
        out_pos.append(in_pos[i])
        out_seq.append(in_seq[i])
    return out_pos, out_seq


# Named tuple type for clarity.
Pair_t = namedtuple("Pair", ["chr", "pos1", "pos2", "base1", "base2"])


def pattern_match(bam, ref_pos_dict, read_count, n_proc):
    # for performance, we prepare a nested lookup structure using
    # dict and set so that variant existence can be quickly assessed.

    def _accept(_read):
        """
        Acceptance criteria for reads. At a minimum, reads must be mapped to be relevant. Supplied
        BAM files may not be filtered to exclude unmapped reads.
        """
        return not _read.is_unmapped

    # for every reference
    ref_lookup = {}
    for ref_name, pair_pos in ref_pos_dict.items():
        d = defaultdict(set)
        # pos1 indexes a map of int-to->set, while pos2 forms the set.
        for p1, p2 in pair_pos:
            d[p1].add(p2)
        ref_lookup[ref_name] = d

    with pysam.AlignmentFile(bam, "rb", threads=n_proc) as bam_file:

        # check BAM file ordering is query-name
        try:
            _header = bam_file.header["HD"]
        except KeyError as ex:
            raise RuntimeError("no header information in bam: {}".format(bam))
        if "SO" not in _header:
            raise RuntimeError("No sorting defined in bam: {}".format(bam))
        elif _header["SO"] != "queryname":
            raise RuntimeError("BAM file must be query-name sorted")

        # quick reference id to name lookup
        id_to_name = {
            i: bam_file.references[i] for i in range(len(bam_file.references))
        }

        pair_table = defaultdict(int)

        with tqdm.tqdm(total=read_count) as progress:

            bam_iter = bam_file.fetch(until_eof=True)
            while True:

                # Get read
                try:
                    r1 = next(bam_iter)
                    progress.update()

                    if r1.is_unmapped:
                        continue

                except StopIteration:
                    break

                # Get reference_id for read
                ref_id = r1.reference_id

                # prepare lookup map of position on reference to read nucleotide
                ref_positions = r1.get_reference_positions(full_length=True)
                query_sequence = r1.query_sequence

                ref_positions, query_sequence = remove_unaligned(
                    ref_positions, query_sequence
                )

                # quick lookup of reference position to read nucleotide
                pos_to_seq = dict(zip(ref_positions, query_sequence))

                ref_name = id_to_name[ref_id]
                # next read-pair if reference has no variants
                if ref_name not in ref_pos_dict:
                    continue

                # get the variant map for this reference
                snp_lookup = ref_lookup[ref_name]

                # check whether the mapped read-pair covers a variant pair
                # keep a tally each time we find one.
                # the nested loops will traverse only unique combinations of pairs
                pos = sorted(ref_positions)
                for i in range(len(pos)):
                    pos1 = pos[i]
                    if pos1 not in snp_lookup:
                        continue
                    for j in range(i + 1, len(pos)):
                        pos2 = pos[j]
                        if pos2 not in snp_lookup[pos1]:
                            continue
                        if pos_to_seq[pos1] != "N" and pos_to_seq[pos2] != "N":
                            base1 = pos_to_seq[pos1]
                            base2 = pos_to_seq[pos2]
                        else:
                            continue
                        # tally is indexed by 5 parameters
                        pair_table[Pair_t(ref_name, pos1, pos2, base1, base2)] += 1

        # convert the dict-based pairs table into a dataframe for downstream tools
        base_combinations = [
            "AA",
            "AC",
            "AG",
            "AT",
            "CA",
            "CC",
            "CG",
            "CT",
            "GA",
            "GC",
            "GG",
            "GT",
            "TA",
            "TC",
            "TG",
            "TT",
        ]
        d = {}
        for _pair, _count in pair_table.items():
            # rows are uniquely indexed by 3 parameters
            ix = (_pair.chr, _pair.pos1, _pair.pos2)
            if ix not in d:
                d[ix] = dict(zip(base_combinations, [0] * 16))
            d[ix][f"{_pair.base1}{_pair.base2}"] = _count
        # initialise the dataframe in one go
        pair_table = pd.DataFrame.from_dict(d, orient="index")

        return pair_table


def main(bam, vcf_file, num_cores, fragment_len):
    window_size = fragment_len

    # get number of aligned reads for progress tracking
    read_count = count_bam_reads(bam, max_cpu=num_cores)
    print("The BAM file contains {:,} reads".format(read_count))

    variant_positions = get_var_pos_from_vcf(vcf_file)
    print(
        "Number of variant positions to analyze: {:,}".format(
            sum(len(v) for v in variant_positions.values())
        )
    )

    reference_pair_positions = get_final_ref_pos_list(variant_positions, window_size)
    print(
        "Number of pair positions across references: {:,}".format(
            sum(len(v) for v in reference_pair_positions.values())
        )
    )

    pairwise_table = pattern_match(bam, reference_pair_positions, read_count, num_cores)

    pairwise_table.to_pickle("pairwise_table.pkl")


if __name__ == "__main__":
    bam_file = sys.argv[1]
    vcf_file = sys.argv[2]
    num_cores = int(sys.argv[3])
    fragment_len = int(sys.argv[4])

    main(bam_file, vcf_file, num_cores, fragment_len)
