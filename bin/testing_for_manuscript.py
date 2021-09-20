#!/usr/bin/env python

import pandas as pd
from ldpop import rhos_from_string

import m_isolate_by_depth
import m_biallelic_filter_pairwise_table
import m_pairwise_lookup_format_pyrho
import m_custom_hap_sets_and_merge
import m_pij_grid_vectorised
import m_pairwise_rho_estimator_intp_rect_biv

if __name__ == '__main__':

    recom_tract_len = 500
    depth_range = "3,250"
    n_resamples = 20
    lookup_table_rho_range = "101,100"
    pairwise_table_file = "Recom_Est_Output/pairwise_table.pkl"
    num_cores = 4

    # load pairwise table
    pairwise_table = pd.read_pickle(pairwise_table_file)

    # isolate a single depth for testing
    pairwise_table_slice = m_isolate_by_depth.main(pairwise_table, 5)

    # perform bi-allelic filtering
    pairwise_biallelic_table = m_biallelic_filter_pairwise_table.main(pairwise_table_slice.copy())