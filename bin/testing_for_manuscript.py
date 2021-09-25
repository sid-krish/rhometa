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
    lookup_table_rho_vals = rhos_from_string(lookup_table_rho_range)
    lookup_table = f"Lookup_tables/lk_downsampled_{5}.csv"

    # load pairwise table
    pairwise_table = pd.read_pickle(pairwise_table_file)

    # isolate a single depth for testing
    pairwise_table_slice = m_isolate_by_depth.main(pairwise_table, 5)

    # perform bi-allelic filtering
    pairwise_biallelic_table = m_biallelic_filter_pairwise_table.main(pairwise_table_slice.copy())

    # convert to lookup format to match against likelihood tables
    lookup_formatted_table = m_pairwise_lookup_format_pyrho.main(pairwise_biallelic_table.copy())

    # Merge lookup formatted table on likelihood table
    merged_eq3_table, table_ids_for_eq3 = m_custom_hap_sets_and_merge.main(pairwise_biallelic_table.copy(),
                                                                           lookup_formatted_table.copy(),
                                                                           lookup_table_rho_vals,
                                                                           lookup_table)

    # Calculate p_ij values for variant pairs
    p_ij_grid = m_pij_grid_vectorised.main(recom_tract_len, lookup_table_rho_vals, merged_eq3_table.copy())

    # Get final pairwise (variant pairs) likelihoods
    interpolated_eq2_df = m_pairwise_rho_estimator_intp_rect_biv.main(merged_eq3_table.copy(),
                                                                      table_ids_for_eq3.copy(),
                                                                      p_ij_grid.copy(),
                                                                      lookup_table,
                                                                      5)
    # Collect pairwise likelihoods across depths

    # Bootstrap to get final results with confidence interval
    bootstrap_results = interpolated_eq2_df.sample(frac=1, replace=True, axis=0)
