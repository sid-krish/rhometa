#!/usr/bin/env python

import sys
import pandas as pd
from ldpop import rhos_from_string

import m_isolate_by_depth
import m_biallelic_filter_pairwise_table
import m_pairwise_lookup_format_pyrho
import m_custom_hap_sets_and_merge
import m_pij_grid_vectorised
import m_pairwise_rho_estimator_intp_rect_biv

if __name__ == '__main__':
    recom_tract_len = int(sys.argv[1])
    depth_range = sys.argv[2]
    n_resamples = int(sys.argv[3])
    lookup_table_rho_range = sys.argv[4]
    pairwise_table_file = sys.argv[5]

    # recom_tract_len = 500
    # depth_range = "3,200"
    # n_resamples = 50
    # lookup_table_rho_range = "101,100"
    # pairwise_table_file = "pairwise_table.pkl"

    pairwise_table = pd.read_pickle(pairwise_table_file)

    lookup_table_rho_vals = rhos_from_string(lookup_table_rho_range)

    depth_lower_limit, depth_upper_limit = [int(i) for i in depth_range.split(',')]

    results_across_depths = pd.DataFrame()  # Pairwise likelihoods across depths

    # Step 1: Isolate pairwise entries based on depth
    # The pipeline isolates pairwise entries based on depth and uses the appropriate lookup table for each depth.
    # note: after isolation some files will be empty
    for depth in range(depth_lower_limit, depth_upper_limit + 1):
        lookup_table = f"../../../Lookup_tables/lk_downsampled_{depth}.csv"
        pairwise_table_slice = m_isolate_by_depth.main(pairwise_table,
                                                       depth).copy()  # copy to create a new df using slice

        # Step 3: Filter pairwise table so that it is bi-allelic pairs only
        # note: after bi-allelic filtering some files will be empty
        if not pairwise_table_slice.empty:
            pairwise_biallelic_table = m_biallelic_filter_pairwise_table.main(pairwise_table_slice).copy()

        else:
            continue  # skip to next iteration

        # Step 4: convert to lookup format to match against likelihood tables
        if not pairwise_biallelic_table.empty:
            lookup_formatted_table = m_pairwise_lookup_format_pyrho.main(pairwise_biallelic_table)

        else:
            continue

        # Step 5: Merge lookup formatted table on likelihood table
        # this step is equation 3 as defined in PIIM paper
        merged_eq3_table, table_ids_for_eq3 = m_custom_hap_sets_and_merge.main(pairwise_biallelic_table,
                                                                               lookup_formatted_table,
                                                                               lookup_table_rho_vals, lookup_table)

        # Step 6: Calculate p_ij values for variant pairs
        p_ij_grid = m_pij_grid_vectorised.main(recom_tract_len, lookup_table_rho_vals, merged_eq3_table)

        # Step 7: Get final pairwise (variant pairs) likelihoods
        interpolated_eq2_df = m_pairwise_rho_estimator_intp_rect_biv.main(merged_eq3_table,
                                                                          table_ids_for_eq3,
                                                                          p_ij_grid, lookup_table,
                                                                          depth)

        # Step 8: Collect pairwise likelihoods across depths
        results_across_depths = results_across_depths.append(interpolated_eq2_df, ignore_index=True)

    # Step 9: Bootstrap to get final results with confidence interval
    # NOTE: the bootstrapping process can be parallelized using dask
    final_results_df = pd.DataFrame()
    final_results_max_rho_and_likelihoods_df = pd.DataFrame()

    # use numpy to generate to generate a array of numbers based on a seed val for sample random_state
    for i in range(n_resamples):
        # the resampled df will have the same num of rows as the original
        bootstrap_results = results_across_depths.sample(frac=1, replace=True, random_state=i, axis=0)

        bootstrap_results_sums = bootstrap_results.sum(axis=0)

        bootstrap_results_finalised = pd.DataFrame()
        bootstrap_results_finalised["rho"] = lookup_table_rho_vals
        bootstrap_results_finalised["likelihood_sums"] = list(bootstrap_results_sums)
        bootstrap_results_finalised["bootstrap_sample"] = i

        final_results_df = final_results_df.append(bootstrap_results_finalised, ignore_index=True)

        # Get max rho and likelihood for each bootstrap
        bootstrap_results_finalised.sort_values(by=["likelihood_sums", "rho"], ascending=[False, True], inplace=True)

        final_results_max_rho_and_likelihoods_df = final_results_max_rho_and_likelihoods_df.append(
            bootstrap_results_finalised.iloc[0], ignore_index=True)

    # Step 10: Export results
    final_results_df.to_csv("final_results.csv", index=False)

    # set appropriate column names indicating they are max values
    final_results_max_rho_and_likelihoods_df.columns = ["bootstrap_sample", "max_lk", "max_rho"]
    # reorder columns
    final_results_max_rho_and_likelihoods_df = final_results_max_rho_and_likelihoods_df[
        ["max_rho", "max_lk", "bootstrap_sample"]]
    final_results_max_rho_and_likelihoods_df.to_csv("final_results_max_vals.csv", index=False)

    # summary statistics on 'max_rho' column
    summary_stats = final_results_max_rho_and_likelihoods_df['max_rho'].describe()
    summary_stats = pd.DataFrame(summary_stats).T

    summary_stats = summary_stats.rename(columns={"count": "bootstrap_samples"})

    summary_stats.to_csv("final_results_summary.csv")
