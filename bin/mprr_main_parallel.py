#!/usr/bin/env python

import sys
import pandas as pd
from multiprocessing import Pool
from ldpop import rhos_from_string

import m_isolate_by_depth
import m_biallelic_filter_pairwise_table
import m_pairwise_lookup_format_pyrho
import m_custom_hap_sets_and_merge
import m_pij_grid_vectorised
import m_pairwise_rho_estimator_intp_rect_biv


def steps_1_to_7(arg_list):
    depth = arg_list[0]
    recom_tract_len = arg_list[1]
    pairwise_table = arg_list[2]
    lookup_table_rho_vals = arg_list[3]

    # Step 1: Isolate pairwise entries based on depth
    # The pipeline isolates pairwise entries based on depth and uses the appropriate lookup table for each depth.
    # note: after isolation some files will be empty

    lookup_table = f"lk_downsampled_{depth}.csv"
    pairwise_table_slice = m_isolate_by_depth.main(pairwise_table,
                                                   depth).copy()  # copy to create a new df using slice

    # Step 3: Filter pairwise table so that it is bi-allelic pairs only
    # note: after bi-allelic filtering some files will be empty
    if not pairwise_table_slice.empty:
        pairwise_biallelic_table = m_biallelic_filter_pairwise_table.main(pairwise_table_slice).copy()

    else:
        return None

    # Step 4: convert to lookup format to match against likelihood tables
    if not pairwise_biallelic_table.empty:
        lookup_formatted_table = m_pairwise_lookup_format_pyrho.main(pairwise_biallelic_table)

    else:
        return None

    # Step 5: Merge lookup formatted table on likelihood table
    # this step is equation 3 as defined in PIIM paper
    merged_eq3_table, table_ids_for_eq3 = m_custom_hap_sets_and_merge.main(pairwise_biallelic_table,
                                                                           lookup_formatted_table,
                                                                           lookup_table_rho_vals,
                                                                           lookup_table)

    # Step 6: Calculate p_ij values for variant pairs
    p_ij_grid = m_pij_grid_vectorised.main(recom_tract_len, lookup_table_rho_vals, merged_eq3_table)

    # Step 7: Get final pairwise (variant pairs) likelihoods
    interpolated_eq2_df = m_pairwise_rho_estimator_intp_rect_biv.main(merged_eq3_table,
                                                                      table_ids_for_eq3,
                                                                      p_ij_grid,
                                                                      lookup_table,
                                                                      depth)

    return interpolated_eq2_df


def step_9(arg_list):
    nth_resample = arg_list[0]
    results_across_depths = arg_list[1]
    lookup_table_rho_vals = arg_list[2]

    # use numpy to generate to generate a array of numbers based on a seed val for sample random_state
    # currently uses the nth_resample value as each one is unique

    #the resampled df will have the same num of rows as the original
    bootstrap_results = results_across_depths.sample(frac=1, replace=True, random_state=nth_resample, axis=0)

    bootstrap_results_sums = bootstrap_results.sum(axis=0)

    bootstrap_results_finalised = pd.DataFrame()
    bootstrap_results_finalised["rho"] = lookup_table_rho_vals
    bootstrap_results_finalised["likelihood_sums"] = list(bootstrap_results_sums)
    bootstrap_results_finalised["bootstrap_sample"] = nth_resample

    return bootstrap_results_finalised


if __name__ == '__main__':
    recom_tract_len = int(sys.argv[1])
    depth_range = sys.argv[2]
    n_resamples = int(sys.argv[3])
    lookup_table_rho_range = sys.argv[4]
    pairwise_table_file = sys.argv[5]
    num_cores = int(sys.argv[6])


    # recom_tract_len = 500
    # depth_range = "3,100"
    # n_resamples = 20
    # lookup_table_rho_range = "101,100"
    # pairwise_table_file = "../Output/rho_10_sam_10_gen_10000/rho_10_sam_10_gen_10000_pairwise_table.pkl"
    # num_cores = 4

    pairwise_table = pd.read_pickle(pairwise_table_file)

    lookup_table_rho_vals = rhos_from_string(lookup_table_rho_range)

    depth_lower_limit, depth_upper_limit = [int(i) for i in depth_range.split(',')]



    depth_vals = [i for i in range(depth_lower_limit, depth_upper_limit + 1)]


    steps_1_to_7_arg_list = []
    for i in depth_vals:
        steps_1_to_7_arg_list.append([i, recom_tract_len, pairwise_table, lookup_table_rho_vals])

    # Parallel execution of step 1 through 7
    with Pool(processes=num_cores) as p:
        steps_1_to_7_results = p.map(steps_1_to_7, steps_1_to_7_arg_list)

    # Step 8: Collect pairwise likelihoods across depths
    results_across_depths = pd.DataFrame()  # Pairwise likelihoods across depths
    for i in steps_1_to_7_results:
        results_across_depths = results_across_depths.append(i, ignore_index=True)

    # Step 9: Bootstrap to get final results with confidence interval
    steps_9_arg_list = []
    resample_range = [i for i in range(n_resamples)]
    for i in resample_range:
        steps_9_arg_list.append([i, results_across_depths, lookup_table_rho_vals])

    # Parallel execution of step 9
    with Pool(processes=num_cores) as p:
        steps_9_results = p.map(step_9, steps_9_arg_list)

    final_results_df = pd.DataFrame()
    final_results_max_rho_and_likelihoods_df = pd.DataFrame()

    for i in steps_9_results:
        final_results_df = final_results_df.append(i, ignore_index=True)

        # Get max rho and likelihood for each bootstrap
        i = i.sort_values(by=["likelihood_sums", "rho"], ascending=[False, True])
        final_results_max_rho_and_likelihoods_df = final_results_max_rho_and_likelihoods_df.append(i.iloc[0],
                                                                                                   ignore_index=True)

    # Step 10: Export results
    final_results_df.to_csv("final_results.csv", index=False)

    # set appropriate column names indicating they are max values
    final_results_max_rho_and_likelihoods_df.columns = ["bootstrap_sample", "max_lk", "max_rho"]
    # reorder columns
    final_results_max_rho_and_likelihoods_df = final_results_max_rho_and_likelihoods_df[
        ["max_rho", "max_lk", "bootstrap_sample"]]
    final_results_max_rho_and_likelihoods_df.to_csv("final_results_max_vals.csv", index=False)

    # summary statistics on 'max_rho' column
    summary_stats = final_results_max_rho_and_likelihoods_df['max_rho'].describe(percentiles = [.05, .25, .5, .75, .95])
    summary_stats = pd.DataFrame(summary_stats).T

    summary_stats = summary_stats.rename(columns={"count": "bootstrap_samples"})

    summary_stats.to_csv("final_results_summary.csv")
