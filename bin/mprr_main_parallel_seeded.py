#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np
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

    print('Processing steps 1-7 for depth: {}'.format(depth))

    lookup_table = f"lk_downsampled_{depth}.csv"
    pairwise_table_slice = m_isolate_by_depth.main(pairwise_table,
                                                   depth)

    # Step 3: Filter pairwise table so that it is bi-allelic pairs only
    # note: after bi-allelic filtering some files will be empty
    # copy to create a new df using slice
    if not pairwise_table_slice.empty:
        print('Slice at {} contained {} pairs'.format(depth, len(pairwise_table_slice)))
        pairwise_biallelic_table = m_biallelic_filter_pairwise_table.main(pairwise_table_slice.copy())

    else:
        print('Slice at {} was empty'.format(depth))
        return None 

    # Step 4: convert to lookup format to match against likelihood tables
    if not pairwise_biallelic_table.empty:
        lookup_formatted_table = m_pairwise_lookup_format_pyrho.main(pairwise_biallelic_table.copy())

    else:
        return None

    # Step 5: Merge lookup formatted table on likelihood table
    merged_eq3_table, table_ids_for_eq3 = m_custom_hap_sets_and_merge.main(pairwise_biallelic_table.copy(),
                                                                           lookup_formatted_table.copy(),
                                                                           lookup_table_rho_vals,
                                                                           lookup_table)

    # Step 6: Calculate p_ij values for variant pairs
    p_ij_grid = m_pij_grid_vectorised.main(recom_tract_len, lookup_table_rho_vals, merged_eq3_table.copy())

    # Step 7: Get final pairwise (variant pairs) likelihoods
    interpolated_eq2_df = m_pairwise_rho_estimator_intp_rect_biv.main(merged_eq3_table.copy(),
                                                                      table_ids_for_eq3.copy(),
                                                                      p_ij_grid.copy(),
                                                                      lookup_table,
                                                                      depth)

    return interpolated_eq2_df


def step_9(arg_list):
    rng_val = arg_list[0]
    results_across_depths = arg_list[1]
    lookup_table_rho_vals = arg_list[2]

    #the resampled df will have the same num of rows as the original
    bootstrap_results = results_across_depths.sample(frac=1, replace=True, random_state=rng_val, axis=0)

    bootstrap_results_sums = bootstrap_results.sum(axis=0)

    bootstrap_results_finalised = pd.DataFrame()
    bootstrap_results_finalised["rho"] = lookup_table_rho_vals
    bootstrap_results_finalised["likelihood_sums"] = list(bootstrap_results_sums)
    bootstrap_results_finalised["bootstrap_sample"] = rng_val

    return bootstrap_results_finalised


if __name__ == '__main__':
    recom_tract_len = int(sys.argv[1])
    depth_range = sys.argv[2]
    n_resamples = int(sys.argv[3])
    lookup_table_rho_range = sys.argv[4]
    pairwise_table_file = sys.argv[5]
    num_cores = int(sys.argv[6])
    seed_val = int(sys.argv[7])


    # recom_tract_len = 1000
    # depth_range = "3,250"
    # n_resamples = 50
    # lookup_table_rho_range = "0,0.01,1,1,100"
    # pairwise_table_file = "Recom_Est_Output/rho_0.001_theta_0.01_sample_size_25_depth_4_genome_size_100000_seed_0_final_pairwise_table.pkl"
    # num_cores = 4

    pairwise_table = pd.read_pickle(pairwise_table_file)

    lookup_table_rho_vals = rhos_from_string(lookup_table_rho_range)

    depth_lower_limit, depth_upper_limit = [int(i) for i in depth_range.split(',')]

    depth_vals = [_dp for _dp in range(depth_lower_limit, depth_upper_limit + 1)]

    steps_1_to_7_arg_list = []
    for _dp in depth_vals:
        print('Setting up for depth: {}'.format(_dp))
        steps_1_to_7_arg_list.append([_dp, recom_tract_len, pairwise_table, lookup_table_rho_vals])

    # Parallel execution of step 1 through 7
    with Pool(processes=num_cores) as p:
        steps_1_to_7_results = p.map(steps_1_to_7, steps_1_to_7_arg_list)

    # Step 8: Collect pairwise likelihoods across depths
    results_across_depths = pd.concat(steps_1_to_7_results, ignore_index=True)

    # Step 9: Bootstrap to get final results with confidence interval
    steps_9_arg_list = []
    
    rng = np.random.default_rng(seed=seed_val)

    resample_vals = rng.integers(low=0,high=100, size=n_resamples)

    for _rv in resample_vals:
        steps_9_arg_list.append([_rv, results_across_depths, lookup_table_rho_vals])

    # Parallel execution of step 9
    with Pool(processes=num_cores) as p:
        steps_9_results = p.map(step_9, steps_9_arg_list)

    # combine bootstrap results
    final_results_df = pd.concat(steps_9_results, ignore_index=True)
    # for each bootstrap, determine rho at max likeilhood
    ix_max = final_results_df.groupby('bootstrap_sample').idxmax()['likelihood_sums']
    final_results_max_rho_and_likelihoods_df = final_results_df.loc[ix_max]

    # Step 10: Export results
    final_results_df.to_csv("final_results.csv", index=False)

    # set appropriate column names indicating they are max values
    final_results_max_rho_and_likelihoods_df.rename({"rho": 'max_rho', "likelihood_sums": 'max_lk'},
                                                    axis=1, inplace=True)

    # reorder columns
    final_results_max_rho_and_likelihoods_df = final_results_max_rho_and_likelihoods_df[
        ["max_rho", "max_lk", "bootstrap_sample"]]
    final_results_max_rho_and_likelihoods_df.to_csv("final_results_max_vals.csv", index=False)

    # summary statistics on 'max_rho' column
    summary_stats = final_results_max_rho_and_likelihoods_df['max_rho'].describe(percentiles = [.5])
    summary_stats = pd.DataFrame(summary_stats).T
    summary_stats = summary_stats.drop(columns=["min", "max", "std"])
    summary_stats = summary_stats.rename(columns={"count": "num_bootstraps", "mean": "mean_rho_est-full_seq", "50%": "median_rho_est-full_seq"})
    summary_stats["specified_tract_len"] = recom_tract_len
    summary_stats["mean_rho_est-per_site-full_seq_rho_div_by_tract"] = summary_stats["mean_rho_est-full_seq"] / recom_tract_len
    summary_stats["median_rho_est-per_site-full_seq_rho_div_by_tract"] = summary_stats["median_rho_est-full_seq"] / recom_tract_len
    summary_stats = pd.DataFrame(summary_stats).T
    summary_stats.to_csv("final_results_summary.csv", header=False)
