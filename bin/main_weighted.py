#!/usr/bin/env python

import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

from ldpop import rhos_from_string
from rmeta_main import biallelic_filter_pairwise_table
from rmeta_main import custom_hap_sets_and_merge
from rmeta_main import isolate_by_depth
from rmeta_main import pairwise_lookup_format_pyrho
from rmeta_main import pairwise_rho_estimator_intp_rect_biv
from rmeta_main import pij_grid_vectorised


def module_calls(arg_list):
    depth = arg_list[0]
    recom_tract_len = arg_list[1]
    pairwise_table = arg_list[2]
    lookup_table_rho_vals = arg_list[3]

    # Step 1: Isolate pairwise entries based on depth
    # The pipeline isolates pairwise entries based on depth and uses the appropriate lookup table for each depth.
    # note: after isolation some files will be empty

    lookup_table = f"lk_downsampled_{depth}.csv"
    # lookup_table = f"/Volumes/Backup/Lookup_tables/Lookup_tables_stp/lk_downsampled_{depth}.csv"
    pairwise_table_slice = isolate_by_depth.main(pairwise_table, depth)

    # Step 2: Filter pairwise table so that it is bi-allelic pairs only
    # note: after bi-allelic filtering some files will be empty
    # copy to create a new df using slice
    if not pairwise_table_slice.empty:
        pairwise_biallelic_table = biallelic_filter_pairwise_table.main(
            pairwise_table_slice.copy()
        )

    else:
        return None

    # Step 3: convert to lookup format to match against likelihood tables
    if not pairwise_biallelic_table.empty:
        lookup_formatted_table = pairwise_lookup_format_pyrho.main(
            pairwise_biallelic_table.copy()
        )

    else:
        return None

    # Step 4: Merge lookup formatted table on likelihood table
    merged_eq3_table, table_ids_for_eq3 = custom_hap_sets_and_merge.main(
        pairwise_biallelic_table.copy(),
        lookup_formatted_table.copy(),
        lookup_table_rho_vals,
        lookup_table,
    )

    # Step 5: Calculate p_ij values for variant pairs
    p_ij_grid = pij_grid_vectorised.main(
        recom_tract_len, lookup_table_rho_vals, merged_eq3_table.copy()
    )

    # Step 6: Get final pairwise (variant pairs) likelihoods
    interpolated_eq2_df = pairwise_rho_estimator_intp_rect_biv.main(
        merged_eq3_table.copy(),
        table_ids_for_eq3.copy(),
        p_ij_grid.copy(),
        lookup_table,
        depth,
    )

    # Step 7: Performing weighting of log_likelihoods
    # 7.1: Convert log-likelihoods to likelihoods
    interpolated_eq2_np = interpolated_eq2_df.to_numpy(dtype=np.float128, copy=True)
    exp_np = np.exp(interpolated_eq2_np)

    # 7.2: We normalise the likelihood values such that the values add up to 1. We do this by dividing each value by the sum of all values.
    exp_norm_np = exp_np / np.sum(exp_np)

    # 7.3: Multiply the values by the depth and the number of observations.
    # This scaling is the “weighting” procedure, the higher the number of observations or depth the more it will be weighted for the final log sum calculation.
    observations = exp_norm_np.shape[0]
    exp_norm_scaled_np = depth * observations * exp_norm_np

    # 7.4: Following this weighting step, we convert back to log space by taking the log of the values
    weighted_log_np = np.log(exp_norm_scaled_np)

    # 7.5: Back to df
    weighted_log_df = pd.DataFrame(weighted_log_np)

    return weighted_log_df


if __name__ == "__main__":
    recom_tract_len = int(sys.argv[1])
    depth_range = sys.argv[2]
    lookup_table_rho_range = sys.argv[3]
    pairwise_table_file = sys.argv[4]
    num_cores = int(sys.argv[5])
    genome_size  = int(sys.argv[6])

    # recom_tract_len = 1000
    # depth_range = "3,250"
    # lookup_table_rho_range = "0,0.01,1,1,100"
    # pairwise_table_file = "../Recom_Est_Output/123_rho_0.005_theta_0.005_sample_size_20_depth_4_genome_size_5000_seed_1_final_pairwise_table.pkl"
    # num_cores = 4

    pairwise_table = pd.read_pickle(pairwise_table_file)

    lookup_table_rho_vals = rhos_from_string(lookup_table_rho_range)

    depth_lower_limit, depth_upper_limit = [int(i) for i in depth_range.split(",")]

    depth_vals = [i for i in range(depth_lower_limit, depth_upper_limit + 1)]

    module_calls_arg_list = [
        [i, recom_tract_len, pairwise_table, lookup_table_rho_vals] for i in depth_vals
    ]

    # Parallel execution of steps 1 through 7
    with Pool(processes=num_cores) as p:
        module_calls_results = p.map(module_calls, module_calls_arg_list)

    # Step 8: Collect pairwise likelihoods across depths
    results_across_depths = pd.concat(module_calls_results, ignore_index=True)

    # Step 9: Export results
    log_sums = results_across_depths.sum(axis=0)

    results = pd.DataFrame()
    results["rho"] = [float("%.5g" % i) for i in lookup_table_rho_vals]
    results["log_likelihood_sum"] = [float("%.10g" % i) for i in log_sums]

    results.to_csv("log_likelihood_sums.csv", index=False)

    rho_estimate = results.iloc[results["log_likelihood_sum"].idxmax()]
    rho_estimate = rho_estimate.to_frame().T
    rho_estimate["genome_length"] = genome_size
    rho_estimate["recom_tract_len"] = recom_tract_len
    rho_estimate["rho_per_site"] = rho_estimate["rho"] / rho_estimate["recom_tract_len"]
    rho_estimate.to_csv("rho_estimate.csv", index=False)
