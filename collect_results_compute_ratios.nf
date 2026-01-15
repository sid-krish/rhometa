#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process COLLECT_RESULTS {
    // debug true

    publishDir params.final_results_out_dir, mode: 'copy'

    input:
    path theta_files
    path rho_files

    output:
    path "collected_merged_results.csv"

    script:
    """
    # Collect files from nextflow input and merge results from theta and rho files
    collect_results.py
    """
}

process COMPUTE_RATIOS {
    publishDir params.final_results_out_dir, mode: 'copy'

    input:
    path results_file

    output:
    path "collected_merged_final.csv"

    script:
    """
    compute_ratios.py ${results_file}
    """
}

workflow {
    params.theta_dir = "Theta_Est_Output/theta_estimate"
    params.rho_dir = "Rho_Est_Output/rho_estimate"

    params.final_results_out_dir = "Final_Output"

    theta_files = Channel.fromPath(params.theta_dir + "/*filtered_angsd_theta_final.csv", checkIfExists: true).collect()
    rho_files = Channel.fromPath(params.rho_dir + "/*rho_estimate.csv", checkIfExists: true).collect()

    COLLECT_RESULTS(theta_files, rho_files)
    COMPUTE_RATIOS(COLLECT_RESULTS.out)
}