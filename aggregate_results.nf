#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process RECOM_RESULTS {
    input:
        path(recom_results)

    output:
        path("recom_results.csv")

    script:
    """
    aggregate_recom_results.py ${recom_results}
    """
}


process THETA_RESULTS {
    input:
        path(theta_results)

    output:
        path("theta_results.csv")

    script:
    """
    aggregate_theta_results.py ${theta_results}
    """
}


process PLOT_RECOM_RESULTS {
    publishDir "agg_results", mode: "copy"

    input:
        path(agg_recom_results)

    output:
        path("agg_recom_results_plot.png")

    script:
    """
    agg_recom_results_plot.py ${agg_recom_results}
    """
}


process PLOT_THETA_RESULTS {
    publishDir "agg_results", mode: "copy"

    input:
        path(agg_theta_results)

    output:
        path("agg_theta_results_plot.png")

    script:
    """
    agg_theta_results_plot.py ${agg_theta_results}
    """
}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // For each process there is a output of tuple with the necessary files/values to move forward until they are no longer need

    // Params
    params.recom_summary = "./s_pnemonia/Recom_Est_Output/*final_results_summary.csv"
    params.theta_summary = "./s_pnemonia/Theta_Est_Output/*Theta_estimate_stats.csv"

    // Channels
    recom_summary_channel = Channel.fromPath( params.recom_summary )
    theta_summary_channel = Channel.fromPath( params.theta_summary )

    // Process execution
    RECOM_RESULTS(recom_summary_channel) | collectFile(name:"aggregate_recom_results.csv", keepHeader:true, storeDir:"agg_results")
    
    PLOT_RECOM_RESULTS(Channel.fromPath( "agg_results/aggregate_recom_results.csv" ))

    THETA_RESULTS(theta_summary_channel) | collectFile(name:"aggregate_theta_results.csv", keepHeader:true, storeDir:"agg_results")

    PLOT_THETA_RESULTS(Channel.fromPath( "agg_results/aggregate_theta_results.csv" ))

}