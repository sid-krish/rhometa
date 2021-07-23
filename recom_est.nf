#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


def helpMessage() {

    // TO DO: add remaining arguments and descriptions    
    log.info"""
    Usage:

    nextflow run main.nf --bam_file in.bam --reference_genome ref.fa [options]

    Options:
    --bam_file [.bam]    Description
    --reference_genome [.fa]  Description
    --ldpop_rho_range [int,int]   Description
    ... etc.

    `nextflow run main.nf --help` to show this help message
    """.stripIndent()

}


process LOFREQ{
    publishDir "Recom_Est_Output", mode: "copy"

    maxForks 1

    cpus {num_cores}
    // memory '1 GB'

    // executor 'local'
    // time '30m'
    // scratch true
    // queue 'i3q'

    input:
        path reference_fa
        path bam
        val num_cores

    output:
        path "lofreqOut.vcf", emit: lofreqOut_vcf

    script:
    """
    samtools faidx ${reference_fa}
    samtools sort --threads ${num_cores} ${bam} -o Aligned.csorted.bam
    samtools index -@ ${num_cores} Aligned.csorted.bam

    #lofreq call -f ${reference_fa} -o lofreqOut.vcf Aligned.csorted.bam
    #lofreq call-parallel --pp-threads ${num_cores} --no-default-filter -f ${reference_fa} -o lofreqOut.vcf Aligned.csorted.bam
    lofreq call-parallel --pp-threads ${num_cores} -f ${reference_fa} -o lofreqOut.vcf Aligned.csorted.bam
    """
}


process PAIRWISE_TABLE{
    publishDir "Recom_Est_Output", mode: "copy"

    maxForks 1

    echo true

    cpus {num_cores}

    input:
        path aligned_bam
        path vcf_file
        val single_end
        val num_cores
        val window_size

    output:
        path "pairwise_table.pkl", emit: pairwise_table_pkl

    script:
    """
    gen_pairwise_table.py ${single_end} ${aligned_bam} ${vcf_file} ${num_cores} ${window_size}
    """

}


process RECOM_RATE_ESTIMATOR {
    publishDir "Recom_Est_Output", mode: "copy"

    maxForks 1

    cpus {num_cores}

    input:
        path downsampled_lookup_tables
        path pairwise_table_pkl
        val recom_tract_len
        val depth_range
        val n_bootstrap_samples
        val ldpop_rho_range
        val num_cores


    output:
        path "final_results.csv", emit: final_results_csv
        path "final_results_max_vals.csv", emit: final_results_max_vals_csv
        path "final_results_summary.csv", emit: final_results_summary_csv

    script:
    """
    mprr_main_parallel.py ${recom_tract_len} ${depth_range} ${n_bootstrap_samples} ${ldpop_rho_range} ${pairwise_table_pkl} ${num_cores}
    """

}


process FINAL_RESULTS_PLOT {
    publishDir "Recom_Est_Output", mode: "copy"

    maxForks 1

    cpus 1

    input:
        path final_results_csv

    output:
        path "final_results_plot.png", emit: final_results_plot_png

    script:
    """
    final_results_plot.py final_results.csv
    """
}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context
    
    params.help = false
    params.num_cores = 4
    params.recom_tract_len = 500
    params.ldpop_rho_range = "101,100"
    params.window_size = 500 // For single end this is the read size, for paired end this is the max fragment length
    params.single_end = false
    params.depth_range = "3,100" // min_depth, max_depth
    params.n_bootstrap_samples = 20 // number of bootstrap samples to get error bars for final results

    params.bam_file = 'none'
    params.reference_genome = 'none'
    params.lookup_tables = "$baseDir/Lookup_tables/lk_downsampled_*.csv"
    bam_file_channel = Channel.fromPath( params.bam_file )
    reference_genome_channel = Channel.fromPath( params.reference_genome )
    downsampled_lookup_tables = Channel.fromPath( params.lookup_tables ).collect()

    if (params.help) {
        // Show help message from helpMessage() function
        // params.help = false by default
        helpMessage()
        exit 0
    }

    if (params.reference_genome == 'none') {
        println "No input .fa specified. Use --reference_genome [.fa]"
        exit 1
    }

    if (params.bam_file == 'none') {
        println "No input .bam specified. Use --bam_file [.bam]"
        exit 1
    }

    LOFREQ(reference_genome_channel, bam_file_channel, params.num_cores)

    // Bams need to be query name sorted.
    PAIRWISE_TABLE(bam_file_channel, LOFREQ.out.lofreqOut_vcf, params.single_end, params.num_cores, params.window_size)

    RECOM_RATE_ESTIMATOR(downsampled_lookup_tables, PAIRWISE_TABLE.out.pairwise_table_pkl, params.recom_tract_len, params.depth_range, params.n_bootstrap_samples, params.ldpop_rho_range, params.num_cores)

    FINAL_RESULTS_PLOT(RECOM_RATE_ESTIMATOR.out.final_results_csv)

}