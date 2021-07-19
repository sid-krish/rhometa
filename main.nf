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

if (params.help) {
    // Show help message above 
    // In .config, params.help = false by default
    helpMessage()
    exit 0
}

// Deal with missing (mandatory) arguments 

if (params.bam_file == 'none') {
    println "No input .bam specified. Use --bam_file [.bam]"
    exit 1
}

process LOFREQ{
    publishDir "Output", mode: "copy"

    maxForks 1

    // cpus 4
    // memory '5 GB'

    // executor 'pbspro'
    // time '30m'
    // scratch true
    // queue 'i3q'

    input:
        path reference_genome_file
        path bam_file
        
    output:
        path "lofreqOut.vcf", emit: lofreqOut_vcf

    script:
    """
    samtools faidx ${reference_genome_file}
    samtools sort --threads 4 ${bam_file} -o chrom_sorted.bam
    samtools index -@ 4 chrom_sorted.bam

    #lofreq call -f ${reference_genome_file} -o lofreqOut.vcf chrom_sorted.bam
    #lofreq call-parallel --pp-threads 4 --no-default-filter -f ${reference_genome_file} -o lofreqOut.vcf chrom_sorted.bam
    lofreq call-parallel --pp-threads 4 -f ${reference_genome_file} -o lofreqOut.vcf chrom_sorted.bam
    """
}


process PAIRWISE_TABLE{
    publishDir "Output", mode: "copy"

    maxForks 1

    input:
        path vcf_file
        path aligned_bam
        

    output:
        path "pairwise_table.pkl", emit: pairwise_table_pkl

    script:
    """
    gen_pairwise_table.py ${params.seq_type} ${aligned_bam} ${vcf_file}
    """

}


process RECOM_RATE_ESTIMATOR {
    publishDir "Output", mode: "copy"

    maxForks 1

    input:
        // TODO: need to use collect files to get lookup tables per depth
        path pairwise_table_pkl

    output:
        path "final_results.csv", emit: final_results_csv
        path "final_results_max_vals.csv", emit: final_results_max_vals_csv
        path "final_results_summary.csv", emit: final_results_summary_csv

    script:
    """
    mprr_main.py ${params.recom_tract_len} ${params.depth_range} ${params.n_bootstrap_samples} ${params.ldpop_rho_range} ${pairwise_table_pkl}
    """

}

process FINAL_RESULTS_PLOT {
    publishDir "Output", mode: "copy"

    maxForks 1

    input:
        path final_results_csv

    output:
        path "final_results_plot.png", emit: final_results_plot_png

    script:
    """
    final_results_plot.py final_results.csv
    """
}

// parameters

// params.bam_file = ""

// input validation
// if (params.bam_file) {

// } else {
//     exit 1, 'Missing bam file'
// }


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // TODO: ALL PARAMS MUST BE SENT IN AS AN INPUT CHANNEL CAUSES ISSUES OTHERWISE
    
    bam_file_channel = Channel.fromPath( params.bam_file )
    reference_genome_channel = Channel.fromPath( params.reference_genome )

    LOFREQ(reference_genome_channel, bam_file_channel)

    PAIRWISE_TABLE(LOFREQ.out.lofreqOut_vcf, bam_file_channel)

    RECOM_RATE_ESTIMATOR(PAIRWISE_TABLE.out.pairwise_table_pkl)

    FINAL_RESULTS_PLOT(RECOM_RATE_ESTIMATOR.out.final_results_csv)

}
