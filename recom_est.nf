#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


def helpMessage() {

    log.info"""
    Description:

    Usage:
    nextflow run recom_est.nf --bam_file in.bam --reference_genome ref.fa [options]

    Help:
    nextflow run recom_est.nf --help

    Required:
    --bam_file [*.bam], Query name sorted bam file
    --reference_genome [*.fa],  Single genome fasta file
    --downsampled_lookup_tables [dir], default:[Lookup_tables], Folder containing downsampled lookup tables generated by lookup_table_gen.nf

    Options:
    --ldpop_rho_range [int,int], default:[101,100], The range of rho values used to generate lookup tables
    --recom_tract_len [int], default:[500], Recombination tract length to use
    --window_size [int], default:[500], Window size for variant pairs. For single end this is the read size, for paired end this is the max fragment length
    --single_end, Used for single end read bams
    --subsample_bam, Used for when read depths are higher than what can be analysed with available lookup tables (downsample bam to match max lookup table depth)
    --depth_range [int,int], default:[3,100], Minimum and maximum depth downsampled lookup tables available. Minimum should be no less than 3
    --n_bootstrap_samples [int], default:[20], Number of bootstrap samples to get confidence interval for recombination rate estimate

    """.stripIndent()

}


process SUBSAMPLE_BAM{
    publishDir "Recom_Est_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

    maxForks 1

    input:
        val depth_range
        path bam
        val prepend_filename
    
    output:
        path "subsampled.bam", emit: subsampled_bam

    script:
    """
    samtools sort ${bam} > Aligned_sorted.bam
    samtools mpileup Aligned_sorted.bam > Aligned_sorted.pileup

    subsample_bam.py Aligned_sorted.pileup ${depth_range} ${bam}
    """
}


process LOFREQ{
    publishDir "Recom_Est_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

    maxForks 1

    input:
        path reference_fa
        path bam
        val prepend_filename

    output:
        path "lofreqOut.vcf", emit: lofreqOut_vcf

    script:
    """
    samtools faidx ${reference_fa}
    samtools sort --threads $task.cpus ${bam} -o Aligned.csorted.bam
    samtools index -@ $task.cpus Aligned.csorted.bam

    #lofreq call -f ${reference_fa} -o lofreqOut.vcf Aligned.csorted.bam
    #lofreq call-parallel --pp-threads $task.cpus --no-default-filter -f ${reference_fa} -o lofreqOut.vcf Aligned.csorted.bam
    lofreq call-parallel --pp-threads $task.cpus -f ${reference_fa} -o lofreqOut.vcf Aligned.csorted.bam
    """
}


process PAIRWISE_TABLE{
    publishDir "Recom_Est_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

    maxForks 1

    input:
        path aligned_bam
        path vcf_file
        val single_end
        val window_size
        val prepend_filename

    output:
        path "pairwise_table.pkl", emit: pairwise_table_pkl

    script:
    """
    gen_pairwise_table.py ${single_end} ${aligned_bam} ${vcf_file} $task.cpus ${window_size}
    """
}


process RECOM_RATE_ESTIMATOR {
    publishDir "Recom_Est_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

    maxForks 1

    input:
        path downsampled_lookup_tables
        path pairwise_table_pkl
        val recom_tract_len
        val depth_range
        val n_bootstrap_samples
        val ldpop_rho_range
        val prepend_filename


    output:
        path "final_results.csv", emit: final_results_csv
        path "final_results_max_vals.csv", emit: final_results_max_vals_csv
        path "final_results_summary.csv", emit: final_results_summary_csv

    script:
    """
    mprr_main_parallel.py ${recom_tract_len} ${depth_range} ${n_bootstrap_samples} ${ldpop_rho_range} ${pairwise_table_pkl} $task.cpus
    """

}


process FINAL_RESULTS_PLOT {
    publishDir "Recom_Est_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

    maxForks 1

    input:
        path final_results_csv
        val prepend_filename

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

    // Params
    params.help = false
    params.subsample_bam = false
    params.prepend_filename = ""
    params.recom_tract_len = 500
    params.ldpop_rho_range = "0,0.01,1,1,100"
    params.window_size = 300 // For single end this is the read size, for paired end this is the max insert length
    params.single_end = false
    params.depth_range = "3,200" // min_depth, max_depth
    params.n_bootstrap_samples = 25 // number of bootstrap samples to get error bars for final results

    params.bam_file = 'none'
    params.reference_genome = 'none'
    params.lookup_tables = "Lookup_tables"

    // Channels
    bam_file_channel = Channel.fromPath( params.bam_file )
    reference_genome_channel = Channel.fromPath( params.reference_genome )
    downsampled_lookup_tables = Channel.fromPath( "${params.lookup_tables}/lk_downsampled_*.csv" ).collect()

    // Input verification
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

    // Process execution
    if (params.subsample_bam) {
        // Bams need to be query name sorted.
        SUBSAMPLE_BAM(params.depth_range, bam_file_channel, params.prepend_filename)

        LOFREQ(reference_genome_channel, SUBSAMPLE_BAM.out.subsampled_bam, params.prepend_filename)

        PAIRWISE_TABLE(SUBSAMPLE_BAM.out.subsampled_bam, LOFREQ.out.lofreqOut_vcf, params.single_end, params.window_size, params.prepend_filename)

    }

    else {
        // Bams need to be query name sorted.
        LOFREQ(reference_genome_channel, bam_file_channel, params.prepend_filename)

        PAIRWISE_TABLE(bam_file_channel, LOFREQ.out.lofreqOut_vcf, params.single_end, params.window_size, params.prepend_filename)
    }

    RECOM_RATE_ESTIMATOR(downsampled_lookup_tables, PAIRWISE_TABLE.out.pairwise_table_pkl, params.recom_tract_len, params.depth_range, params.n_bootstrap_samples, params.ldpop_rho_range, params.prepend_filename)

    FINAL_RESULTS_PLOT(RECOM_RATE_ESTIMATOR.out.final_results_csv, params.prepend_filename)

}