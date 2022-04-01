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
    --bam_file [*.bam], Multi bam support via glob input e.g. "*.bam", quotes but be included for glob. Use with one fasta file only
    --reference_genome [*.fa],  Single/Multi genome fasta file
    --downsampled_lookup_tables [dir], default:[Lookup_tables], Folder containing downsampled lookup tables generated by lookup_table_gen.nf

    Options:
    --ldpop_rho_range [int,int], default:[101,100], The range of rho values used to generate lookup tables
    --recom_tract_len [int], default:[500], Recombination tract length to use
    --window_size [int], default:[500], Window size for variant pairs. For single end this is the read size, for paired end this is the max fragment length
    --single_end, Used for single end read bams
    --depth_range [int,int], default:[3,100], Minimum and maximum depth downsampled lookup tables available. Minimum should be no less than 3
    --n_bootstrap_samples [int], default:[20], Number of bootstrap samples to get confidence interval for recombination rate estimate
    --seed [int] , default:[123], Seed value for samtools subsamping and final bootstrap algorithm. The seed value will be displayed at the start of the output file names

    """.stripIndent()

}


process PREFIX_FILENAME {

    // maxForks 1

    // echo true

    input:
        tuple path(bam),
            path(fasta)
        
        val prefix_fn
        each seed

    output:
        tuple stdout,
            path(bam),
            path(fasta),
            val(seed)

    script:
    """
    prefix_filename_seeded.py ${bam} ${prefix_fn} ${seed}
    """
}


process SUBSAMPLE_BAM{
    // publishDir "Recom_Est_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    // echo true

    input:
        tuple val(prefix_filename),
            path(bam),
            path(fasta),
            val(seed)

            
        val depth_range
    
    output:
        tuple val(prefix_filename),
            path("subsampled.bam"),
            path(fasta),
            val(seed)

    script:
    """
    samtools sort ${bam} > Aligned_sorted.bam
    samtools mpileup Aligned_sorted.bam > Aligned_sorted.pileup

    subsample_bam_seeded.py Aligned_sorted.pileup ${depth_range} Aligned_sorted.bam ${seed}
    """
}


process LOFREQ{
    // publishDir "Recom_Est_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    input:
        tuple val(prefix_filename),
            path(bam),
            path(fasta)

    output:
        tuple val(prefix_filename),
            path(bam),
            path(fasta),
            path("lofreqOut.vcf")

    script:
    """
    samtools faidx ${fasta}
    samtools sort --threads $task.cpus ${bam} -o Aligned.csorted.bam
    samtools index -@ $task.cpus Aligned.csorted.bam

    #lofreq call -f ${fasta} -o lofreqOut.vcf Aligned.csorted.bam
    #lofreq call-parallel --pp-threads $task.cpus --no-default-filter -f ${fasta} -o lofreqOut.vcf Aligned.csorted.bam
    lofreq call-parallel --pp-threads $task.cpus -f ${fasta} -o lofreqOut.vcf Aligned.csorted.bam
    """
}


process FREEBAYES {
    // publishDir "Recom_Est_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    input:
        tuple val(prefix_filename),
            path(bam),
            path(fasta),
            val(seed)

    output:
        tuple val(prefix_filename),
            path(bam),
            path(fasta),
            path("freeBayesOut.vcf"),
            val(seed)

    script:
    """
    samtools faidx ${fasta}
    samtools sort --threads $task.cpus ${bam} -o Aligned.csorted.bam
    samtools index -@ $task.cpus Aligned.csorted.bam

    # only keep SNP type entries
    freebayes -f ${fasta} -p 1 Aligned.csorted.bam | grep -e '^#' -e 'TYPE=snp' > freeBayesOut.vcf
    """
}


process PAIRWISE_TABLE{
    publishDir "Recom_Est_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    input:
        tuple val(prefix_filename),
            path(bam),
            path(fasta),
            path(vcf_file),
            val(seed)

        val single_end
        val window_size

    output:
        tuple val(prefix_filename),
            path("pairwise_table.pkl"),
            val(seed)

    script:
    // -n Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates.
    """
    samtools sort -n --threads $task.cpus ${bam} -o qsorted.bam
    gen_pairwise_table.py ${single_end} qsorted.bam ${vcf_file} $task.cpus ${window_size}
    """
}


process RECOM_RATE_ESTIMATOR {
    publishDir "Recom_Est_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    echo true

    // maxForks 1

    input:
        tuple val(prefix_filename),
            path(pairwise_table_pkl),
            val(seed)

        path downsampled_lookup_tables
        val recom_tract_len
        val depth_range
        val n_bootstrap_samples
        val ldpop_rho_range


    output:
        tuple val(prefix_filename),
            path("final_results.csv"),
            path("final_results_max_vals.csv")
            // path("final_results_summary.csv")

    script:
    """
    main_weighted.py ${recom_tract_len} ${depth_range} ${n_bootstrap_samples} ${ldpop_rho_range} ${pairwise_table_pkl} $task.cpus ${seed}
    """

}


process FINAL_RESULTS_PLOT {
    publishDir "Recom_Est_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    input:
        tuple val(prefix_filename),
            path("final_results.csv"),
            path("final_results_max_vals.csv")
            // path("final_results_summary.csv")

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

    // For each process there is a output of tuple with the necessary files/values to move forward until they are no longer need

    // Params
    params.help = false
    params.seed = [123] // used for samtools subsamping and final bootstrap algorithm
    params.prefix_filename = "none"
    params.recom_tract_len = 1000
    params.ldpop_rho_range = "0,0.01,1,1,100"
    params.window_size = 1000 // For single end this is the read size, for paired end this is the max insert length (1000bp is a practical upper limit)
    params.single_end = false
    params.depth_range = "3,200" // min_depth, max_depth
    params.n_bootstrap_samples = 50 // number of bootstrap samples to get error bars for final results

    params.bam_file = 'none'
    params.reference_genome = 'none'
    
    params.lookup_tables = "/Volumes/Backup/Lookup_tables/Lookup_tables_stp"
    // params.lookup_tables = "/shared/homes/11849395/Lookup_tables/Lookup_tables_stp"
    // params.lookup_tables = "/shared/homes/11849395/lookup_table_gen/Lookup_tables(0.00126)" // hpylori
    // params.lookup_tables = "/shared/homes/11849395/lookup_table_gen/Lookup_tables(0.00002)" // s_pne 5ng
    // params.lookup_tables = "/shared/homes/11849395/lookup_table_gen/Lookup_tables(0.00003)" // s_pne exp1_500ng
    // params.lookup_tables = "/shared/homes/11849395/lookup_table_gen/Lookup_tables(0.00002)" // 84 samples
    // params.lookup_tables = "Lookup_tables"


    // Channels
    bam_file_channel = Channel.fromPath( params.bam_file, checkIfExists: true )
    reference_genome_channel = Channel.fromPath( params.reference_genome, checkIfExists: true )
    downsampled_lookup_tables = Channel.fromPath( "${params.lookup_tables}/lk_downsampled_*.csv", checkIfExists: true ).collect()

    bam_and_fa = bam_file_channel.combine(reference_genome_channel)

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

    PREFIX_FILENAME(bam_and_fa, params.prefix_filename, params.seed)

    SUBSAMPLE_BAM(PREFIX_FILENAME.out, params.depth_range)

    FREEBAYES(SUBSAMPLE_BAM.out)

    PAIRWISE_TABLE(FREEBAYES.out, params.single_end, params.window_size)

    RECOM_RATE_ESTIMATOR(PAIRWISE_TABLE.out, downsampled_lookup_tables, params.recom_tract_len, params.depth_range, params.n_bootstrap_samples, params.ldpop_rho_range)

    // FINAL_RESULTS_PLOT(RECOM_RATE_ESTIMATOR.out)

}