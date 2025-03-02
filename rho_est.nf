#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


def helpMessage() {

    log.info"""
    Usage:
    nextflow run rho_est.nf --bam in.bam --fa ref.fa [options]

    Help:
    nextflow run rho_est.nf --help

    Required:
    --bam [*.bam], Multi bam support via glob input e.g. "*.bam", quotes but be included for glob. Use with one fasta file only
    --fa [*.fa],  Single/Multi genome fasta file
    --lookup_tables [dir], default:[Lookup_tables], Folder containing downsampled lookup tables generated by lookup_table_gen.nf

    Options:
    --lookup_grid [int,int], default:[101,100], The range of rho values used to generate lookup tables
    --tract_len [int], default:[1000], Recombination tract length to use
    --window_size [int], default:[1000], Window size for variant pairs. For single end this is the read size, for paired end this is the max fragment length
    --single_end, Toggle used for single end read bams
    --depth_range [int,int], default:[3,85], Minimum and maximum depth downsampled lookup tables available. Minimum should be no less than 3
    --seed [int] , default:[123], Seed value for samtools subsamping. 
                                The seed value will be displayed at the start of the output file names.
                                A list of seed values can be used, a run will be performed for each seed.
    --filename_prefix [str], prefix string to output filenames to help distinguish runs
    --output_dir [str], default:[Rho_Est_Output], Directory to save results in
    --snp_qual [int], default:[20], Minimum phred-scaled quality score to filter vcf by
    --min_snp_depth [int], default:[10], Minimum read depth to filter vcf by
    """.stripIndent()

}


process FILENAME_PREFIX {

    input:
        tuple path(bam),
            path(fasta)
        
        val(prefix_fn)
        each seed

    output:
        tuple stdout,
            path(bam),
            path(fasta),
            val(seed)

    script:
    """
    filename_prefix_seeded.py ${bam} ${prefix_fn} ${seed}
    """
}


process SORT_BAM {
    /**
      * Simply sort the BAM in coordinate order.
      **/
    
    input:
    tuple val(filename_prefix), 
        path(bam), 
        path(fasta), 
        val(seed)

    output:
    tuple val(filename_prefix), 
        path('aligned_sorted.bam'), 
        path(fasta), 
        val(seed)

    """
    samtools sort -@ ${task.cpus} -o aligned_sorted.bam ${bam}
    """
}


process MAKE_PILEUP {
    /**
      * Create a mpileup file, a pre-req for down-sampling
      **/

    input:
    tuple val(filename_prefix), 
        path('aligned_sorted.bam'), 
        path(fasta), 
        val(seed)

    output:
    tuple val(filename_prefix), 
        path('aligned_sorted.bam'), 
        path('aligned_sorted.pileup'), 
        path(fasta), 
        val(seed)

    """
    samtools mpileup -o aligned_sorted.pileup aligned_sorted.bam
    """
}


process SUBSAMPLE {
    /**
      * Down-sample the BAM file to a given maximum depth.
      **/

    publishDir params.output_dir, mode: 'copy', pattern: '*.txt', saveAs: {filename -> "subsample/${filename_prefix}${filename}"}
    
    input:
    tuple val(filename_prefix), 
        path('aligned_sorted.bam'), 
        path('aligned_sorted.pileup'), 
        path(fasta), 
        val(seed)

    val depth_range

    output:
    tuple val(filename_prefix), 
        path("subsampled.bam"), 
        path(fasta), 
        val(seed)

    path("subsample_fraction.txt")

    """
    subsample_bam_seeded.py aligned_sorted.pileup ${depth_range} aligned_sorted.bam ${seed} > subsample_fraction.txt
    """
}


process FREEBAYES {
    publishDir params.output_dir, mode: 'copy', pattern: '*.vcf', saveAs: {filename -> "freebayes/${filename_prefix}${filename}"}

    input:
        tuple val(filename_prefix),
            path(bam),
            path(fasta),
            val(seed)

    output:
        tuple val(filename_prefix),
            path(bam),
            path(fasta),
            path("freebayes_raw.vcf"),
            val(seed)

    script:
    """
    samtools faidx ${fasta}
    samtools sort -@ $task.cpus ${bam} -o Aligned.csorted.bam
    samtools index -@ $task.cpus Aligned.csorted.bam

    # call variants with freebayes
    freebayes -f ${fasta} -p 1 Aligned.csorted.bam > freebayes_raw.vcf
    """
}


process VCF_FILTER {
    /**
      * Filter the VCF.
      * - SNPS Only
      * - Minimum variant quality
      * - Depths: "DP, AO and RO" minimums
      * - Outlier depth cutoff
      **/
    publishDir params.output_dir, mode: 'copy', pattern: '*.vcf', saveAs: {filename -> "freebayes/${filename_prefix}${filename}"}

    // maxForks 1

    input:
        tuple val(filename_prefix),
            path(bam),
            path(fasta),
            path("freebayes_raw.vcf"),
            val(seed)

        val snp_qual
        val min_snp_depth
        val top_depth_cutoff_percentage

    output:
        tuple val(filename_prefix),
            path(bam),
            path(fasta),
            path("freebayes_filt.vcf"),
            val(seed)

    script:
    """
    vcf_filter.py ${task.cpus} freebayes_raw.vcf ${snp_qual} ${min_snp_depth} ${top_depth_cutoff_percentage}
    """
}


process PAIRWISE_TABLE_SINGLE_END{
    /**
      * Create pair-wise table for final stage of rhometa analysis.
      **/
    // publishDir params.output_dir, mode: 'copy', pattern: '*.pkl', saveAs: {filename -> "pairwise_table/${filename_prefix}${filename}"}

    // echo true

    input:
        tuple val(filename_prefix),
            path(bam),
            path(fasta),
            path(vcf_file),
            val(seed)

        val single_end
        val window_size

    output:
        tuple val(filename_prefix),
            path("pairwise_table.pkl"),
            val(seed),
            path(vcf_file),
            env(genome_size)

    script:
    // -n Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates.
    """
    samtools sort -n -@ $task.cpus -o qsorted.bam ${bam}
    pairwise_table_single_end.py qsorted.bam ${vcf_file} $task.cpus ${window_size}

    #New solution works for bams aligned to single/multiple sequences by summing the lengths of all @SQ records.
    genome_size=\$(samtools view -H ${bam} |  awk '/^@SQ/ {l+=substr(\$3,4)}END{print l}')
    """
}


process PAIRWISE_TABLE_PAIRED_END{
    /**
      * Create pair-wise table for final stage of rhometa analysis.
      **/
    // publishDir params.output_dir, mode: 'copy', pattern: '*.pkl', saveAs: {filename -> "pairwise_table/${filename_prefix}${filename}"}

    // echo true

    input:
        tuple val(filename_prefix),
            path(bam),
            path(fasta),
            path(vcf_file),
            val(seed)

        val single_end
        val window_size

    output:
        tuple val(filename_prefix),
            path("pairwise_table.pkl"),
            val(seed),
            path(vcf_file),
            env(genome_size)

    script:
    // -n Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates.
    """
    samtools sort -n -@ $task.cpus -o qsorted.bam ${bam}
    pairwise_table_paired_end.py qsorted.bam ${vcf_file} $task.cpus ${window_size}

    #New solution works for bams aligned to single/multiple sequences by summing the lengths of all @SQ records.
    genome_size=\$(samtools view -H ${bam} |  awk '/^@SQ/ {l+=substr(\$3,4)}END{print l}')
    """
}


process RHO_ESTIMATE {
    /**
      * Composite likelihood estimation of recombination rate.
      **/

    // debug true

    publishDir params.output_dir, mode: 'copy', saveAs: {filename -> "rho_estimate/${filename_prefix}${filename}"}

    input:
        tuple val(filename_prefix),
            path(pairwise_table_pkl),
            val(seed),
            path(vcf_file),
            val(genome_size)

        path downsampled_lookup_tables
        val tract_len
        val depth_range
        val lookup_grid


    output:
        tuple val(filename_prefix),
            path("log_likelihood_sums.csv"),
            path("rho_estimate.csv")

    script:
    """
    subst_probability=\$(subst_probability.py $genome_size ${vcf_file})
    main_weighted.py ${tract_len} ${depth_range} ${lookup_grid} ${pairwise_table_pkl} $task.cpus $genome_size \$subst_probability
    """
}


process RESULTS_PLOT {
    /**
      * Plot of maximum likelihood search over requested rho range
      **/

    publishDir params.output_dir, mode: 'copy', saveAs: {filename -> "results_plot/${filename_prefix}${filename}"}

    input:
        tuple val(filename_prefix),
            path("log_likelihood_sums.csv"),
            path("rho_estimate.csv")

    output:
        path "results_plot.png"

    script:
    """
    results_plot.py log_likelihood_sums.csv
    """
}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // For each process there is a output of tuple with the necessary files/values to move forward until they are no longer need

    // Params
    params.help = false
    params.seed = [123] // used for samtools subsamping and final bootstrap algorithm
    params.filename_prefix = "none"
    params.tract_len = 1000
    params.lookup_grid = "101,100" // The range of rho values used to generate lookup tables
    params.window_size = 1000 // For single end this is the read size, for paired end this is the max insert length (1000bp is a practical upper limit)
    params.single_end = false
    params.depth_range = "3,85" // min_depth, max_depth
    
    // VCF filter settings
    params.snp_qual = 20 // Minimum phred-scaled quality score to filter vcf by
    params.min_snp_depth = 10 // Minimum read depth to filter vcf by
    params.top_depth_cutoff_percentage = 5 // Top n percent of depth to cutoff from the vcf file

    params.output_dir = 'Rho_Est_Output'
    params.bam = 'none'
    params.fa = 'none'
    params.lookup_tables = "Lookup_tables"


    // Channels
    bam_channel = Channel.fromPath( params.bam, checkIfExists: true )
    fa_channel = Channel.fromPath( params.fa, checkIfExists: true )
    downsampled_lookup_tables = Channel.fromPath( "${params.lookup_tables}/lk_downsampled_*.csv", checkIfExists: true ).collect()

    bam_and_fa = bam_channel.combine(fa_channel)

    // Input verification
    if (params.help) {
        // Show help message from helpMessage() function
        // params.help = false by default
        helpMessage()
        exit 0
    }

    if (params.fa == 'none') {
        println "No input .fa specified. Use --fa [.fa]"
        exit 1
    }

    if (params.bam == 'none') {
        println "No input .bam specified. Use --bam [.bam]"
        exit 1
    }

    // Process execution
    FILENAME_PREFIX(bam_and_fa, 
                    params.filename_prefix, 
                    params.seed)
    
    SORT_BAM(FILENAME_PREFIX.out)
    
    MAKE_PILEUP(SORT_BAM.out)
    
    SUBSAMPLE(MAKE_PILEUP.out, 
              params.depth_range)
    
    FREEBAYES(SUBSAMPLE.out[0])

    VCF_FILTER(FREEBAYES.out, params.snp_qual, params.min_snp_depth, params.top_depth_cutoff_percentage)

    if (params.single_end == true) {
        PAIRWISE_TABLE_SINGLE_END(VCF_FILTER.out, 
                    params.single_end, 
                    params.window_size)

        RHO_ESTIMATE(PAIRWISE_TABLE_SINGLE_END.out, 
                            downsampled_lookup_tables, 
                            params.tract_len, 
                            params.depth_range,
                            params.lookup_grid)
    }
    
    else if (params.single_end == false) {
        PAIRWISE_TABLE_PAIRED_END(VCF_FILTER.out, 
                    params.single_end, 
                    params.window_size)

        RHO_ESTIMATE(PAIRWISE_TABLE_PAIRED_END.out, 
                            downsampled_lookup_tables, 
                            params.tract_len, 
                            params.depth_range,
                            params.lookup_grid)
    }

    // RESULTS_PLOT(RHO_ESTIMATE.out)

}
