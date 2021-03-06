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
    """.stripIndent()
    // Testing
    // --snp_qual [int], default:[20], Minimum phred-scaled quality score to filter vcf by
    // --min_snp_depth [int], default:[10], Minimum read depth to filter vcf by

}


process FILENAME_PREFIX {

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
    filename_prefix_seeded.py ${bam} ${prefix_fn} ${seed}
    """
}

process FILTER_BAM {
    /**
      * Filter input BAM files to assure that reported alignments will be of acceptable quality.
      * - filters on mapping quality and alignment-score relative to read-length.
      * - publishes a flagstat report of before and after.
      **/

    publishDir params.output_dir, mode: 'copy', pattern: 'flagstat.*.txt', saveAs: {filename -> "filter_bam/${filename_prefix}${filename}"}

    input:
        tuple val(filename_prefix),
            path(bam),
            path(fasta),
            val(seed)
    
    output:
        tuple val(filename_prefix),
            path("filtered.bam"),
            path(fasta),
            val(seed)

        path 'flagstat.*.txt'

    script:
    """
    samtools view -b -e "mapq>=40 && [AS]/rlen>0.75" $bam > filtered.bam
    samtools flagstat $bam > flagstat.before.txt
    samtools flagstat filtered.bam > flagstat.after.txt
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
        path('Aligned_sorted.bam'), 
        path(fasta), 
        val(seed)

    """
    samtools sort -@${task.cpus} -o Aligned_sorted.bam ${bam}
    """
}

process MAKE_PILEUP {
    /**
      * Create a mpileup file, a pre-req for down-sampling
      **/

    input:
    tuple val(filename_prefix), 
        path('Aligned_sorted.bam'), 
        path(fasta), 
        val(seed)

    output:
    tuple val(filename_prefix), 
        path('Aligned_sorted.bam'), 
        path('Aligned_sorted.pileup'), 
        path(fasta), 
        val(seed)

    """
    samtools mpileup -o Aligned_sorted.pileup Aligned_sorted.bam
    """
}

process SUBSAMPLE {
    /**
      * Down-sample the BAM file to a given maximum depth.
      **/
    
    input:
    tuple val(filename_prefix), 
        path('Aligned_sorted.bam'), 
        path('Aligned_sorted.pileup'), 
        path(fasta), 
        val(seed)

    val depth_range

    output:
    tuple val(filename_prefix), 
        path("subsampled.bam"), 
        path(fasta), 
        val(seed)

    """
    subsample_bam_seeded.py Aligned_sorted.pileup ${depth_range} Aligned_sorted.bam ${seed}
    """
}

process FREEBAYES {
    /**
      * Call variants using FreeBayes and filter the result VCF for quality follow FreeBayes reccomendations.
      * - Only snps
      * - Variant quality minimum
      * - Depths: "DP, AO and RO" minimums
      **/
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
            path("freebayes_filt.vcf"),
            val(seed)
        
        path 'freebayes_raw.vcf'

    script:
    """
    samtools faidx ${fasta}
    samtools sort --threads $task.cpus ${bam} -o Aligned.csorted.bam
    samtools index -@ $task.cpus Aligned.csorted.bam

    # call variants with freebayes
    freebayes -f ${fasta} -p 1 Aligned.csorted.bam > freebayes_raw.vcf

    # keep only SNPs
    bcftools filter --threads ${task.cpus} -i 'TYPE="snp"' freebayes_raw.vcf > freebayes_filt.vcf
    """
    // Testing
    // keep SNPs and remove low quality and low depth calls
    // bcftools filter --threads ${task.cpus} \
    //     -i 'TYPE="snp" && QUAL>=${params.snp_qual} && FORMAT/DP>=${params.min_snp_depth} && FORMAT/RO>=2 && FORMAT/AO>=2' freebayes_raw.vcf > freebayes_filt.vcf
}


process PAIRWISE_TABLE_SINGLE_END{
    /**
      * Create pair-wise table for final stage of rhometa analysis.
      **/
    // publishDir params.output_dir, mode: 'copy', saveAs: {filename -> "pairwise_table/${filename_prefix}${filename}"}

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
            val(seed)

    script:
    // -n Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates.
    """
    samtools sort -n --threads $task.cpus -o qsorted.bam ${bam}
    pairwise_table_single_end.py qsorted.bam ${vcf_file} $task.cpus ${window_size}
    """
}


process PAIRWISE_TABLE_PAIRED_END{
    /**
      * Create pair-wise table for final stage of rhometa analysis.
      **/
    // publishDir params.output_dir, mode: 'copy', saveAs: {filename -> "pairwise_table/${filename_prefix}${filename}"}

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
            val(seed)

    script:
    // -n Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates.
    """
    samtools sort -n --threads $task.cpus -o qsorted.bam ${bam}
    pairwise_table_paired_end.py qsorted.bam ${vcf_file} $task.cpus ${window_size}
    """
}


process RECOM_RATE_ESTIMATOR {
    /**
      * Maximum likelihood estimation of recombination rate.
      **/
    publishDir params.output_dir, mode: 'copy', saveAs: {filename -> "${filename_prefix}${filename}"}

    input:
        tuple val(filename_prefix),
            path(pairwise_table_pkl),
            val(seed)

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
    main_weighted.py ${tract_len} ${depth_range} ${lookup_grid} ${pairwise_table_pkl} $task.cpus
    """

}


process RESULTS_PLOT {
    /**
      * Plot of maximum likelihood search over requested rho range
      **/
    publishDir params.output_dir, mode: 'copy', saveAs: {filename -> "${filename_prefix}${filename}"}

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
    // params.snp_qual = 20 // Minimum phred-scaled quality score to filter vcf by
    // params.min_snp_depth = 10 // Minimum read depth to filter vcf by

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
    
    FILTER_BAM(FILENAME_PREFIX.out)
    
    SORT_BAM(FILTER_BAM.out[0])
    
    MAKE_PILEUP(SORT_BAM.out)
    
    SUBSAMPLE(MAKE_PILEUP.out, 
              params.depth_range)
    
    FREEBAYES(SUBSAMPLE.out) // freebayes returns two channels, we just need the first

    if (params.single_end == true) {
        PAIRWISE_TABLE_SINGLE_END(FREEBAYES.out[0], 
                    params.single_end, 
                    params.window_size)

        RECOM_RATE_ESTIMATOR(PAIRWISE_TABLE_SINGLE_END.out, 
                            downsampled_lookup_tables, 
                            params.tract_len, 
                            params.depth_range,
                            params.lookup_grid)
    }
    
    else if (params.single_end == false) {
        PAIRWISE_TABLE_PAIRED_END(FREEBAYES.out[0], 
                    params.single_end, 
                    params.window_size)

        RECOM_RATE_ESTIMATOR(PAIRWISE_TABLE_PAIRED_END.out, 
                            downsampled_lookup_tables, 
                            params.tract_len, 
                            params.depth_range,
                            params.lookup_grid)
    }

    RESULTS_PLOT(RECOM_RATE_ESTIMATOR.out)

}
