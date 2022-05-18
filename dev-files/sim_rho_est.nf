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
    --lookup_tables [dir], default:[Lookup_tables], Folder containing downsampled lookup tables generated by lookup_table_gen.nf

    Options:
    --ldpop_rho_range [int,int], default:[101,100], The range of rho values used to generate lookup tables
    --recom_tract_len [int], default:[500], Recombination tract length to use
    --window_size [int], default:[500], Window size for variant pairs. For single end this is the read size, for paired end this is the max fragment length
    --single_end, Toggle used for single end read bams
    --depth_range [int,int], default:[3,100], Minimum and maximum depth downsampled lookup tables available. Minimum should be no less than 3
    --prefix_filename [str], prefix string to output filenames to help distinguish runs
    --output_dir [str], default:[Rho_Est_Output], Directory to save results in
    --snp_qual [int], default:[20], Minimum phred-scaled quality score to filter vcf by
    --min_snp_depth [int], default:[10], Minimum read depth to filter vcf by

    """.stripIndent()

}

process PREFIX_FILENAME {

    // maxForks 1

    // echo true

    input:
        tuple val(sample_id),
            path(bam_and_fa)
        
        val(prefix_fn)

    output:
        tuple stdout,
            path("${bam_and_fa[0]}"),
            path("${bam_and_fa[1]}")

    script:
    """
    filename_prefix.py ${bam_and_fa[0]} ${prefix_fn}
    """
}

process GET_SEED_VAL{
    // publishDir "Recom_Est_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    // echo true

    input:
        tuple val(prefix_filename),
            path(bam),
            path(fasta)
    
    output:
        tuple val(prefix_filename),
            path(bam),
            path(fasta),
            stdout

    script:
    """
    get_seed_val.py ${prefix_filename}
    """
}

process FILTER_BAM {
    /**
      * Filter input BAM files to assure that reported alignments will be of acceptable quality.
      * - filters on mapping quality and alignment-score relative to read-length.
      * - publishes a flagstat report of before and after.
      **/

    publishDir params.output_dir, mode: 'copy', pattern: 'flagstat.*.txt', saveAs: {filename -> "filter_bam/${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            path(bam),
            path(fasta),
            val(seed)
    
    output:
        tuple val(prefix_filename),
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
    tuple val(prefix_filename), 
        path(bam), 
        path(fasta), 
        val(seed)

    output:
    tuple val(prefix_filename), 
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
    tuple val(prefix_filename), 
        path('Aligned_sorted.bam'), 
        path(fasta), 
        val(seed)

    output:
    tuple val(prefix_filename), 
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
    tuple val(prefix_filename), 
        path('Aligned_sorted.bam'), 
        path('Aligned_sorted.pileup'), 
        path(fasta), 
        val(seed)

    val depth_range

    output:
    tuple val(prefix_filename), 
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
    publishDir params.output_dir, mode: 'copy', pattern: '*.vcf', saveAs: {filename -> "freebayes/${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            path(bam),
            path(fasta),
            val(seed)

    output:
        tuple val(prefix_filename),
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

    # keep only SNPs and remove low quality calls
    #bcftools filter --threads ${task.cpus} \
    #    -i 'TYPE="snp" && QUAL>=${params.snp_qual} && FORMAT/DP>=${params.min_snp_depth} && FORMAT/RO>=2 && FORMAT/AO>=2' freebayes_raw.vcf > freebayes_filt.vcf

    bcftools filter --threads ${task.cpus} \
        -i 'TYPE="snp"' freebayes_raw.vcf > freebayes_filt.vcf
    """
}


process PAIRWISE_TABLE_SINGLE_END{
    /**
      * Create pair-wise table for final stage of rhometa analysis.
      **/
    // publishDir params.output_dir, mode: 'copy', saveAs: {filename -> "pairwise_table/${prefix_filename}${filename}"}

    // echo true

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
    samtools sort -n --threads $task.cpus -o qsorted.bam ${bam}
    pairwise_table_single_end.py qsorted.bam ${vcf_file} $task.cpus ${window_size}
    """
}


process PAIRWISE_TABLE_PAIRED_END{
    /**
      * Create pair-wise table for final stage of rhometa analysis.
      **/
    // publishDir params.output_dir, mode: 'copy', saveAs: {filename -> "pairwise_table/${prefix_filename}${filename}"}

    // echo true

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
    samtools sort -n --threads $task.cpus -o qsorted.bam ${bam}
    pairwise_table_paired_end.py qsorted.bam ${vcf_file} $task.cpus ${window_size}
    """
}


process RECOM_RATE_ESTIMATOR {
    /**
      * Maximum likelihood estimation of recombination rate.
      **/
    publishDir params.output_dir, mode: 'copy', saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            path(pairwise_table_pkl),
            val(seed)

        path downsampled_lookup_tables
        val recom_tract_len
        val depth_range
        val ldpop_rho_range


    output:
        tuple val(prefix_filename),
            path("log_likelihood_sums.csv"),
            path("rho_estimate.csv")

    script:
    """
    main_weighted.py ${recom_tract_len} ${depth_range} ${ldpop_rho_range} ${pairwise_table_pkl} $task.cpus
    """

}


process RESULTS_PLOT {
    /**
      * Plot of maximum likelihood search over requested rho range
      **/
    publishDir params.output_dir, mode: 'copy', saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
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
    params.prefix_filename = "none"
    params.recom_tract_len = 1000
    params.ldpop_rho_range = "201,2"
    params.window_size = 1000 // For single end this is the read size, for paired end this is the max insert length (1000bp is a practical upper limit)
    params.single_end = false
    params.depth_range = "3,200" // min_depth, max_depth
    // VCF filter settings
    // params.snp_qual = 20 // Minimum phred-scaled quality score to filter vcf by
    // params.min_snp_depth = 10 // Minimum read depth to filter vcf by

    params.output_dir = 'Rho_Est_Output'
    params.lookup_tables = "/shared/homes/11849395/Lookup_tables/Lookup_tables_0-2/Lookup_tables"

    // Channels
    downsampled_lookup_tables = Channel.fromPath( "${params.lookup_tables}/lk_downsampled_*.csv", checkIfExists: true ).collect()

    bam_and_fa = Channel.fromFilePairs('./Sim_Gen_Output/*.{bam,fa}', checkIfExists: true)

    // Input verification
    if (params.help) {
        // Show help message from helpMessage() function
        // params.help = false by default
        helpMessage()
        exit 0
    }

    // Process execution

    PREFIX_FILENAME(bam_and_fa, 
                    params.prefix_filename)

    GET_SEED_VAL(PREFIX_FILENAME.out)
    
    FILTER_BAM(GET_SEED_VAL.out)
    
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
                            params.recom_tract_len, 
                            params.depth_range,
                            params.ldpop_rho_range)
    }
    
    else if (params.single_end == false) {
        PAIRWISE_TABLE_PAIRED_END(FREEBAYES.out[0], 
                    params.single_end, 
                    params.window_size)

        RECOM_RATE_ESTIMATOR(PAIRWISE_TABLE_PAIRED_END.out, 
                            downsampled_lookup_tables, 
                            params.recom_tract_len, 
                            params.depth_range,
                            params.ldpop_rho_range)
    }

    RESULTS_PLOT(RECOM_RATE_ESTIMATOR.out)

}
