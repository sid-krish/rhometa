#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


def helpMessage() {

    log.info"""
    Usage:
    nextflow run theta_est.nf --bam in.bam --fa ref.fa [options]

    Help:
    nextflow run theta_est.nf --help

    Required:
    --bam [*.bam], Query name sorted bam file. Multi bam support via glob input e.g. "*.bam", quotes must be included for glob. Use with one fasta file only
    --fa [*.fa],  Single/Multi genome fasta file

    Options:
    --filename_prefix [str], prefix string to output filenames to help distinguish runs
    --output_dir [str], default:[Theta_Est_Output], Directory to save results in
    --snp_qual [int], default:[20], Minimum phred-scaled quality score to filter vcf by
    --min_snp_depth [int], default:[10], Minimum read depth to filter vcf by
    """.stripIndent()
}


process PREFIX_FILENAME {

    // maxForks 1

    // echo true

    label 'LOCAL_EXEC'

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
    #!/usr/bin/env python

    bam = "${bam_and_fa[0]}"
    prepend_file = "${prefix_fn}"

    # bam = "a.s_bs_ds.bam"
    # prepend_file = "none"

    if prepend_file == "none":
        file_name = bam.rsplit(".",1)[0]
        print(f"{file_name}_", end = '')

    else:
        print(prepend_file, end = '')
    """
}


process SORT_BAM {
    /**
      * Simply sort the BAM in coordinate order.
      **/
    
    input:
    tuple val(filename_prefix), 
        path(bam), 
        path(fasta)

    output:
    tuple val(filename_prefix), 
        path('Aligned_sorted.bam'), 
        path(fasta)

    """
    samtools sort -@ $task.cpus -o Aligned_sorted.bam ${bam}
    """
}

process FREEBAYES {
    publishDir params.output_dir, mode: 'copy', pattern: '*.vcf', saveAs: {filename -> "freebayes/${filename_prefix}${filename}"}

    // maxForks 1

    input:
        tuple val(filename_prefix),
            path(bam),
            path(fasta)

    output:
        tuple val(filename_prefix),
            path(bam),
            path(fasta),
            path("freebayes_raw.vcf")

    script:
    """
    samtools faidx ${fasta}
    samtools index -@ $task.cpus ${bam}

    # call variants with freebayes
    # --haplotype-length=-1: Allow haplotype calls with contiguous embedded matches of up to this length. Set N=-1 to disable clumping.
    freebayes -f ${fasta} -p 1 ${bam} --haplotype-length=-1 --min-alternate-count 2 > freebayes_raw.vcf
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
            path("freebayes_raw.vcf")

        val snp_qual
        val min_snp_depth
        val top_depth_cutoff_percentage

    output:
        tuple val(filename_prefix),
            path(bam),
            path(fasta),
            path(params.VCF_FILTERED_FILE)

    script:
    """
    vcf_filter.py ${task.cpus} freebayes_raw.vcf ${snp_qual} ${min_snp_depth} ${top_depth_cutoff_percentage} ${params.VCF_FILTERED_FILE}
    """
}

process THETA_ESTIMATE {
    publishDir params.output_dir, mode: "copy", saveAs: {filename -> "theta_estimate/${filename_prefix}${filename}"}

    // maxForks 1

    debug true

    input:
        tuple val(filename_prefix),
            path(bam),
            path(fasta),
            path(vcf)

    output:
        // path "Aligned_sorted.pileup"
        path "Theta_estimate_stats.csv"
        
    script:
    """
    samtools mpileup ${bam} > Aligned_sorted.pileup

    #old solution only worked for bams aligned to single sequence
    #genome_size=\$(samtools view -H ${bam} | grep "@SQ" | awk '{ print \$3 }' | cut -c 4-)

    #New solution works for bams aligned to single/multiple sequences by summing the lengths of all @SQ records.
    genome_size=\$(samtools view -H ${bam} |  awk '/^@SQ/ {l+=substr(\$3,4)}END{print l}')

    theta_est.py \$genome_size Aligned_sorted.pileup ${vcf}
    """
}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // For each process there is a output of tuple with the necessary files/values to move forward until they are no longer need

    // Params
    params.help = false
    params.filename_prefix = "none"
    // VCF filter settings
    params.snp_qual = 20 // Minimum phred-scaled quality score to filter vcf by
    params.min_snp_depth = 10 // Minimum read depth to filter vcf by
    params.top_depth_cutoff_percentage = 5 // Top n percent of depth to cutoff from the vcf file

    params.output_dir = 'Theta_Est_Output'
    params.bam = 'none'
    params.fa = 'none'

    // Output file names
    params.VCF_FILTERED_FILE = 'freebayes_filt.vcf'
    
    // Channels
    bam_and_fa = Channel.fromFilePairs('./Sim_Gen_Output/*.{bam,fa}', checkIfExists: true)

    // Input verification
    if (params.help) {
        // Show help message from helpMessage() function
        // params.help = false by default
        helpMessage()
        exit 0
    }

    // Process execution
    PREFIX_FILENAME(bam_and_fa, params.filename_prefix)

    SORT_BAM(PREFIX_FILENAME.out)

    FREEBAYES(SORT_BAM.out)

    VCF_FILTER(FREEBAYES.out, params.snp_qual, params.min_snp_depth, params.top_depth_cutoff_percentage)

    THETA_ESTIMATE(VCF_FILTER.out)

    // THETA_EST_PLOT(VCF_FILTER.out)

}