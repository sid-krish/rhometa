#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


def helpMessage() {

    log.info"""
    Usage:
    nextflow run theta_est.nf --bam in.bam --fa ref.fa [options]

    Help:
    nextflow run theta_est.nf --help

    Required:
    --input_csv, csv file with paths to bams and references, "bam" and "reference" need to be the header

    Options:
    --filename_prefix [str], prefix string to output filenames to help distinguish runs
    --output_dir [str], default:[Theta_Est_Output], Directory to save results in
    --snp_qual [int], default:[20], Minimum phred-scaled quality score to filter vcf by
    --min_snp_depth [int], default:[10], Minimum read depth to filter vcf by
    """.stripIndent()


}


process FILENAME_PREFIX {
    
    // maxForks 1

    input:
        tuple path(bam),
            path(fasta)
        
        val(prefix_fn)

    output:
        tuple stdout,
            path(bam),
            path(fasta)

    script:
    """
    filename_prefix.py ${bam} ${prefix_fn} 
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
    freebayes -f ${fasta} -p 1 ${bam} > freebayes_raw.vcf
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
            path("freebayes_filt.vcf")

    script:
    """
    vcf_filter.py ${task.cpus} freebayes_raw.vcf ${snp_qual} ${min_snp_depth} ${top_depth_cutoff_percentage}
    """
}


process THETA_ESTIMATE {
    publishDir params.output_dir, mode: "copy", saveAs: {filename -> "theta_estimate/${filename_prefix}${filename}"}

    // maxForks 1

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


process THETA_EST_PLOT {
    publishDir params.output_dir, mode: "copy", saveAs: {filename -> "theta_estimate_plots/${filename_prefix}${filename}"}

    // maxForks 1

    input:
        tuple val(filename_prefix),
            path(bam),
            path(fasta),
            path(vcf)

    output:
        path "theta_estimates.png"
        path "depth_distribution.png"
        
    script:
    """
    samtools mpileup ${bam} > Aligned_sorted.pileup

    #old solution only worked for bams aligned to single sequence
    #genome_size=\$(samtools view -H ${bam} | grep "@SQ" | awk '{ print \$3 }' | cut -c 4-)

    #New solution works for bams aligned to single/multiple sequences by summing the lengths of all @SQ records.
    genome_size=\$(samtools view -H ${bam} |  awk '/^@SQ/ {l+=substr(\$3,4)}END{print l}')

    theta_plot.py \$genome_size Aligned_sorted.pileup ${vcf}
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

    // Input verification
    if (params.help) {
        // Show help message from helpMessage() function
        // params.help = false by default
        helpMessage()
        exit 0
    }

    params.input_csv = "" // csv file with paths to bams and references, "bam" and "reference" need to the header

    input_list = Channel.fromPath(params.input_csv)
                .splitCsv(header:true)
                .map { row-> tuple(file(row.bam), file(row.reference)) }

    // Process execution
    FILENAME_PREFIX(input_list,
                    params.filename_prefix)

    SORT_BAM(FILENAME_PREFIX.out)

    FREEBAYES(SORT_BAM.out)

    VCF_FILTER(FREEBAYES.out, params.snp_qual, params.min_snp_depth, params.top_depth_cutoff_percentage)

    THETA_ESTIMATE(VCF_FILTER.out)

    // THETA_EST_PLOT(VCF_FILTER.out)

}