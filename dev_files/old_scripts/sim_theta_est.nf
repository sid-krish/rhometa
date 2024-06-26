#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


def helpMessage() {

    log.info"""
    Description:

    Usage:
    nextflow run theta_est.nf --bam_file in.bam --reference_genome ref.fa [options]

    Help:
    nextflow run theta_est.nf --help

    Required:
    --bam_file [*.bam], Query name sorted bam file. Multi bam support via glob input e.g. "*.bam", quotes but be included for glob. Use with one fasta file only
    --reference_genome [*.fa],  Single/Multi genome fasta file

    Options:
    --prefix_filename [str], prefix string to output filenames to help distinguish runs
    --output_dir [str], default:[Theta_Est_Output], Directory to save results in
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
    prefix_filename.py ${bam_and_fa[0]} ${prefix_fn}
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
            path(fasta)
    
    output:
        tuple val(prefix_filename),
            path("filtered.bam"),
            path(fasta)

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
        path(fasta)

    output:
    tuple val(prefix_filename), 
        path('Aligned_sorted.bam'), 
        path(fasta)

    """
    samtools sort --threads $task.cpus -o Aligned_sorted.bam ${bam}
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

    // maxForks 1

    input:
        tuple val(prefix_filename),
            path(bam),
            path(fasta)

    output:
        tuple val(prefix_filename),
            path(bam),
            path(fasta),
            path("freebayes_filt.vcf")
        
        path 'freebayes_raw.vcf'

    script:
    """
    samtools faidx ${fasta}
    samtools index -@ $task.cpus ${bam}

    # call variants with freebayes
    freebayes -f ${fasta} -p 1 ${bam} > freebayes_raw.vcf

    # keep only SNPs and remove low quality calls
    #bcftools filter --threads ${task.cpus} \
    #    -i 'TYPE="snp" && QUAL>=${params.snp_qual} && FORMAT/DP>=${params.min_snp_depth} && FORMAT/RO>=2 && FORMAT/AO>=2' freebayes_raw.vcf > freebayes_filt.vcf

    bcftools filter --threads ${task.cpus} \
        -i 'TYPE="snp"' freebayes_raw.vcf > freebayes_filt.vcf
    """
}


process THETA_ESTIMATE {
    publishDir params.output_dir, mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    input:
        tuple val(prefix_filename),
            path(bam),
            path(fasta),
            path(vcf)

    output:
        // path "Aligned_sorted.pileup"
        path "theta_estimates.png"
        path "depth_distribution.png"
        path "Theta_estimate_stats.csv"
        
    script:
    """
    samtools mpileup ${bam} > Aligned_sorted.pileup
    genome_size=\$(samtools view -H ${bam} | grep "@SQ" | awk '{ print \$3 }' | cut -c 4-)

    theta_est.py \$genome_size Aligned_sorted.pileup ${vcf}
    """
}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // For each process there is a output of tuple with the necessary files/values to move forward until they are no longer need

    // Params
    params.help = false
    params.prefix_filename = "none"
    // VCF filter settings
    params.snp_qual = 20 // Minimum phred-scaled quality score to filter vcf by
    params.min_snp_depth = 10 // Minimum read depth to filter vcf by

    params.output_dir = 'Theta_Est_Output'
    params.bam_file = 'none'
    params.reference_genome = 'none'
    
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
    PREFIX_FILENAME(bam_and_fa, params.prefix_filename)

    FILTER_BAM(PREFIX_FILENAME.out)

    SORT_BAM(FILTER_BAM.out[0])

    FREEBAYES(SORT_BAM.out)

    // freebayes returns two channels, we just need the first
    THETA_ESTIMATE(FREEBAYES.out[0])

}