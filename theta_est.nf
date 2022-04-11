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
    --min_snp_depth [int], Phred-scaled quality score to filter vcf by
    --output_dir [str], default:[Theta_Est_Output], Directory to save results in

    """.stripIndent()

}


process PREFIX_FILENAME {
    
    // maxForks 1

    // echo true

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
    prefix_filename.py ${bam} ${prefix_fn} 
    """
}


process FREEBAYES {
    /**
      * Call variants using FreeBayes and filter the result VCF for quality follow FreeBayes reccomendations.
      * - Only snps
      * - Variant quality minimum
      * - Depths: "DP, AO and RO" minimums
      **/

    publishDir params.output_dir, mode: 'copy', pattern: '*.vcf', saveAs: {filename -> "freebayes/${filename}"}

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
    samtools sort --threads $task.cpus ${bam} -o Aligned.csorted.bam
    samtools index -@ $task.cpus Aligned.csorted.bam

    # call variants with freebayes
    freebayes -f ${fasta} -p 1 Aligned.csorted.bam > freebayes_raw.vcf

    # keep only SNPs and remove low quality calls
    bcftools filter --threads ${task.cpus} \
        -i 'TYPE="snp" && QUAL>=${params.min_snp_depth} && FORMAT/DP>20 && FORMAT/RO>=2 && FORMAT/AO>=2' freebayes_raw.vcf > freebayes_filt.vcf
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
    samtools sort ${bam} > Aligned_sorted.bam
    samtools mpileup Aligned_sorted.bam > Aligned_sorted.pileup
    genome_size=\$(samtools view -H Aligned_sorted.bam | grep "@SQ" | awk '{ print \$3 }' | cut -c 4-)

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
    params.min_snp_depth = 20

    params.output_dir = 'Theta_Est_Output'
    params.bam_file = 'none'
    params.reference_genome = 'none'

    // Channels
    bam_file_channel = Channel.fromPath( params.bam_file, checkIfExists: true )
    reference_genome_channel = Channel.fromPath( params.reference_genome, checkIfExists: true )

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
    PREFIX_FILENAME(bam_and_fa, params.prefix_filename)

    FREEBAYES(PREFIX_FILENAME.out)

    // freebayes returns two channels, we just need the first
    THETA_ESTIMATE(FREEBAYES.out[0])

}