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
    --bam_file [*.bam], Query name sorted bam file
    --reference_genome [*.fa],  Single genome fasta file

    Options:
    --prepend_filename [str], Prepend string to output filenames to help distinguish runs

    """.stripIndent()

}

process LOFREQ{
    // publishDir "Theta_Est_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

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


process THETA_ESTIMATE {
    publishDir "Theta_Est_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

    maxForks 1

    input:
        path bam
        path vcf
        val prepend_filename

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

    m_theta.py \$genome_size Aligned_sorted.pileup lofreqOut.vcf
    """
}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // Params
    params.help = false
    params.prepend_filename = ""
    params.bam_file = 'none'
    params.reference_genome = 'none'

    // Channels
    bam_file_channel = Channel.fromPath( params.bam_file )
    reference_genome_channel = Channel.fromPath( params.reference_genome )

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
    LOFREQ(reference_genome_channel, bam_file_channel, params.prepend_filename)

    THETA_ESTIMATE(bam_file_channel, LOFREQ.out.lofreqOut_vcf, params.prepend_filename)

}