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


process LOFREQ{
    // publishDir "Theta_Est_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

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
    // publishDir "Theta_Est_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    input:
        tuple val(prefix_filename),
            path(bam),
            path(fasta)

    output:
        tuple val(prefix_filename),
            path(bam),
            path(fasta),
            path("freeBayesOut.vcf")

    script:
    """
    samtools faidx ${fasta}
    samtools sort --threads $task.cpus ${bam} -o Aligned.csorted.bam
    samtools index -@ $task.cpus Aligned.csorted.bam

    # only keep SNP type entries
    freebayes -f ${fasta} -p 1 Aligned.csorted.bam | grep -e '^#' -e 'TYPE=snp' > freeBayesOut.vcf
    """
}


process THETA_ESTIMATE {
    publishDir "Theta_Est_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

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

    m_theta_final.py \$genome_size Aligned_sorted.pileup ${vcf}
    """
}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // For each process there is a output of tuple with the necessary files/values to move forward until they are no longer need

    // Params
    params.help = false
    params.prefix_filename = "none"

    params.bam_file = 'none'
    params.reference_genome = 'none'

    // Channels
    bam_file_channel = Channel.fromPath( params.bam_file )
    reference_genome_channel = Channel.fromPath( params.reference_genome )

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

    THETA_ESTIMATE(FREEBAYES.out)

}