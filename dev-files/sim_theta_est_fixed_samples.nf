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


process GET_SAMPLE_SIZE{
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
    get_sample_size.py ${prefix_filename}
    """
}


process FREEBAYES {
    // publishDir "Theta_Est_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    input:
        tuple val(prefix_filename),
            path(bam),
            path(fasta),
            val(samples)

    output:
        tuple val(prefix_filename),
            path(bam),
            path(fasta),
            path("freeBayesOut.vcf"),
            val(samples)

    script:
    """
    samtools faidx ${fasta}
    samtools sort --threads $task.cpus ${bam} -o Aligned.csorted.bam
    samtools index -@ $task.cpus Aligned.csorted.bam

    freebayes -f ${fasta} -p 1 Aligned.csorted.bam > freeBayesOut.vcf
    """
}


process THETA_ESTIMATE {
    publishDir "Theta_Est_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    // echo true

    input:
        tuple val(prefix_filename),
            path(bam),
            path(fasta),
            path(vcf),
            val(samples)

    output:
        path "theta_est.csv"
        path "num_var_sites.csv"
        
    script:
    """
    samtools sort ${bam} > Aligned_sorted.bam
    genome_size=\$(samtools view -H Aligned_sorted.bam | grep "@SQ" | awk '{ print \$3 }' | cut -c 4-)

    snps=\$(num_var_sites.py ${vcf})
    echo \$snps > num_var_sites.csv
    theta=\$(original_watterson.py \$genome_size \$snps ${samples})
    echo \$theta > theta_est.csv
    """
}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // For each process there is a output of tuple with the necessary files/values to move forward until they are no longer need

    // Params
    params.help = false
    params.prefix_filename = "none"

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

    GET_SAMPLE_SIZE(PREFIX_FILENAME.out)

    FREEBAYES(GET_SAMPLE_SIZE.out)

    THETA_ESTIMATE(FREEBAYES.out)

}