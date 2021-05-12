#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

process LOFREQ{
    publishDir "Output", mode: "copy"

    maxForks 1

    cpus 4
    memory '1 GB'

    executor 'local'
    time '30m'
    scratch true
    // queue 'i3q'

    input:
        path reference_genome_file
        path bam_file
        
    output:
        path "lofreqOut.vcf", emit: lofreqOut_vcf

    script:
    """
    samtools faidx ${reference_genome_file}
    samtools sort --threads 4 ${bam_file} -o chrom_sorted.bam
    samtools index -@ 4 chrom_sorted.bam

    #lofreq call -f ${reference_genome_file} -o lofreqOut.vcf chrom_sorted.bam
    lofreq call-parallel --pp-threads 4 --no-default-filter -f ${reference_genome_file} -o lofreqOut.vcf chrom_sorted.bam
    """
}


process WATTERSON_ESTIMATE {
    publishDir "Output", mode: "copy"
    
    maxForks 1

    cpus 1
    memory '512 MB'

    executor 'local'
    time '10m'
    scratch true
    // queue 'i3q'

    input:
        path lofreqOut_vcf
        // val sample_size
        // val genome_size

    output:
        path "theta.txt", emit: theta_txt

    script:
    """
    watterson_estimate.py lofreqOut.vcf ${params.genome_size} ${params.sample_size} > theta.txt
    """
}


workflow {
    params.bam_file = ""
    params.reference_genome = ""

    bam_file_channel = Channel.fromPath( params.bam_file )
    reference_genome_channel = Channel.fromPath( params.reference_genome )

    params.sample_size = 180 // Needs to programatically determied
    params.genome_size = 810781 // Needs to programatically determied
    params.recom_tract_len = 500 // Needs to programatically determied

    params.ldpop_rho_range = "101,100"

    LOFREQ(reference_genome_channel, bam_file_channel)

    WATTERSON_ESTIMATE(LOFREQ.out.lofreqOut_vcf)
}