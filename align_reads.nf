#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process BWA_MEM_SINGLE_END {
    publishDir "Align_Reads_Output", mode: 'copy', pattern: '*.bam', saveAs: {filename -> "${filename}"}

    // echo true

    input:
        tuple path(fastq),
        path(reference_fa)

    output:
        tuple path("${fastq.getSimpleName()}_${reference_fa.getSimpleName()}.bam"),
            path(reference_fa)

    script:
    """
    #echo ${fastq}
    #echo ${reference_fa}

    bwa index ${reference_fa}

    #Single end
    bwa mem -t $task.cpus ${reference_fa} ${fastq} | samtools sort -@ $task.cpus -o ${fastq.getSimpleName()}_${reference_fa.getSimpleName()}.bam
    """
}


process BWA_MEM_PAIRED_END {
    publishDir "Align_Reads_Output", mode: 'copy', pattern: '*.bam', saveAs: {filename -> "${filename}"}

    // echo true

    input:
        tuple val(sample_id),
        path(fastqs),
        path(reference_fa)

    output:
        tuple path("${sample_id}_${reference_fa.getSimpleName()}.bam"),
            path(reference_fa)

    script:
    """
    #echo ${sample_id}
    #echo ${fastqs[0]}
    #echo ${fastqs[1]}
    #echo ${reference_fa}

    bwa index ${reference_fa}

    #Paired end
    bwa mem -t $task.cpus ${reference_fa} ${fastqs[0]} ${fastqs[1]} | samtools sort -@ $task.cpus -o ${sample_id}_${reference_fa.getSimpleName()}.bam
    """
}


process FILTER_BAM {
    publishDir "Align_Reads_Output", mode: 'copy', pattern: 'flagstat.*.txt', saveAs: {filename -> "filter_bam/${filename}"}
    publishDir "Align_Reads_Output", mode: 'copy', pattern: '*.bam', saveAs: {filename -> "${filename}"}

    input:
        tuple path(bam),
            path(fasta)
    
    output:
        tuple path("${bam.getSimpleName()}_filtered.bam"),
            path(fasta)

        path 'flagstat.*.txt'

    script:
    """
    #Filter and keep mapped only
    samtools view -h -F 4 -e "mapq>=40 && [AS]/rlen>0.50" $bam | samtools sort -@ $task.cpus -o ${bam.getSimpleName()}_filtered.bam

    samtools flagstat $bam > flagstat.before.txt
    samtools flagstat ${bam.getSimpleName()}_filtered.bam > flagstat.after.txt
    """
}


process SAMTOOLS_COVERAGE {
    publishDir "Align_Reads_Output", mode: "copy", pattern: '*.cov', saveAs: {filename -> "coverage/${filename}"}

    input:
        tuple path(bam),
            path(fa)

    output:
        tuple path(bam),
            path(fa)

        path("${bam.getSimpleName()}.cov")

    script:
    """
    samtools coverage ${bam} > ${bam.getSimpleName()}.cov
    """
}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // For each process there is a output of tuple with the necessary files/values to move forward until they are no longer need

    // Params
    params.fq = "none" // e.g. with quotes "*{1,2}.fq" for paired end
    params.fa = "none"
    params.single_end = false

    if (params.fq == 'none') {
        println "No input .fq specified. Use --fq [.fq]"
        exit 1
    }

    if (params.fa == 'none') {
        println "No input .fa specified. Use --fa [.fa]"
        exit 1
    }


    if (params.single_end == true) {
        // Channels
        fastqs_channel = Channel.fromPath(params.fq)
        reference_genome_channel = Channel.fromPath(params.fa)
        combined_inputs = fastqs_channel.combine(reference_genome_channel)

        // Process execution
        BWA_MEM_SINGLE_END(combined_inputs)

        FILTER_BAM(BWA_MEM_SINGLE_END.out)

        SAMTOOLS_COVERAGE(FILTER_BAM.out[0])
    }


    else if (params.single_end == false) {
        // Channels
        fastqs_channel = Channel.fromFilePairs(params.fq)
        reference_genome_channel = Channel.fromPath(params.fa)
        combined_inputs = fastqs_channel.combine(reference_genome_channel)

        // Process execution
        BWA_MEM_PAIRED_END(combined_inputs)

        FILTER_BAM(BWA_MEM_PAIRED_END.out)

        SAMTOOLS_COVERAGE(FILTER_BAM.out[0])
    }

}