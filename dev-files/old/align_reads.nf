#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process BWA_MEM_PAIRED_END {
    publishDir "Align_Reads_Output", mode: "copy", saveAs: {filename -> "${filename}"}

    // echo true

    input:
        tuple val(sample_id),
        path(fastqs),
        path(reference_fa)

    output:
        path("${sample_id}_Aligned_Sorted.bam")

    script:
    """
    #echo ${sample_id}
    #echo ${fastqs[0]}
    #echo ${fastqs[1]}
    #echo ${reference_fa}

    bwa index ${reference_fa}

    #Paired end
    bwa mem -t $task.cpus ${reference_fa} ${fastqs[0]} ${fastqs[1]} > ${sample_id}_Aligned.sam

    samtools view -bS ${sample_id}_Aligned.sam > ${sample_id}_Aligned.bam

    samtools sort --threads $task.cpus ${sample_id}_Aligned.bam -o ${sample_id}_Aligned_Sorted.bam
    """
}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // For each process there is a output of tuple with the necessary files/values to move forward until they are no longer need

    // Params
    params.fastqs = "*_{1,2}.fastq"
    params.reference_genome = 'none'

    // Channels
    fastqs_channel = Channel.fromFilePairs( params.fastqs)
    reference_genome_channel = Channel.fromPath( params.reference_genome )

    combined_inputs = fastqs_channel.combine(reference_genome_channel)

    // Process execution
    BWA_MEM_PAIRED_END(combined_inputs)

}