#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process FILTER_FASTQ_SINGLE_END {
    // publishDir params.align_out_dir, mode: "copy", pattern: '*.gz', saveAs: {filename -> "fastp_out/${filename}"}
    publishDir params.align_out_dir, mode: "copy", pattern: '*.html', saveAs: {filename -> "fastp_out/${filename}"}
    publishDir params.align_out_dir, mode: "copy", pattern: '*.json', saveAs: {filename -> "fastp_out/${filename}"}

    input:
        tuple path(fastq),
        path(reference_fa)

    output:
        tuple path(reference_fa),
            path("${fastq.getSimpleName()}_fp.fq.gz")

        path("${fastq.getSimpleName()}.json")
        path("${fastq.getSimpleName()}.html")

    script:
    """
    fastp --thread 1 --dedup -i ${fastq} -o ${fastq.getSimpleName()}_fp.fq.gz \
    --json ${fastq.getSimpleName()}.json --html ${fastq.getSimpleName()}.html
    """
}


process FILTER_FASTQ_PAIRED_END {
    publishDir params.align_out_dir, mode: 'copy', pattern: '*.html', saveAs: {filename -> "fastp_out/${filename}"}
    publishDir params.align_out_dir, mode: 'copy', pattern: '*.json', saveAs: {filename -> "fastp_out/${filename}"}

    input:
        tuple val(sample_id),
              path(fastq_1),
              path(fastq_2),
              path(reference_fa)

    output:
        tuple val(sample_id),
              path(reference_fa),
              path("${fastq_1.getSimpleName()}_fp.fq.gz"),
              path("${fastq_2.getSimpleName()}_fp.fq.gz")

        path("${fastq_1.getSimpleName()[0..-2]}.json")
        path("${fastq_1.getSimpleName()[0..-2]}.html")

    script:
    """
    fastp --thread 1 --dedup --in1 ${fastq_1} --in2 ${fastq_2} \
    --out1 ${fastq_1.getSimpleName()}_fp.fq.gz --out2 ${fastq_2.getSimpleName()}_fp.fq.gz \
    --json ${fastq_1.getSimpleName()[0..-2]}.json --html ${fastq_1.getSimpleName()[0..-2]}.html
    """
}


process BWA_MEM_SINGLE_END {
    // publishDir params.align_out_dir, mode: 'copy', pattern: '*.bam', saveAs: {filename -> "${filename}"}

    // echo true

    input:
        tuple path(reference_fa),
        path(fastq)

    output:
        tuple path("${fastq.getSimpleName()}_${reference_fa.getSimpleName()}.bam"),
            path(reference_fa)

    script:
    """
    #echo ${fastq}
    #echo ${reference_fa}

    bwa-mem2 index ${reference_fa}

    #Single end
    bwa-mem2 mem -t $task.cpus ${reference_fa} ${fastq} | samtools sort -@ $task.cpus -o ${fastq.getSimpleName()}_${reference_fa.getSimpleName()}.bam
    """
}


process BWA_MEM_PAIRED_END {
    // publishDir params.align_out_dir, mode: 'copy', pattern: '*.bam', saveAs: {filename -> "${filename}"}

    // echo true

    input:
        tuple val(sample_id),
            path(reference_fa),
            path(fastq_0),
            path(fastq_1)

    output:
        tuple path("${sample_id}_${reference_fa.getSimpleName()}.bam"),
            path(reference_fa)

    script:
    """
    #echo ${sample_id}
    #echo ${fastq_0}
    #echo ${fastq_1}
    #echo ${reference_fa}

    bwa-mem2 index ${reference_fa}

    #Paired end
    bwa-mem2 mem -t $task.cpus ${reference_fa} ${fastq_0} ${fastq_1} | samtools sort -@ $task.cpus -o ${sample_id}_${reference_fa.getSimpleName()}.bam
    """
}


process FILTER_BAM {
    publishDir params.align_out_dir, mode: 'copy', pattern: '*flagstat.*.txt', saveAs: {filename -> "filter_bam/${filename}"}
    publishDir params.align_out_dir, mode: 'copy', pattern: '*.bam', saveAs: {filename -> "${filename}"}

    input:
        tuple path(bam),
            path(fasta)
    
    output:
        tuple path("${bam.getSimpleName()}_filtered.bam"),
            path(fasta)

        path("${bam.getSimpleName()}_flagstat.*.txt")

    script:
    """
    #Filter and keep mapped only
    samtools view -h -F 4 -e "mapq>=40 && [AS]/rlen>0.50" $bam | samtools sort -@ $task.cpus -o ${bam.getSimpleName()}_filtered.bam

    samtools flagstat $bam > ${bam.getSimpleName()}_flagstat.before.txt
    samtools flagstat ${bam.getSimpleName()}_filtered.bam > ${bam.getSimpleName()}_flagstat.after.txt
    """
}


process SAMTOOLS_COVERAGE {
    publishDir params.align_out_dir, mode: "copy", pattern: '*.cov', saveAs: {filename -> "coverage/${filename}"}

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
    params.input_csv = 'none' // csv file with paths to fastqs and references, "sample_id","fastq_1","fastq_2" and "reference" or "fastq" and "reference" needs to be in the header
    params.fq = '' // e.g. with quotes "*{1,2}.fq" for paired end
    params.fa = ''
    params.single_end = false

    params.align_out_dir = "Align_Reads_Output"

    if (params.input_csv == 'none') {
        if (params.fq == 'none' || params.fa == 'none') {
            println "Error: --fq and --fa parameters are required (or use --input_csv). Use --help for more information."
            exit 1
        }

        else if (params.single_end == true) {
            fastqs_channel = Channel.fromPath(params.fq)
            reference_genome_channel = Channel.fromPath(params.fa)
            combined_inputs = fastqs_channel.combine(reference_genome_channel)
        }
        
        else if (params.single_end == false) {
            fastqs_channel = Channel.fromFilePairs(params.fq)
            reference_genome_channel = Channel.fromPath(params.fa)
            combined_inputs = fastqs_channel.combine(reference_genome_channel).map { sample_id, fastqs, ref -> tuple(sample_id, fastqs[0], fastqs[1], ref) }
        }
    }

    else if (params.input_csv != 'none') {
        println "Using input csv file: ${params.input_csv}"

        if (params.single_end == true) {
            combined_inputs = Channel.fromPath(params.input_csv)
                .splitCsv(header:true)
                .map { row-> tuple(file(row.fastq), file(row.reference)) }
        }

        else if (params.single_end == false) {
            combined_inputs = Channel.fromPath(params.input_csv)
                .splitCsv(header:true)
                .map { row-> tuple(row.sample_id,file(row.fastq_1),file(row.fastq_2), file(row.reference)) }
        }
    }


    if (params.single_end == true) {
        // Process execution
        FILTER_FASTQ_SINGLE_END(combined_inputs)

        BWA_MEM_SINGLE_END(FILTER_FASTQ_SINGLE_END.out[0])

        FILTER_BAM(BWA_MEM_SINGLE_END.out)

        SAMTOOLS_COVERAGE(FILTER_BAM.out[0])
    }


    else if (params.single_end == false) {
        // Process execution
        FILTER_FASTQ_PAIRED_END(combined_inputs)

        BWA_MEM_PAIRED_END(FILTER_FASTQ_PAIRED_END.out[0])

        FILTER_BAM(BWA_MEM_PAIRED_END.out)

        SAMTOOLS_COVERAGE(FILTER_BAM.out[0])
    }

}