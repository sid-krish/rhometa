#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


def helpMessage() {

    log.info"""
    Usage:
    nextflow run theta_est.nf --bam in.bam --fa ref.fa [options]

    Help:
    nextflow run theta_est.nf --help

    Required:
    --bam [*.bam], Query name sorted bam file. Multi bam support via glob input e.g. "*.bam", quotes must be included for glob. Use with one fasta file only
    --fa [*.fa],  Single/Multi genome fasta file

    Options:
    --filename_prefix [str], prefix string to output filenames to help distinguish runs
    --output_dir [str], default:[Theta_Est_Output], Directory to save results in
    """.stripIndent()
    // Testing
    // --snp_qual [int], default:[20], Minimum phred-scaled quality score to filter vcf by
    // --min_snp_depth [int], default:[10], Minimum read depth to filter vcf by

}


process FILENAME_PREFIX {
    
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
    filename_prefix.py ${bam} ${prefix_fn} 
    """
}


process FILTER_BAM {
    /**
      * Filter input BAM files to assure that reported alignments will be of acceptable quality.
      * - filters on mapping quality and alignment-score relative to read-length.
      * - publishes a flagstat report of before and after.
      **/

    publishDir params.output_dir, mode: 'copy', pattern: 'flagstat.*.txt', saveAs: {filename -> "filter_bam/${filename_prefix}${filename}"}

    input:
        tuple val(filename_prefix),
            path(bam),
            path(fasta)
    
    output:
        tuple val(filename_prefix),
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
    tuple val(filename_prefix), 
        path(bam), 
        path(fasta)

    output:
    tuple val(filename_prefix), 
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
            path("freebayes_filt.vcf")
        
        path 'freebayes_raw.vcf'

    script:
    """
    samtools faidx ${fasta}
    samtools index -@ $task.cpus ${bam}

    # call variants with freebayes
    freebayes -f ${fasta} -p 1 ${bam} > freebayes_raw.vcf

    bcftools filter --threads ${task.cpus} -i 'TYPE="snp"' freebayes_raw.vcf > freebayes_filt.vcf
    """
    // Testing
    // keep SNPs and remove low quality and low depth calls
    // bcftools filter --threads ${task.cpus} \
    //     -i 'TYPE="snp" && QUAL>=${params.snp_qual} && FORMAT/DP>=${params.min_snp_depth} && FORMAT/RO>=2 && FORMAT/AO>=2' freebayes_raw.vcf > freebayes_filt.vcf
}


process THETA_ESTIMATE {
    publishDir params.output_dir, mode: "copy", saveAs: {filename -> "${filename_prefix}${filename}"}

    // maxForks 1

    input:
        tuple val(filename_prefix),
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
    params.filename_prefix = "none"
    // VCF filter settings
    // params.snp_qual = 20 // Minimum phred-scaled quality score to filter vcf by
    // params.min_snp_depth = 10 // Minimum read depth to filter vcf by

    params.output_dir = 'Theta_Est_Output'
    params.bam = 'none'
    params.fa = 'none'

    // Channels
    bam_channel = Channel.fromPath( params.bam, checkIfExists: true )
    fa_channel = Channel.fromPath( params.fa, checkIfExists: true )

    bam_and_fa = bam_channel.combine(fa_channel)

    // Input verification
    if (params.help) {
        // Show help message from helpMessage() function
        // params.help = false by default
        helpMessage()
        exit 0
    }

    if (params.fa == 'none') {
        println "No input .fa specified. Use --fa [.fa]"
        exit 1
    }

    if (params.bam == 'none') {
        println "No input .bam specified. Use --bam [.bam]"
        exit 1
    }

    // Process execution
    FILENAME_PREFIX(bam_and_fa, params.filename_prefix)

    FILTER_BAM(FILENAME_PREFIX.out)

    SORT_BAM(FILTER_BAM.out[0])

    FREEBAYES(SORT_BAM.out)

    // freebayes returns two channels, we just need the first
    THETA_ESTIMATE(FREEBAYES.out[0])

}