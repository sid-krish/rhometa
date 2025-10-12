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
    --theta_est_out_dir [str], default:[Theta_Est_Output], Directory to save results in
    --snp_qual [int], default:[20], Minimum phred-scaled quality score to filter vcf by
    --min_snp_depth [int], default:[10], Minimum read depth to filter vcf by
    """.stripIndent()
}


process FILENAME_PREFIX {
    
    // maxForks 1

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

    script:
    """
    samtools sort -@ $task.cpus -o Aligned_sorted.bam ${bam}
    """
}


process ANGSD_THETA_ESTIMATE {
    publishDir params.theta_est_out_dir, mode: "copy", pattern: "*angsd_theta_*", saveAs: {filename -> "theta_estimate/${filename}"}

    // maxForks 1

    input:
        tuple val(filename_prefix),
            path(bam),
            path(fasta)

    output:
        tuple val(filename_prefix),
            path(bam),
            path(fasta),
            path("${filename_prefix}angsd_theta_raw.tsv"),
            path("${filename_prefix}angsd_theta_final.csv")

    script:
    """
    samtools faidx ${fasta}
    samtools index -@ $task.cpus ${bam}

    # Commands follow steps from https://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests

    #First estimate the site allele frequency likelihood
    angsd -i ${bam} -doSaf 1 -anc ${fasta} -GL 1 -P $task.cpus -out out #output is log scale

    #For folded " If you do not have the ancestral state you can simply use the assembly you have mapped agains, but remember to add -fold 1 in the 'realSFS' and 'realSFS sf2theta' step."

    # Obtain the maximum likelihood estimate of the SFS using the realSFS program
    realSFS out.saf.idx -P $task.cpus -fold 1 > out.sfs

    # Calculate the thetas for each site
    realSFS saf2theta out.saf.idx -outname out -sfs out.sfs -fold 1 

    # Estimate Tajimas D and other statistics
    thetaStat do_stat out.thetas.idx # no longer log scale. Final vals. pestPG file are the sum of the per site estimates for a region -> chromosome wide vals

    # pestPG file has the required output which is converted to TSV format for nextsteps
    # cat out.thetas.idx.pestPG | tr '\t' ',' > angsd_theta_raw.csv
    mv out.thetas.idx.pestPG ${filename_prefix}angsd_theta_raw.tsv

    #New solution works for bams aligned to single/multiple sequences by summing the lengths of all @SQ records.
    genome_size=\$(samtools view -H ${bam} |  awk '/^@SQ/ {l+=substr(\$3,4)}END{print l}')

    # Python script for genome-wide per-site Watterson theta from ANGSD output.
    compute_final_angsd_theta.py ${filename_prefix}angsd_theta_raw.tsv \$genome_size > ${filename_prefix}angsd_theta_final.csv
    
    # Sliding window analysis. Not needed. Testing
    # thetaStat do_stat out.thetas.idx -win 50000 -step 10000 -outnames theta.thetasWindow.gz
    """
}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // For each process there is a output of tuple with the necessary files/values to move forward until they are no longer need

    // Params
    params.help = false
    params.filename_prefix = "none"

    params.theta_est_out_dir = 'Theta_Est_Output'
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
    FILENAME_PREFIX(bam_and_fa,
                    params.filename_prefix)

    SORT_BAM(FILENAME_PREFIX.out)

    ANGSD_THETA_ESTIMATE(SORT_BAM.out)

}