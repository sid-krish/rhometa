#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


def helpMessage() {

    log.info"""
    Description:

    Usage:
    nextflow run sim_gen.nf [options]

    Help:
    nextflow run sim_gen.nf --help

    Options:
    --recom_tract_len [int], default:[500], Recombination tract length to use
    --single_end, Used for single end read bams
    --read_len [int], default:[150], Read length of each individual read
    --paired_end_mean_frag_len [int], default:[300], The mean size of DNA fragments for paired-end simulations 
    --paired_end_std_dev [int], default:[50], The standard deviation of DNA fragment size for paired-end simulations 
    --seed [int], default:[123], Seed value to use for simulation
    --mutation_rate [int], default:[0.01], Population mutation rate, theta
    --rho_rates [int], default:[10], Population recombiation rate, rho 
    --sample_sizes [int], default:[10], Number of haplotypes to use for generating reads
    --genome_sizes [int], default:[10000], Genome size of haplotypes
    --fold_cov [int], default:[10], The fold of read coverage to be simulated or number of reads/read pairs generated for each haplotype genome
    --prepend_filename [str], Prepend string to output filenames to help distinguish runs

    """.stripIndent()

}


process RATE_SELECTOR {
    // This process prevents the need to use each in every process, which can be confusing
    // Perhaps this could be handled in a more elegant way using some DSL2 technique
    
    maxForks 1 // Run sequentially

    input:
        val rho_rate
        val sample_size
        val genome_size

    output:
        val "${rho_rate}", emit: p_val
        val "${sample_size}", emit: sample_size
        val "${genome_size}", emit: genome_size

    script:
    """
    """

}  


process MS {
    // publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

    maxForks 1

    input:
        val rho_rate
        val mutation_rate
        val seed
        val sample_size
        val genome_size
        val recom_tract_len
        val prepend_filename

    output:
        path "trees.txt", emit: trees_txt
  
    script:
    """
    ms ${sample_size} 1 -T -seeds ${seed} -t ${mutation_rate} -r ${rho_rate} ${genome_size} -c 10 ${recom_tract_len} > trees.txt
    """
}


process FAST_SIM_BAC {
    // publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

    maxForks 1
    
    input:
        val rho_rate
        val sample_size
        val seed
        val mutation_rate
        val recom_tract_len
        val genome_size
        val prepend_filename

    output:
        // path "rho_calc.txt", emit: rho_rho_calc_txt
        path "trees.txt", emit: trees_txt
             
    script:
    """
    fastSimBac ${sample_size} ${genome_size} -s ${seed} -T -t ${mutation_rate} -r ${rho_rate} ${recom_tract_len} > trees.txt
    """
}


process MS_PRIME {
    // publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

    maxForks 1
    
    input:
        val rho_rate
        val sample_size
        val genome_size
        val effective_pop_size
        val prepend_filename

    output:
        path "rho_calculation.txt", emit: rho_rho_calculation_txt
        path "trees.txt", emit: trees_txt
             
    script:
    """
    run_msprime.py
    """
}


process CLEAN_TREES {
    // publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

    maxForks 1

    input:
        path trees
        val prepend_filename

    output:
        path "cleanTrees.txt", emit: cleanTrees_txt

    script:
    """
    clean_trees.py trees.txt
    """
}


process SEQ_GEN {
    // publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

    maxForks 1

    input:
        path cleanTrees
        val genome_size
        val seed
        val prepend_filename

    output:
        path "seqgenOut.fa", emit: seqgenout_fa

    script:
    // 1 partition per tree
    // program crashes if seq length is not as the one set for fastsimbac
    """
    numTrees=\$(wc -l cleanTrees.txt | awk '{ print \$1 }')
    seq-gen -m HKY -t 4 -l ${genome_size} -z ${seed} -s 0.01 -p \$numTrees -of cleanTrees.txt > seqgenOut.fa
    """
}


process REFORMAT_FASTA {
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

    maxForks 1

    input:
        path seqgenOut
        val prepend_filename

    output:
        path "reformatted.fa", emit: reformatted_fa

    script:
    """
    samtools faidx seqgenOut.fa
    reformat_fasta_pysam.py seqgenOut.fa
    """
}


process ISOLATE_GENOME {
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

    maxForks 1

    input:
        path reformatted_fa
        val prepend_filename

    output:
        path "firstGenome.fa", emit: firstGenome_fa

    script:
    """
    #!/bin/bash
    head -2 reformatted.fa > firstGenome.fa
    """
}


process INSILICO_SEQ {
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

    maxForks 1

    input:
        path reformatted_fa
        val seed
        val prepend_filename
        
    output:
        path "*.fastq", emit: iss_out_fastq
        path "*abundance.txt"

    script:
    """
    iss generate --cpus $task.cpus --seed ${seed} --genomes reformatted.fa --model HiSeq --abundance lognormal --n_reads 20k  -o iss_out
    """
}


process BWA_MEM_PAIRED_END {
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

    maxForks 1

    input:
        path firstGenome_fa
        path art_out_fq
        val prepend_filename

    output:
        path "Aligned.bam", emit: aligned_bam

    script:
    // using first fa entry only (one genome)
    """
    bwa index firstGenome.fa

    #Paired end
    bwa mem -t $task.cpus firstGenome.fa iss_out_R1.fastq iss_out_R2.fastq > Aligned.sam

    samtools view -bS Aligned.sam > Aligned.bam
    """
}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // Params
    params.help = false
    params.prepend_filename = ""

    params.single_end = false
    params.read_len = 150
    params.paired_end_mean_frag_len = 300
    params.paired_end_std_dev = 25 // +- mean frag len

    params.seed = 123
    params.mutation_rate = 0.01
    params.recom_tract_len = 500
    params.effective_pop_size = 1 // only for msprime
    params.rho_rates = 0.05
    params.sample_sizes = 40
    params.genome_sizes = 25000
    params.fold_cov = 10

    // Input verification
    if (params.help) {
        // Show help message from helpMessage() function
        // params.help = false by default
        helpMessage()
        exit 0
    }
    
    // Process execution
    RATE_SELECTOR(params.rho_rates, params.sample_sizes, params.genome_sizes)

    // MS(RATE_SELECTOR.out.p_val, params.mutation_rate, params.seed, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, params.recom_tract_len, params.prepend_filename)

    FAST_SIM_BAC(RATE_SELECTOR.out.p_val, RATE_SELECTOR.out.sample_size, params.seed, params.mutation_rate, params.recom_tract_len, RATE_SELECTOR.out.genome_size, params.prepend_filename)

    // MS_PRIME(RATE_SELECTOR.out.p_val, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, params.effective_pop_size, params.prepend_filename)

    // CLEAN_TREES(MS.out.trees_txt, params.prepend_filename)

    CLEAN_TREES(FAST_SIM_BAC.out.trees_txt, params.prepend_filename)

    SEQ_GEN(CLEAN_TREES.out.cleanTrees_txt, RATE_SELECTOR.out.genome_size, params.seed, params.prepend_filename)

    REFORMAT_FASTA(SEQ_GEN.out.seqgenout_fa, params.prepend_filename)

    ISOLATE_GENOME(REFORMAT_FASTA.out.reformatted_fa, params.prepend_filename)

    INSILICO_SEQ(REFORMAT_FASTA.out.reformatted_fa, params.seed, params.prepend_filename)
    
    // Query name sorted bam
    BWA_MEM_PAIRED_END(ISOLATE_GENOME.out.firstGenome_fa, INSILICO_SEQ.out.iss_out_fastq, params.prepend_filename)

}