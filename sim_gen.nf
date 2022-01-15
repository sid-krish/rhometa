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

    """.stripIndent()

}


process RATE_SELECTOR {
    // This process prevents the need to use each in every process, which can be confusing
    // Perhaps this could be handled in a more elegant way using some DSL2 technique
    
    // maxForks 1 // Run sequentially
    
    // echo true

    input:
        each rho
        each theta
        each sample_size
        each depth
        each genome_size
        each seed

    output:
        tuple stdout,
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed)

    script:
    """
    printf rho_${rho}_theta_${theta}_sample_size_${sample_size}_depth_${depth}_genome_size_${genome_size}_seed_${seed}_
    """

}  


process FAST_SIM_BAC {
    // publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1
    
    input:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed)

        val(recom_tract_len)

    output:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("trees.txt")
             
    script:
    """
    fastSimBac ${sample_size} ${genome_size} -s ${seed} -T -t ${theta} -r ${rho} ${recom_tract_len} > trees.txt
    """
}


process CLEAN_TREES {
    // publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    input:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("trees.txt")

    output:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("cleanTrees.txt")

    script:
    """
    clean_trees.py trees.txt
    """
}


process SEQ_GEN {
    // publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    input:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("cleanTrees.txt")

    output:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("seqgenOut.fa")

    script:
    // 1 partition per tree
    // program crashes if seq length is not as the one set for fastsimbac
    """
    numTrees=\$(wc -l cleanTrees.txt | awk '{ print \$1 }')
    seq-gen -m HKY -t 4 -l ${genome_size} -z ${seed} -s 0.01 -p \$numTrees -of cleanTrees.txt > seqgenOut.fa
    """
}


process REFORMAT_FASTA {
    // publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    input:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("seqgenOut.fa")

    output:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("reformatted.fa")

    script:
    """
    samtools faidx seqgenOut.fa
    reformat_fasta_pysam.py seqgenOut.fa
    """
}


process ISOLATE_GENOME {
    // publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    input:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("reformatted.fa")

    output:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("reformatted.fa"),
            path("firstGenome.fa")

    script:
    """
    #!/bin/bash
    head -2 reformatted.fa > firstGenome.fa
    """
}


process ART_ILLUMINA_SINGLE_END {
    // publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    input:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("reformatted.fa"),
            path("firstGenome.fa")

        val(read_len)
        val(paired_end_std_dev)
        val(paired_end_mean_frag_len)

    output:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("firstGenome.fa"),
            path("art_out.fq")

    script:
    """
    #Single end
    art_illumina --seqSys HSXt --rndSeed ${seed} --noALN --quiet \
    --in reformatted.fa --len ${read_len} --fcov ${depth} --out art_out
    """
}


process ART_ILLUMINA_PAIRED_END {
    // publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    input:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("reformatted.fa"),
            path("firstGenome.fa")

        val(read_len)
        val(paired_end_std_dev)
        val(paired_end_mean_frag_len)

    output:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("firstGenome.fa"),
            path("art_out_1.fq"),
            path("art_out_2.fq")
             

    script:
    """
    #Paired end
    art_illumina --seqSys HSXt --rndSeed ${seed} --noALN --quiet \
    --in reformatted.fa -p --len ${read_len} --sdev ${paired_end_std_dev} \
    -m ${paired_end_mean_frag_len} --fcov ${depth} --out art_out_
    """
}


process BWA_MEM_SINGLE_END {
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    input:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("firstGenome.fa"),
            path("art_out.fq")

    output:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("final.fa"),
            path("final.bam")

    script:
    // using first fa entry only (one genome)
    """
    bwa index firstGenome.fa

    #Single end
    bwa mem -t $task.cpus firstGenome.fa art_out.fq > Aligned.sam

    samtools view -bS Aligned.sam > Aligned.bam

    mv Aligned.bam final.bam
    mv firstGenome.fa final.fa
    """
}


process BWA_MEM_PAIRED_END {
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    // maxForks 1

    input:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("firstGenome.fa"),
            path("art_out_1.fq"),
            path("art_out_2.fq")

    output:
        tuple val(prefix_filename),
            val(rho),
            val(theta),
            val(sample_size),
            val(depth),
            val(genome_size),
            val(seed),
            path("final.fa"),
            path("final.bam")

    script:
    // using first fa entry only (one genome)
    """
    bwa index firstGenome.fa

    #Paired end
    bwa mem -t $task.cpus firstGenome.fa art_out_1.fq art_out_2.fq > Aligned.sam

    samtools view -bS Aligned.sam > Aligned.bam

    mv Aligned.bam final.bam
    mv firstGenome.fa final.fa
    """
}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // For each process there is a output of tuple with the necessary files/values to move forward until they are no longer need

    // Params
    params.help = false

    params.recom_tract_len = 1000
    params.single_end = false
    params.read_len = 150
    params.paired_end_mean_frag_len = 300
    params.paired_end_std_dev = 25 // +- mean frag len

    params.rho_rates = [0.0001, 0.00025, 0.0005, 0.00075]
    params.theta_rates = [0.01]
    params.sample_sizes = [10, 20, 40, 80, 160, 200]
    params.fold_cov_rates = [1, 2, 4, 8]
    params.genome_sizes = [100000]
    params.seed_vals = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    // Input verification
    if (params.help) {
        // Show help message from helpMessage() function
        // params.help = false by default
        helpMessage()
        exit 0
    }
    
    // Process execution
    RATE_SELECTOR(params.rho_rates, params.theta_rates, params.sample_sizes, params.fold_cov_rates, params.genome_sizes, params.seed_vals)

    FAST_SIM_BAC(RATE_SELECTOR.out, params.recom_tract_len)

    CLEAN_TREES(FAST_SIM_BAC.out)

    SEQ_GEN(CLEAN_TREES.out)

    REFORMAT_FASTA(SEQ_GEN.out)

    ISOLATE_GENOME(REFORMAT_FASTA.out)

    if (params.single_end == true) {
        ART_ILLUMINA_SINGLE_END(ISOLATE_GENOME.out, params.read_len, params.paired_end_std_dev, params.paired_end_mean_frag_len)
        // Query name sorted bam
        BWA_MEM_SINGLE_END(ART_ILLUMINA_SINGLE_END.out)
    }
    
    else if (params.single_end == false) {
        ART_ILLUMINA_PAIRED_END(ISOLATE_GENOME.out, params.read_len, params.paired_end_std_dev, params.paired_end_mean_frag_len)
        // Query name sorted bam
        BWA_MEM_PAIRED_END(ART_ILLUMINA_PAIRED_END.out)
    }

}