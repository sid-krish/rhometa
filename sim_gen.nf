#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process RATE_SELECTOR {
    // This process prevents the need to use each in every process, which can be confusing
    // Perhaps this could be handled in a more elegant way using some DSL2 technique
    
    maxForks 1 // Run sequentially

    cpus 1
    // memory '100 MB'

    // executor 'local'
    // time '5m'
    // scratch true
    // queue 'i3q'

    input:
        val rho_rate
        val sample_size
        val genome_size

    output:
        val "${rho_rate}", emit: p_val
        val "${sample_size}", emit: sample_size
        val "${genome_size}", emit: genome_size

        val "rho_${rho_rate}_sam_${sample_size}_gen_${genome_size}/rho_${rho_rate}_sam_${sample_size}_gen_${genome_size}", emit: path_fn_modifier

    script:
    """
    """

}  


process MS {
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    // memory '128 MB'

    // executor 'local'
    // time '15m'
    // scratch true
    // queue 'i3q'

    input:
        val rho_rate
        val mutation_rate
        val sample_size
        val genome_size
        val recom_tract_len
        val path_fn_modifier

    output:
        path "trees.txt", emit: trees_txt
  
    script:
    """
    echo 123 456 789 > seedms
    ms ${sample_size} 1 -T -t ${mutation_rate} -r ${rho_rate} ${genome_size} -c 10 ${recom_tract_len} > trees.txt
    """
}


process FAST_SIM_BAC {
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    // memory '128 MB'

    // executor 'local'
    // time '1h'
    // scratch true
    // queue 'i3q'
    
    input:
        val rho_rate
        val sample_size
        val seed
        val mutation_rate
        val recom_tract_len
        val genome_size
        val path_fn_modifier

    output:
        // path "rho_calc.txt", emit: rho_rho_calc_txt
        path "trees.txt", emit: trees_txt
             
    script:
    """
    fastSimBac ${sample_size} ${genome_size} -s ${seed} -T -t ${mutation_rate} -r ${rho_rate} ${recom_tract_len} > trees.txt
    """
}


process MS_PRIME {
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    // memory '128 MB'

    // executor 'local'
    // time '1h'
    // scratch true
    // queue 'i3q'
    
    input:
        val rho_rate
        val sample_size
        val genome_size
        val effective_pop_size
        val path_fn_modifier

    output:
        path "rho_calculation.txt", emit: rho_rho_calculation_txt
        path "trees.txt", emit: trees_txt
             
    script:
    """
    run_msprime.py
    """
}


process CLEAN_TREES {
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    // memory '128 MB'

    // executor 'local'
    // time '10m'
    // scratch true
    // queue 'i3q'

    input:
        path trees
        val path_fn_modifier


    output:
        path "cleanTrees.txt", emit: cleanTrees_txt

    script:
    """
    clean_trees.py trees.txt
    """
}


process SEQ_GEN {
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    // memory '128 MB'

    // executor 'local'
    // time '10m'
    // scratch true
    // queue 'i3q'

    input:
        path cleanTrees
        val genome_size
        val seed
        val path_fn_modifier

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
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    // memory '128 MB'

    // executor 'local'
    // time '10m'
    // scratch true
    // queue 'i3q'

    input:
        path seqgenOut
        val path_fn_modifier

    output:
        path "reformatted.fa", emit: reformatted_fa

    script:
    """
    samtools faidx seqgenOut.fa
    reformat_fasta_pysam.py seqgenOut.fa
    """
}


process ISOLATE_GENOME {
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    // memory '128 MB'

    // executor 'local'
    // time '5m'
    // scratch true
    // queue 'i3q'

    input:
        path reformatted_fa
        val path_fn_modifier

    output:
        path "firstGenome.fa", emit: firstGenome_fa

    script:
    """
    #!/bin/bash
    head -2 reformatted.fa > firstGenome.fa
    """
}


process ART_ILLUMINA_SINGLE_END {
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    // memory '128 MB'

    // executor 'local'
    // time '1h'
    // scratch true
    // queue 'i3q'

    input:
        path reformatted_fa
        val seed
        val read_len
        val path_fn_modifier

    output:
        path "*.fq", emit: art_out_fq // only files we are interested in
        // path "*.aln" // for testing

    script:
    """
    #Single end
    art_illumina --seqSys HSXt --rndSeed ${seed} --noALN --quiet \
    --in reformatted.fa --len ${read_len} --fcov 10 --out art_out
    """
}


process ART_ILLUMINA_PAIRED_END {
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    // memory '128 MB'

    // executor 'local'
    // time '1h'
    // scratch true
    // queue 'i3q'

    input:
        path reformatted_fa
        val seed
        val read_len
        val paired_end_std_dev
        val paired_end_mean_frag_len
        val path_fn_modifier

    output:
        path "*.fq", emit: art_out_fq // only files we are interested in
        // path "*.aln" // for testing

    script:
    """
    #Paired end
    art_illumina --seqSys HSXt --rndSeed ${seed} --noALN --quiet \
    --in reformatted.fa -p --len ${read_len} --sdev ${paired_end_std_dev} \
    -m ${paired_end_mean_frag_len} --fcov 10 --out art_out
    """
}


process BWA_MEM_SINGLE_END {
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus {num_cores}
    // memory '1 GB'

    // executor 'local'
    // time '30m'
    // scratch true
    // queue 'i3q'

    input:
        path firstGenome_fa
        path art_out_fq
        val num_cores
        val path_fn_modifier

    output:
        path "Aligned.bam", emit: aligned_bam

    script:
    // using first fa entry only (one genome)
    """
    bwa index firstGenome.fa

    #Single end
    bwa mem -t ${num_cores} firstGenome.fa art_out.fq > Aligned.sam

    samtools view -bS Aligned.sam > Aligned.bam
    """
}


process BWA_MEM_PAIRED_END {
    publishDir "Sim_Gen_Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus {num_cores}
    // memory '1 GB'

    // executor 'local'
    // time '30m'
    // scratch true
    // queue 'i3q'

    input:
        path firstGenome_fa
        path art_out_fq
        val num_cores
        val path_fn_modifier

    output:
        path "Aligned.bam", emit: aligned_bam

    script:
    // using first fa entry only (one genome)
    """
    bwa index firstGenome.fa

    #Paired end
    bwa mem -t ${num_cores} firstGenome.fa art_out1.fq art_out2.fq > Aligned.sam

    samtools view -bS Aligned.sam > Aligned.bam
    """
}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    params.num_cores = 4
   
    params.seed = 123
    params.mutation_rate = 0.01
    params.recom_tract_len = 500
    params.effective_pop_size = 1

    params.seq_type = 1  // 0 - single end, 1 - paired end
    params.read_len = 150
    params.paired_end_mean_frag_len = 300
    params.paired_end_std_dev = 50 // +- mean frag len

    params.rho_rates = 10  // Using each multiple times, causes nextflow to work improperly so removed it here for use with depth
    params.sample_sizes = 10
    params.genome_sizes = 10000
    
    RATE_SELECTOR(params.rho_rates, params.sample_sizes, params.genome_sizes)

    MS(RATE_SELECTOR.out.p_val, params.mutation_rate, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, params.recom_tract_len, RATE_SELECTOR.out.path_fn_modifier)

    // FAST_SIM_BAC(RATE_SELECTOR.out.p_val, RATE_SELECTOR.out.sample_size, params.seed, params.mutation_rate, params.recom_tract_len, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    // MS_PRIME(RATE_SELECTOR.out.p_val, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, params.effective_pop_size, RATE_SELECTOR.out.path_fn_modifier)

    CLEAN_TREES(MS.out.trees_txt, RATE_SELECTOR.out.path_fn_modifier)

    // CLEAN_TREES(FAST_SIM_BAC.out.trees_txt, RATE_SELECTOR.out.path_fn_modifier)

    SEQ_GEN(CLEAN_TREES.out.cleanTrees_txt, RATE_SELECTOR.out.genome_size, params.seed, RATE_SELECTOR.out.path_fn_modifier)

    REFORMAT_FASTA(SEQ_GEN.out.seqgenout_fa, RATE_SELECTOR.out.path_fn_modifier)

    ISOLATE_GENOME(REFORMAT_FASTA.out.reformatted_fa, RATE_SELECTOR.out.path_fn_modifier)

    if (params.seq_type == 0) {
        ART_ILLUMINA_SINGLE_END(REFORMAT_FASTA.out.reformatted_fa, params.seed, params.read_len, RATE_SELECTOR.out.path_fn_modifier)

        BWA_MEM_SINGLE_END(ISOLATE_GENOME.out.firstGenome_fa, ART_ILLUMINA_SINGLE_END.out.art_out_fq, params.num_cores, RATE_SELECTOR.out.path_fn_modifier)
    }
    
    else if (params.seq_type == 1) {
        ART_ILLUMINA_PAIRED_END(REFORMAT_FASTA.out.reformatted_fa, params.seed, params.read_len, params.paired_end_std_dev, params.paired_end_mean_frag_len, RATE_SELECTOR.out.path_fn_modifier)

        BWA_MEM_PAIRED_END(ISOLATE_GENOME.out.firstGenome_fa, ART_ILLUMINA_PAIRED_END.out.art_out_fq, params.num_cores, RATE_SELECTOR.out.path_fn_modifier)
    }

}