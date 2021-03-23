#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process RATE_SELECTOR {
    // This process prevents the need to use each in every process, which can be confusing
    // Perhaps this could be handled in a more elegant way using some DSL2 technique
    
    maxForks 1 // Run sequentially

    input:
        each rho_rate
        each sample_size
        each genome_size

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
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1 

    input:
        val rho_rate
        val sample_size
        val genome_size
        val path_fn_modifier

    output:
        path "trees.txt", emit: trees_txt
  
    script:
    """
    echo 123 456 789 > seedms
    ms ${sample_size} 1 -T -t ${params.mutation_rate} -r ${rho_rate} ${genome_size} -c 10 ${params.recom_tract_len} > trees.txt
    """
}


process FAST_SIM_BAC {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1 
    
    input:
        val rho_rate
        val sample_size
        val genome_size
        val path_fn_modifier

    output:
        path "rho_calc.txt", emit: rho_rho_calc_txt
        path "trees.txt", emit: trees_txt
             
    script:
    """
    fastSimBac ${sample_size} ${genome_size} -s ${params.seed} -T -t ${params.mutation_rate} -r ${rho_rate} ${params.recom_tract_len} > trees.txt
    calc_rho.py ${params.effective_pop_size} ${rho_rate} ${genome_size}
    """
}


process MS_PRIME {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1 
    
    input:
        val rho_rate
        val sample_size
        val genome_size
        val path_fn_modifier

    output:
        path "rho_calculation.txt", emit: rho_rho_calculation_txt
        path "trees.txt", emit: trees_txt
             
    script:
    """
    run_msprime.py
    calc_rho.py ${params.effective_pop_size} ${rho_rate} ${genome_size}
    """
}


process CLEAN_TREES {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}
    maxForks 1

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
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path cleanTrees
        val genome_size
        val path_fn_modifier

    output:
        path "seqgenOut.fa", emit: seqgenout_fa

    script:
    // 1 partiion per tree
    // program crashes if seq length is not as the one set for fastsimbac
    """
    numTrees=\$(wc -l cleanTrees.txt | awk '{ print \$1 }')
    seq-gen -m HKY -t 4 -l ${genome_size} -z ${params.seed} -s 0.01 -p \$numTrees -of cleanTrees.txt > seqgenOut.fa
    """
}


process REFORMAT_FASTA {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path seqgenOut
        val path_fn_modifier

    output:
        path "reformatted.fa", emit: reformatted_fa

    script:
    """
    reformat_fasta.py seqgenOut.fa
    """
}


process ISOLATE_GENOME{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path reformatted_fa
        val path_fn_modifier

    output:
        path "firstGenome.fa", emit: firstGenome_fa

    script:
    """
    isolateGenome.py reformatted.fa
    """
}


process ART_ILLUMINA{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path reformatted_fa
        val path_fn_modifier

    output:
        path "art_fastSimBac.fq", emit: art_fastSimBac_fq // only file we are interested in

    script:
    """
    art_illumina --seqSys HSXt --rndSeed ${params.seed} --noALN \
    --in reformatted.fa --len ${params.meanFragmentLen} --fcov 2 --out art_fastSimBac
    """
    // go with single end reads initially to make things easier
    // mflen should be around 500, sdev around 50-60
}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

   
    // params.meanFragmentLen = 150
    params.seed = 123
    params.mutation_rate = 0.01
    params.recom_tract_len = 500
    params.ldpop_rho_range = "101,100"
    params.effective_pop_size = 1
    params.meanFragmentLen = 150
    
    // precomputed likelihood table
    // lookup_Table = Channel.fromPath("$baseDir/lookupTable.txt")
    
    trees = Channel.fromPath("$baseDir/trees.txt")
    
    rho_rates = Channel.from(15) // For fastsimbac use this for recom rate (it doesn't accept rho)
    sample_sizes = Channel.from(10)
    genome_sizes = Channel.from(20000)
    
    RATE_SELECTOR(rho_rates, sample_sizes, genome_sizes)

    MS(RATE_SELECTOR.out.p_val, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    // FAST_SIM_BAC(RATE_SELECTOR.out.p_val, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    // MS_PRIME(RATE_SELECTOR.out.p_val, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    CLEAN_TREES(MS.out.trees_txt, RATE_SELECTOR.out.path_fn_modifier)

    // CLEAN_TREES(FAST_SIM_BAC.out.trees_txt, RATE_SELECTOR.out.path_fn_modifier)

    // CLEAN_TREES(trees, RATE_SELECTOR.out.path_fn_modifier)

    SEQ_GEN(CLEAN_TREES.out.cleanTrees_txt, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    REFORMAT_FASTA(SEQ_GEN.out.seqgenout_fa, RATE_SELECTOR.out.path_fn_modifier)

    ISOLATE_GENOME(REFORMAT_FASTA.out.reformatted_fa, RATE_SELECTOR.out.path_fn_modifier)

    ART_ILLUMINA(REFORMAT_FASTA.out.reformatted_fa, RATE_SELECTOR.out.path_fn_modifier)

}