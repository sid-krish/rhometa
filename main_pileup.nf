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
        // path "rho_calc.txt", emit: rho_rho_calc_txt
        path "trees.txt", emit: trees_txt
             
    script:
    """
    fastSimBac ${sample_size} ${genome_size} -s ${params.seed} -T -t ${params.mutation_rate} -r ${rho_rate} ${params.recom_tract_len} > trees.txt
    # calc_rho.py ${params.effective_pop_size} ${rho_rate} ${genome_size}
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


process ISOLATE_GENOME {
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


process ART_ILLUMINA {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path reformatted_fa
        val path_fn_modifier

    output:
        path "art_fastSimBac.fq", emit: art_fastSimBac_fq // only file we are interested in
        // path "*.aln" // for testing

    script:
    // can try --rcount as an alternative to --fcov
    // number of reads/read pairs to be generated per sequence/amplicon (not be used together with -f/--fcov)
    """
    art_illumina --seqSys HSXt --rndSeed ${params.seed} --noALN \
    --in reformatted.fa --len ${params.meanFragmentLen} --fcov 10 --maxIndel 0 --out art_fastSimBac
    #art_illumina --rndSeed ${params.seed} --noALN \
    #--in reformatted.fa --len ${params.meanFragmentLen} --fcov 20 --out art_fastSimBac
    """
    // go with single end reads initially to make things easier
    // mflen should be around 500, sdev around 50-60
}


process BWA_MEM {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    // cpus 8

    maxForks 1

    input:
        path firstGenome_fa
        path art_fastSimBac_fq
        val path_fn_modifier

    output:
        path "Aligned.bam", emit: aligned_bam

    script:
    // using first fa entry only (one genome)
    """
    bwa index firstGenome.fa
    bwa mem -t 4 firstGenome.fa art_fastSimBac.fq > Aligned.sam
    samtools view -bS Aligned.sam > Aligned.bam
    """
    // HPC will consider this to be using 8 threads
}


process PROCESS_SORT_INDEX{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1
    
    input:
        path aligned_bam
        val path_fn_modifier

    output:
        path "Aligned.csorted_fm_md.bam", emit: processed_bam
        path "Aligned.csorted_fm_md.bam.bai", emit: processed_index
        // path "Aligned.csorted.bam", emit: processed_bam
        // path "Aligned.csorted.bam.bai", emit: processed_index
        path "bam_stats.txt", emit: bam_stats_txt

    script:
    // samtools fixmate requires queryname sorting (-n option)
    // samtools fixmate -r removes secondary and unmapped reads
    // samtools fixmate -m adds ms (mate score) tags. These are used by markdup to select the best reads to keep.
    
    // samtools markdup requires chromosomal coordinate sorting
    // samtools markdup -r remove duplicates
    // samtools markdup -l expected maximum read length of INT bases.
    """
    #!/bin/bash
    bam_file_name=\$(echo ${aligned_bam} | cut -d. -f1)

    samtools stats --threads 4 "\$bam_file_name".bam  > bam_stats.txt
    max_read_len=\$(grep "maximum length" bam_stats.txt | cut -f 3)

    samtools fixmate --threads 4 -r -m  "\$bam_file_name".bam "\$bam_file_name".qsorted_fm.bam

    samtools sort --threads 4 "\$bam_file_name".qsorted_fm.bam -o "\$bam_file_name".csorted_fm.bam
    samtools markdup --threads 4 -r -l \$max_read_len "\$bam_file_name".csorted_fm.bam "\$bam_file_name".csorted_fm_md.bam

    samtools index -@ 4 "\$bam_file_name".csorted_fm_md.bam

    #samtools sort --threads 4 "\$bam_file_name".bam -o "\$bam_file_name".csorted.bam
    #samtools index -@ 4 "\$bam_file_name".csorted.bam
    """
}


process LOFREQ{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path firstGenome_fa
        path processed_bam
        path processed_index
        val path_fn_modifier

    output:
        path "lofreqOut.vcf", emit: lofreqOut_vcf

    script:
    """
    samtools faidx firstGenome.fa
    #lofreq call -f firstGenome.fa -o lofreqOut.vcf Aligned.csorted_fm_md.bam
    lofreq call-parallel --pp-threads 4 --no-default-filter -f firstGenome.fa -o lofreqOut.vcf Aligned.csorted_fm_md.bam
    """
}


process PILEUP_PAIRWISE_TABLE{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path lofreqOut_vcf 
        path bam_file
        path bam_index
        val path_fn_modifier

    output:
        path "pairwise_table.csv", emit: pairwise_table_csv

    script:
    """
    samtools view -h Aligned.csorted_fm_md.bam > Aligned.csorted_fm_md.sam
    pileup_table.py Aligned.csorted_fm_md.sam
    pileup_table_var_sites.py pileup_table.csv
    pileup_pairwise_table.py ${params.meanFragmentLen} variants_from_pileup.csv ${lofreqOut_vcf}
    """
}


process PAIRWISE_RESAMPLE{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    echo true

    input:
        path pairwise_table_csv
        val sample_size
        val path_fn_modifier

    output:
        path "pairwise_resampled.csv", emit: pairwise_resampled_csv

    script:
    """
    pairwise_resample_v2.py pairwise_table.csv ${params.seed} ${sample_size}
    """
}


process PAIRWISE_BIALLELIC_TABLE{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path pairwise_resampled_csv
        val path_fn_modifier

    output:
        path "pairwise_biallelic_table.csv", emit: pairwise_biallelic_table_csv

    script:
    // NOTE
    // There was serious error which prevented it from working properly, but is now fixed.
    // Need to fix in other projects that use this also
    """
    biallelic_filter.py pairwise_resampled.csv
    """
}


process LOOKUP_TABLE_LDPOP {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1 
    
    input:
        val sample_size
        val path_fn_modifier

    output:
        path "lookupTable.txt", emit: lookupTable_txt

    script:
    // There are other parameters that can be adjusted, I've left them out for the time being
    // also they mention twice muation and recom rate, for the mutation and recom parameters which I am unsure how to interpret
    """
    ldtable.py --cores 4 -n ${sample_size} -th ${params.mutation_rate} -rh ${params.ldpop_rho_range} --approx > lookupTable.txt
    """
}


process PAIRWISE_LOOKUP_FORMAT {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path pairwise_biallelic_table_csv
        val path_fn_modifier

    output:
        path "lookup_format.csv", emit: lookup_format_csv

    script:
    """
    pairwise_lookup_format_pyrho.py pairwise_biallelic_table.csv
    """
}


process CUSTOM_HAP_SETS_AND_MERGE {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path lookup_table_txt
        path pairwise_biallelic_table_csv
        path lookup_format_csv
        val path_fn_modifier

    output:
        path "table_ids_for_eq3.csv", emit: table_ids_for_eq3_csv
        path "eq3.csv", emit: eq3_csv

    script:
    """
    custom_hap_sets_and_merge.py lookupTable.txt pairwise_biallelic_table.csv lookup_format.csv ${params.ldpop_rho_range} > table_ids_for_eq3.csv
    """
}


process WATTERSON_ESTIMATE {
    
    maxForks 1

    input:
        path lofreqOut_vcf
        val sample_size
        val genome_size

    output:
        stdout emit: theta

    script:
    """
    watterson_estimate.py lofreqOut.vcf ${genome_size} ${sample_size}
    """
}


process P_IJ_GRID {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path eq3_csv
        val genome_size
        val path_fn_modifier

    output:
        path "p_ij_grid.csv", emit: p_ij_grid_csv

    script:
    """
    pij_grid_vectorised.py ${genome_size} ${params.recom_tract_len} ${params.ldpop_rho_range} eq3.csv
    """
}


process PAIRWISE_ESTIMATOR {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    // errorStrategy 'ignore'

    // echo true
    
    maxForks 1
    
    input:
        path eq3_csv
        path p_ij_grid_csv
        path table_ids_for_eq3_csv
        path lookup_table
        val path_fn_modifier

    output:
        path "collected_likelihoods.csv", emit: collected_likelihoods_csv
    
    script:
    """
    pairwise_rho_estimator_intp_rect_biv.py eq3.csv table_ids_for_eq3.csv p_ij_grid.csv lookupTable.txt ${params.ldpop_rho_range}
    """

}


process FINAL_RESULTS {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    // echo true

    maxForks 1

    input:
        path collectedFile
        val theta
        val path_fn_modifier

    output:
        path "final_results.txt", emit: final_results_txt

    script:
    """
    final_results.py collected_likelihoods.csv ${theta}
    """
}


process PROCESS_OUTPUT{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    input:
        path final_results_txt
        val rho_rate
        val sample_size
        val genome_size
        val path_fn_modifier

    output:
        path "processed_results.csv", emit: processed_results_csv

    script:
        """
        custom_est_process_output.py final_results.txt ${rho_rate} ${sample_size} ${genome_size}
        """

}


process PLOT_RESULTS{
    publishDir "Output/Results", mode: "copy"

    maxForks 1

    input:
        path collectedFile


    output:
        path "rho_comparision.png", emit: rho_comparision_png
        path "max_lk_comparision.png", emit: max_lk_comparision_png

    script:
        """
        plot_results.py collected_results.csv
        """

}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    params.seed = 123
    params.mutation_rate = 0.01
    params.recom_tract_len = 500
    params.ldpop_rho_range = "101,100"
    params.effective_pop_size = 1
    params.meanFragmentLen = 150
    
    // precomputed likelihood table
    // lookup_Table = Channel.fromPath("$baseDir/lookupTable.txt")
    
    // trees = Channel.fromPath("$baseDir/trees.txt")

    rho_rates = Channel.from(10) // For fastsimbac use this for recom rate (it doesn't accept rho)
    sample_sizes = Channel.from(10)
    genome_sizes = Channel.from(10000)
    
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

    BWA_MEM(ISOLATE_GENOME.out.firstGenome_fa, ART_ILLUMINA.out.art_fastSimBac_fq, RATE_SELECTOR.out.path_fn_modifier)

    PROCESS_SORT_INDEX(BWA_MEM.out.aligned_bam, RATE_SELECTOR.out.path_fn_modifier)

    LOFREQ(ISOLATE_GENOME.out.firstGenome_fa, PROCESS_SORT_INDEX.out.processed_bam, PROCESS_SORT_INDEX.out.processed_index, RATE_SELECTOR.out.path_fn_modifier)

    PILEUP_PAIRWISE_TABLE(LOFREQ.out.lofreqOut_vcf,PROCESS_SORT_INDEX.out.processed_bam,PROCESS_SORT_INDEX.out.processed_index, RATE_SELECTOR.out.path_fn_modifier)

    PAIRWISE_RESAMPLE(PILEUP_PAIRWISE_TABLE.out.pairwise_table_csv, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.path_fn_modifier)

    PAIRWISE_BIALLELIC_TABLE(PAIRWISE_RESAMPLE.out.pairwise_resampled_csv, RATE_SELECTOR.out.path_fn_modifier)

    LOOKUP_TABLE_LDPOP(RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.path_fn_modifier)

    PAIRWISE_LOOKUP_FORMAT(PAIRWISE_BIALLELIC_TABLE.out.pairwise_biallelic_table_csv, RATE_SELECTOR.out.path_fn_modifier)

    CUSTOM_HAP_SETS_AND_MERGE(LOOKUP_TABLE_LDPOP.out.lookupTable_txt, PAIRWISE_BIALLELIC_TABLE.out.pairwise_biallelic_table_csv, PAIRWISE_LOOKUP_FORMAT.out.lookup_format_csv, RATE_SELECTOR.out.path_fn_modifier)

    WATTERSON_ESTIMATE(LOFREQ.out.lofreqOut_vcf, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size)

    P_IJ_GRID(CUSTOM_HAP_SETS_AND_MERGE.out.eq3_csv, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    PAIRWISE_ESTIMATOR(CUSTOM_HAP_SETS_AND_MERGE.out.eq3_csv, CUSTOM_HAP_SETS_AND_MERGE.out.table_ids_for_eq3_csv, P_IJ_GRID.out.p_ij_grid_csv, LOOKUP_TABLE_LDPOP.out.lookupTable_txt, RATE_SELECTOR.out.path_fn_modifier)

    FINAL_RESULTS(PAIRWISE_ESTIMATOR.out.collected_likelihoods_csv, WATTERSON_ESTIMATE.out.theta, RATE_SELECTOR.out.path_fn_modifier)

    PROCESS_OUTPUT(FINAL_RESULTS.out.final_results_txt, RATE_SELECTOR.out.p_val, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    collectedFile = PROCESS_OUTPUT.out.processed_results_csv.collectFile(name:"collected_results.csv",storeDir:"Output/Results", keepHeader:true)

    PLOT_RESULTS(collectedFile)

}