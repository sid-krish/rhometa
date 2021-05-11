#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process RATE_SELECTOR {
    // This process prevents the need to use each in every process, which can be confusing
    // Perhaps this could be handled in a more elegant way using some DSL2 technique
    
    maxForks 1 // Run sequentially

    // cpus 1
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
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    memory '128 MB'

    executor 'local'
    time '15m'
    scratch true
    // queue 'i3q'

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

    cpus 1
    memory '128 MB'

    executor 'local'
    time '1h'
    scratch true
    // queue 'i3q'
    
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

    cpus 1
    memory '128 MB'

    executor 'local'
    time '1h'
    scratch true
    // queue 'i3q'
    
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

    cpus 1
    memory '128 MB'

    executor 'local'
    time '10m'
    scratch true
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
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    memory '128 MB'

    executor 'local'
    time '10m'
    scratch true
    // queue 'i3q'

    input:
        path cleanTrees
        val genome_size
        val path_fn_modifier

    output:
        path "seqgenOut.fa", emit: seqgenout_fa

    script:
    // 1 partition per tree
    // program crashes if seq length is not as the one set for fastsimbac
    """
    numTrees=\$(wc -l cleanTrees.txt | awk '{ print \$1 }')
    seq-gen -m HKY -t 4 -l ${genome_size} -z ${params.seed} -s 0.01 -p \$numTrees -of cleanTrees.txt > seqgenOut.fa
    """
}


process REFORMAT_FASTA {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    memory '128 MB'

    executor 'local'
    time '10m'
    scratch true
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
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    memory '128 MB'

    executor 'local'
    time '5m'
    scratch true
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


process ART_ILLUMINA {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    memory '128 MB'

    executor 'local'
    time '1h'
    scratch true
    // queue 'i3q'

    input:
        path reformatted_fa
        val path_fn_modifier

    output:
        path "*.fq", emit: art_out_fq // only files we are interested in
        // path "*.aln" // for testing

    script:
    // can try --rcount as an alternative to --fcov
    // number of reads/read pairs to be generated per sequence/amplicon (not be used together with -f/--fcov)
    """
    #Single end
    #art_illumina --seqSys HSXt --rndSeed ${params.seed} --noALN --quiet \
    #--in reformatted.fa --len ${params.read_len} --fcov 20 --out art_out

    #Paired end
    art_illumina --seqSys HSXt --rndSeed ${params.seed} --noALN --quiet \
    --in reformatted.fa -p --len ${params.read_len} --sdev ${params.paired_end_std_dev} \
    -m ${params.paired_end_mean_frag_len} --fcov 20 --out art_out
    """
}


process BWA_MEM {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 4
    memory '1 GB'

    executor 'local'
    time '30m'
    scratch true
    // queue 'i3q'

    input:
        path firstGenome_fa
        path art_out_fq
        val path_fn_modifier

    output:
        path "Aligned.bam", emit: aligned_bam

    script:
    // using first fa entry only (one genome)
    """
    bwa index firstGenome.fa

    #Single end
    #bwa mem -t 4 firstGenome.fa art_out.fq > Aligned.sam

    #Paired end
    bwa mem -t 4 firstGenome.fa art_out1.fq art_out2.fq > Aligned.sam

    samtools view -bS Aligned.sam > Aligned.bam
    """
}


process PROCESS_SORT_INDEX{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 4
    memory '10 GB'

    executor 'local'
    time '30m'
    scratch true
    // queue 'i3q'
    
    input:
        path aligned_bam
        val path_fn_modifier

    output:
        path "Aligned.csorted_fm_md.bam", emit: processed_bam
        path "Aligned.csorted_fm_md.bam.bai", emit: processed_index
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
    """
}


process LOFREQ{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 4
    memory '1 GB'

    executor 'local'
    time '30m'
    scratch true
    // queue 'i3q'

    input:
        path firstGenome_fa
        path bam
        val path_fn_modifier

    output:
        path "lofreqOut.vcf", emit: lofreqOut_vcf

    script:
    """
    samtools faidx firstGenome.fa
    samtools sort --threads 4 Aligned.bam -o Aligned.csorted.bam
    samtools index -@ 4 Aligned.csorted.bam

    #lofreq call -f firstGenome.fa -o lofreqOut.vcf Aligned.csorted.bam
    lofreq call-parallel --pp-threads 4 --no-default-filter -f firstGenome.fa -o lofreqOut.vcf Aligned.csorted.bam
    """
}


process PAIRWISE_TABLE_SINGLE_END{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    memory '1 GB'

    executor 'local'
    time '30m'
    scratch true
    // queue 'i3q'

    input:
        path lofreqOut_vcf 
        path bam_stats_txt
        path bam_file
        path bam_index
        val path_fn_modifier

    output:
        path "pairwise_table.csv", emit: pairwise_table_csv

    script:
    """
    #!/bin/bash
    max_read_len=\$(grep "maximum length" bam_stats.txt | cut -f 3)
    pairwise_table.py \$max_read_len ${bam_file} ${lofreqOut_vcf}
    """
}


process PAIRWISE_TABLE_PAIRED_END{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    memory '1 GB'

    executor 'local'
    time '30m'
    scratch true
    // queue 'i3q'

    input:
        path lofreqOut_vcf 
        path bam_file
        val path_fn_modifier

    output:
        path "pairwise_table.csv", emit: pairwise_table_csv

    script:
    """
    pairwise_table_paired_end.py ${bam_file} ${lofreqOut_vcf}
    """
}


process PAIRWISE_BIALLELIC_TABLE{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    // maxForks 4

    cpus 1
    memory '1 GB'

    executor 'local'
    time '30m'
    scratch true
    // queue 'i3q'

    input:
        path pairwise_table_csv
        each depth
        val path_fn_modifier

    output:
        path 'pairwise_biallelic_table_depth_*.csv', optional: true, emit: pairwise_biallelic_table_csv // after filtering some files will be empty

    script:
    // NOTE
    // There was serious error which prevented it from working properly, but is now fixed.
    // Need to fix in other projects that use this also
    """
    biallelic_filter.py pairwise_table.csv ${depth}
    """
}


process PAIRWISE_LOOKUP_FORMAT {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    // maxForks 4

    cpus 1
    memory '1 GB'

    executor 'local'
    time '30m'
    scratch true
    // queue 'i3q'

    input:
        path pairwise_biallelic_table_csv
        val path_fn_modifier

    output:
        path "lookup_format_depth_*.csv", emit: lookup_format_csv
        path pairwise_biallelic_table_csv, emit: pairwise_biallelic_table_csv // just resending file

    script:
    """
    depth=\$(echo ${pairwise_biallelic_table_csv} | cut -d. -f1 | cut -d_ -f5)
    pairwise_lookup_format_pyrho.py ${pairwise_biallelic_table_csv} "\$depth"
    """
}


process CUSTOM_HAP_SETS_AND_MERGE {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    // maxForks 1

    cpus 1
    memory '4 GB'

    executor 'local'
    time '1h'
    scratch true
    // queue 'i3q'

    input:
        path pairwise_biallelic_table_csv
        path lookup_format_csv
        val path_fn_modifier

    output:
        path "table_ids_for_eq3_depth_*.csv", emit: table_ids_for_eq3_csv
        path "eq3_depth_*.csv", emit: eq3_csv

    script:
    // lookup table location hard coded. Need to fix
    """
    depth=\$(echo ${pairwise_biallelic_table_csv} | cut -d. -f1 | cut -d_ -f5)
    custom_hap_sets_and_merge.py ${pairwise_biallelic_table_csv} ${lookup_format_csv} ${params.ldpop_rho_range} "\$depth" $PWD/Lookup_tables/lk_downsampled_"\$depth".csv  > table_ids_for_eq3_depth_"\$depth".csv
    """
}


process WATTERSON_ESTIMATE {
    
    maxForks 1

    cpus 1
    memory '512 MB'

    executor 'local'
    time '10m'
    scratch true
    // queue 'i3q'

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

    // maxForks 4

    cpus 1
    memory '1 GB'

    executor 'local'
    time '30m'
    scratch true
    // queue 'i3q'

    input:
        path eq3_csv
        path table_ids_for_eq3_csv
        val genome_size
        val path_fn_modifier

    output:
        path "p_ij_grid_depth_*.csv", emit: p_ij_grid_csv
        path table_ids_for_eq3_csv, emit: table_ids_for_eq3_csv // just resending file
        path eq3_csv, emit: eq3_csv // just resending file

    script:
    """
    depth=\$(echo ${eq3_csv} | cut -d. -f1 | cut -d_ -f3)
    pij_grid_vectorised.py ${genome_size} ${params.recom_tract_len} ${params.ldpop_rho_range} ${eq3_csv} "\$depth"
    """
}


process PAIRWISE_ESTIMATOR {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}
    
    // maxForks 4

    cpus 1
    memory '4 GB'

    executor 'local'
    time '30m'
    scratch true
    // queue 'i3q'
    
    input:
        path eq3_csv
        path table_ids_for_eq3_csv
        path p_ij_grid_csv
        val path_fn_modifier

    output:
        path "collected_likelihoods_depth_*.csv", emit: collected_likelihoods_csv
    
    script:
    """
    depth=\$(echo ${eq3_csv} | cut -d. -f1 | cut -d_ -f3)
    pairwise_rho_estimator_intp_rect_biv.py ${eq3_csv} ${table_ids_for_eq3_csv} ${p_ij_grid_csv} ${params.ldpop_rho_range} "\$depth" $PWD/Lookup_tables/lk_downsampled_"\$depth".csv
    """

}


process FINAL_RESULTS {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    // maxForks 4

    cpus 1
    memory '256 MB'

    executor 'local'
    time '10m'
    scratch true
    // queue 'i3q'

    input:
        path collectedFile
        val theta
        val path_fn_modifier

    output:
        path "final_results_depth_*.txt", emit: final_results_txt

    script:
    """
    depth=\$(echo ${collectedFile} | cut -d. -f1 | cut -d_ -f4)
    final_results.py ${collectedFile} "\$depth" ${theta}
    """
}


process AGGREGATE_RESULTS{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${path_fn_modifier}_${filename}"}

    maxForks 1

    cpus 1
    memory '256 MB'

    executor 'local'
    time '10m'
    scratch true
    // queue 'i3q'

    input:
        path collected_results
        val theta
        val path_fn_modifier

    output:
        path "Aggregate_results_for_full_dataset.out"
        
    script:
        """
        merge_and_get_final_result.py ${theta}
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

    params.read_len = 150

    params.paired_end_mean_frag_len = 300
    params.paired_end_std_dev = 50 // +- mean frag len
    
    depth = Channel.from(10..200)
    
    // trees = Channel.fromPath("$baseDir/trees.txt")
    // custom_pairwise_pairwise_table = Channel.fromPath("$baseDir/pairwise_table.csv") // for testing
    // custom_pairwise_pairwise_biallelic_table = Channel.fromPath("$baseDir/pairwise_biallelic_table.csv") // for testing

    params.rho_rates = 15  // Using each multiple times, causes nextflow to work improperly so removed it here for use with depth
    params.sample_sizes = 15
    params.genome_sizes = 10000
    
    RATE_SELECTOR(params.rho_rates, params.sample_sizes, params.genome_sizes)

    MS(RATE_SELECTOR.out.p_val, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    // FAST_SIM_BAC(RATE_SELECTOR.out.p_val, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    // MS_PRIME(RATE_SELECTOR.out.p_val, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    CLEAN_TREES(MS.out.trees_txt, RATE_SELECTOR.out.path_fn_modifier)

    // CLEAN_TREES(FAST_SIM_BAC.out.trees_txt, RATE_SELECTOR.out.path_fn_modifier)

    SEQ_GEN(CLEAN_TREES.out.cleanTrees_txt, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    REFORMAT_FASTA(SEQ_GEN.out.seqgenout_fa, RATE_SELECTOR.out.path_fn_modifier)

    ISOLATE_GENOME(REFORMAT_FASTA.out.reformatted_fa, RATE_SELECTOR.out.path_fn_modifier)

    ART_ILLUMINA(REFORMAT_FASTA.out.reformatted_fa, RATE_SELECTOR.out.path_fn_modifier)

    BWA_MEM(ISOLATE_GENOME.out.firstGenome_fa, ART_ILLUMINA.out.art_out_fq, RATE_SELECTOR.out.path_fn_modifier)

    // PROCESS_SORT_INDEX(BWA_MEM.out.aligned_bam, RATE_SELECTOR.out.path_fn_modifier)

    LOFREQ(ISOLATE_GENOME.out.firstGenome_fa, BWA_MEM.out.aligned_bam, RATE_SELECTOR.out.path_fn_modifier)

    // PAIRWISE_TABLE_SINGLE_END(LOFREQ.out.lofreqOut_vcf,PROCESS_SORT_INDEX.out.bam_stats_txt,PROCESS_SORT_INDEX.out.processed_bam,PROCESS_SORT_INDEX.out.processed_index, RATE_SELECTOR.out.path_fn_modifier)

    // PAIRWISE_BIALLELIC_TABLE(PAIRWISE_TABLE_SINGLE_END.out.pairwise_table_csv, depth, RATE_SELECTOR.out.path_fn_modifier)

    PAIRWISE_TABLE_PAIRED_END(LOFREQ.out.lofreqOut_vcf, BWA_MEM.out.aligned_bam, RATE_SELECTOR.out.path_fn_modifier)

    PAIRWISE_BIALLELIC_TABLE(PAIRWISE_TABLE_PAIRED_END.out.pairwise_table_csv, depth, RATE_SELECTOR.out.path_fn_modifier)

    PAIRWISE_LOOKUP_FORMAT(PAIRWISE_BIALLELIC_TABLE.out.pairwise_biallelic_table_csv, RATE_SELECTOR.out.path_fn_modifier)

    CUSTOM_HAP_SETS_AND_MERGE(PAIRWISE_LOOKUP_FORMAT.out.pairwise_biallelic_table_csv, PAIRWISE_LOOKUP_FORMAT.out.lookup_format_csv, RATE_SELECTOR.out.path_fn_modifier)

    WATTERSON_ESTIMATE(LOFREQ.out.lofreqOut_vcf, RATE_SELECTOR.out.sample_size, RATE_SELECTOR.out.genome_size)

    P_IJ_GRID(CUSTOM_HAP_SETS_AND_MERGE.out.eq3_csv, CUSTOM_HAP_SETS_AND_MERGE.out.table_ids_for_eq3_csv, RATE_SELECTOR.out.genome_size, RATE_SELECTOR.out.path_fn_modifier)

    PAIRWISE_ESTIMATOR(P_IJ_GRID.out.eq3_csv, P_IJ_GRID.out.table_ids_for_eq3_csv, P_IJ_GRID.out.p_ij_grid_csv, RATE_SELECTOR.out.path_fn_modifier)

    FINAL_RESULTS(PAIRWISE_ESTIMATOR.out.collected_likelihoods_csv, WATTERSON_ESTIMATE.out.theta, RATE_SELECTOR.out.path_fn_modifier)

    collected_results = FINAL_RESULTS.out.final_results_txt.collect()

    AGGREGATE_RESULTS(collected_results, WATTERSON_ESTIMATE.out.theta, RATE_SELECTOR.out.path_fn_modifier)

}