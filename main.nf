#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

process LOFREQ{
    publishDir "Output", mode: "copy"

    maxForks 1

    // cpus 4
    // memory '5 GB'

    // executor 'pbspro'
    // time '30m'
    // scratch true
    // queue 'i3q'

    input:
        path reference_genome_file
        path bam_file
        
    output:
        path "lofreqOut.vcf", emit: lofreqOut_vcf

    script:
    """
    samtools faidx ${reference_genome_file}
    samtools sort --threads 4 ${bam_file} -o chrom_sorted.bam
    samtools index -@ 4 chrom_sorted.bam

    #lofreq call -f ${reference_genome_file} -o lofreqOut.vcf chrom_sorted.bam
    lofreq call-parallel --pp-threads 4 --no-default-filter -f ${reference_genome_file} -o lofreqOut.vcf chrom_sorted.bam
    """
}


process PAIRWISE_TABLE_SINGLE_END{
    publishDir "Output", mode: "copy"

    maxForks 1

    // cpus 1
    // memory '5 GB'

    // executor 'pbspro'
    // time '30m'
    // scratch true
    // queue 'i3q'

    input:
        path lofreqOut_vcf 
        path bam_stats_txt
        path bam_file
        path bam_index
        

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
    publishDir "Output", mode: "copy"

    maxForks 1

    // cpus 4
    // memory '25 GB'

    // executor 'pbspro'
    // time '2d'
    // scratch true
    // queue 'i3q'

    input:
        path lofreqOut_vcf 
        path bam_file
        

    output:
        path "pairwise_table.csv", emit: pairwise_table_csv

    script:
    """
    pairwise_table_paired_end.py ${bam_file} ${lofreqOut_vcf}
    """
}


process PAIRWISE_BIALLELIC_TABLE{
    publishDir "Output", mode: "copy"

    // maxForks 4

    // cpus 1
    // memory '10 GB'

    // executor 'pbspro'
    // time '5h'
    // scratch true
    // queue 'i3q'

    input:
        path pairwise_table_csv
        each depth
        

    output:
        path 'pairwise_biallelic_table_depth_*.csv', optional: true, emit: pairwise_biallelic_table_csv // after filtering some files will be empty

    script:
    """
    biallelic_filter.py pairwise_table.csv ${depth}
    """
}


process PAIRWISE_LOOKUP_FORMAT {
    publishDir "Output", mode: "copy"

    // maxForks 4

    // cpus 1
    // memory '10 GB'

    // executor 'pbspro'
    // time '5h'
    // scratch true
    // queue 'i3q'

    input:
        path pairwise_biallelic_table_csv
        

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
    publishDir "Output", mode: "copy"

    // maxForks 1

    // cpus 1
    // memory '30 GB'

    // executor 'pbspro'
    // time '5h'
    // scratch true
    // queue 'i3q'

    input:
        path pairwise_biallelic_table_csv
        path lookup_format_csv
        

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


process P_IJ_GRID {
    publishDir "Output", mode: "copy"

    // maxForks 4

    // cpus 1
    // memory '10 GB'

    // executor 'pbspro'
    // time '5h'
    // scratch true
    // queue 'i3q'

    input:
        path eq3_csv
        path table_ids_for_eq3_csv
        val genome_size
        

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
    publishDir "Output", mode: "copy"

    // maxForks 4

    // cpus 1
    // memory '30 GB'

    // executor 'pbspro'
    // time '5h'
    // scratch true
    // queue 'i3q'
    
    input:
        path eq3_csv
        path table_ids_for_eq3_csv
        path p_ij_grid_csv
        

    output:
        path "collected_likelihoods_depth_*.csv", emit: collected_likelihoods_csv
    
    script:
    """
    depth=\$(echo ${eq3_csv} | cut -d. -f1 | cut -d_ -f3)
    pairwise_rho_estimator_intp_rect_biv.py ${eq3_csv} ${table_ids_for_eq3_csv} ${p_ij_grid_csv} ${params.ldpop_rho_range} "\$depth" $PWD/Lookup_tables/lk_downsampled_"\$depth".csv
    """

}


process FINAL_RESULTS {
    publishDir "Output", mode: "copy"

    // maxForks 4

    // cpus 1
    // memory '1 GB'

    // executor 'pbspro'
    // time '1h'
    // scratch true
    // queue 'i3q'

    input:
        path collectedFile
        

    output:
        path "final_results_depth_*.txt", emit: final_results_txt

    script:
    """
    depth=\$(echo ${collectedFile} | cut -d. -f1 | cut -d_ -f4)
    final_results.py ${collectedFile} "\$depth" ${params.theta}
    """
}


process AGGREGATE_RESULTS{
    publishDir "Output", mode: "copy"

    maxForks 1

    // cpus 1
    // memory '1 GB'

    // executor 'pbspro'
    // time '1h'
    // scratch true
    // queue 'i3q'

    input:
        path collected_results
        

    output:
        path "Aggregate_results_for_full_dataset.out"
        
    script:
        """
        merge_and_get_final_result.py ${params.theta}
        """

}

// parameters

// params.bam_file = ""

// input validation
// if (params.bam_file) {

// } else {
//     exit 1, 'Missing bam file'
// }


workflow {

    params.bam_file = ""
    params.reference_genome = ""
    params.theta = 0.005

    bam_file_channel = Channel.fromPath( params.bam_file )
    reference_genome_channel = Channel.fromPath( params.reference_genome )

    params.genome_size = 1911230 // Needs to be programatically determied
    params.recom_tract_len = 500 // Needs to be programatically determied

    params.ldpop_rho_range = "101,1"
    
    depth = Channel.from(2..200) // Needs to be determined based on min and max depth and used in the pipeline accordingly, should be $min..$max, for example

    // Feature idea: automatic single end / paired end read detection

    LOFREQ(reference_genome_channel, bam_file_channel)

    // PAIRWISE_TABLE_SINGLE_END(LOFREQ.out.lofreqOut_vcf,PROCESS_SORT_INDEX.out.bam_stats_txt,PROCESS_SORT_INDEX.out.processed_bam,PROCESS_SORT_INDEX.out.processed_index)

    // PAIRWISE_BIALLELIC_TABLE(PAIRWISE_TABLE_SINGLE_END.out.pairwise_table_csv, depth)

    PAIRWISE_TABLE_PAIRED_END(LOFREQ.out.lofreqOut_vcf, bam_file_channel)

    PAIRWISE_BIALLELIC_TABLE(PAIRWISE_TABLE_PAIRED_END.out.pairwise_table_csv, depth)

    PAIRWISE_LOOKUP_FORMAT(PAIRWISE_BIALLELIC_TABLE.out.pairwise_biallelic_table_csv)

    CUSTOM_HAP_SETS_AND_MERGE(PAIRWISE_LOOKUP_FORMAT.out.pairwise_biallelic_table_csv, PAIRWISE_LOOKUP_FORMAT.out.lookup_format_csv)

    P_IJ_GRID(CUSTOM_HAP_SETS_AND_MERGE.out.eq3_csv, CUSTOM_HAP_SETS_AND_MERGE.out.table_ids_for_eq3_csv, params.genome_size)

    PAIRWISE_ESTIMATOR(P_IJ_GRID.out.eq3_csv, P_IJ_GRID.out.table_ids_for_eq3_csv, P_IJ_GRID.out.p_ij_grid_csv)

    FINAL_RESULTS(PAIRWISE_ESTIMATOR.out.collected_likelihoods_csv)

    collected_results = FINAL_RESULTS.out.final_results_txt.collect()

    AGGREGATE_RESULTS(collected_results)

}