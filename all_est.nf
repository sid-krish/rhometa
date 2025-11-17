#! /usr/bin/env nextflow
nextflow.enable.dsl=2

// Params
// Input files
params.input_csv = 'none' // csv file with paths to bams and references, "bam" and "reference" need to the header
params.bam = 'none'
params.fa = 'none'

// rho_est
params.seed = [123] // used for samtools subsamping
params.filename_prefix = "none"
params.tract_len = 1000
params.lookup_grid = "101,100" // The range of rho values used to generate lookup tables
params.window_size = 1000 // For single end this is the read size, for paired end this is the max insert length (1000bp is a practical upper limit)
params.depth_range = "3,85" // min_depth, max_depth
params.lookup_tables = "Lookup_tables"
params.vcf_filtered_file = 'freebayes_filt.vcf'

// align_reads and rho_est
params.single_end = false

// VCF filter settings
params.snp_qual = 20 // Minimum phred-scaled quality score to filter vcf by
params.min_snp_depth = 10 // Minimum read depth to filter vcf by
params.top_depth_cutoff_percentage = 5 // Top n percent of depth to cutoff from the vcf file

// Align reads output directory
params.align_out_dir = "Align_Reads_Output"

// Rho and theta estimate output directories
params.rho_est_out_dir = "Rho_Est_Output"
params.theta_est_out_dir = "Theta_Est_Output"

// Final results output directory
params.final_results_out_dir = "Final_Output"


// Import modules
include {
    SAMTOOLS_COVERAGE
} from './align_reads'


include {
    FILENAME_PREFIX as TE_FILENAME_PREFIX;
    SORT_BAM as TE_SORT_BAM;
    ANGSD_THETA_ESTIMATE
} from './theta_est'


include {
    FILENAME_PREFIX as RE_FILENAME_PREFIX;
    SORT_BAM as RE_SORT_BAM;
    GET_READ_DEPTH as RE_GET_READ_DEPTH;
    SUBSAMPLE as RE_SUBSAMPLE;
    FREEBAYES as RE_FREEBAYES;
    VCF_FILTER as RE_VCF_FILTER;
    PAIRWISE_TABLE_SINGLE_END as RE_PAIRWISE_TABLE_SINGLE_END;
    PAIRWISE_TABLE_PAIRED_END as RE_PAIRWISE_TABLE_PAIRED_END;
    RHO_ESTIMATE;
    RESULTS_PLOT as RE_RESULTS_PLOT
} from './rho_est'


include {
    COLLECT_RESULTS;
    COMPUTE_RATIOS
} from './collect_results_compute_ratios'


workflow align_reads {
    take:
        bam_and_fa

    main:
        SAMTOOLS_COVERAGE(bam_and_fa)
}


workflow theta_est {
    take:
        bam_and_fa

    main:
        TE_FILENAME_PREFIX(bam_and_fa, params.filename_prefix)

        TE_SORT_BAM(TE_FILENAME_PREFIX.out)

        ANGSD_THETA_ESTIMATE(TE_SORT_BAM.out)

    emit:
        ANGSD_THETA_ESTIMATE.out[0].map { it[4] } // emit the angsd_theta_final.csv file
}


workflow rho_est {

    take:
        bam_and_fa
        downsampled_lookup_tables
    
    main:
        RE_FILENAME_PREFIX(bam_and_fa,
            params.filename_prefix,
            params.seed)

        RE_SORT_BAM(RE_FILENAME_PREFIX.out)

        RE_GET_READ_DEPTH(RE_SORT_BAM.out)

        RE_SUBSAMPLE(RE_GET_READ_DEPTH.out,
            params.depth_range)

        RE_FREEBAYES(RE_SUBSAMPLE.out[0])

        RE_VCF_FILTER(RE_FREEBAYES.out, params.snp_qual, params.min_snp_depth, params.top_depth_cutoff_percentage)

        if (params.single_end == true) {
            RE_PAIRWISE_TABLE_SINGLE_END(RE_VCF_FILTER.out,
                params.single_end,
                params.window_size)

            RHO_ESTIMATE(RE_PAIRWISE_TABLE_SINGLE_END.out,
                downsampled_lookup_tables,
                params.tract_len,
                params.depth_range,
                params.lookup_grid)
        }

        else if (params.single_end == false) {
            RE_PAIRWISE_TABLE_PAIRED_END(RE_VCF_FILTER.out,
                params.single_end,
                params.window_size)

            RHO_ESTIMATE(RE_PAIRWISE_TABLE_PAIRED_END.out,
                downsampled_lookup_tables,
                params.tract_len,
                params.depth_range,
                params.lookup_grid)
        }

    // RE_RESULTS_PLOT(RHO_ESTIMATE.out)

    emit:
        RHO_ESTIMATE.out[0].map { it[2] } // emit the rho_estimate.csv file
}


workflow collect_results_compute_ratios {
    take:
        theta_ch
        rho_ch

    main:
        COLLECT_RESULTS(theta_ch, rho_ch)
        COMPUTE_RATIOS(COLLECT_RESULTS.out)
}


workflow {
    main:
         // Channels
        downsampled_lookup_tables = Channel.fromPath( "${params.lookup_tables}/lk_downsampled_*.csv", checkIfExists: true ).collect()

        if (params.input_csv == 'none') {
            
            if (params.bam == 'none' || params.fa == 'none') {
            println "Error: --bam and --fa parameters are required (or use --input_csv). Use --help for more information."
            exit 1
            }

            else {
            bam_channel = Channel.fromPath( params.bam, checkIfExists: true )
            fa_channel = Channel.fromPath( params.fa, checkIfExists: true )
            bam_and_fa = bam_channel.combine(fa_channel)
            }
        }

        else if (params.input_csv != 'none') {
            println "Using input csv file: ${params.input_csv}"
            bam_and_fa = Channel.fromPath( params.input_csv, checkIfExists: true ).splitCsv(header:true).map { row -> tuple( file(row.bam), file(row.reference) ) }
        }

        align_reads(bam_and_fa)
        theta_est(bam_and_fa)
        rho_est(bam_and_fa, downsampled_lookup_tables)
        collect_results_compute_ratios(theta_est.out.collect(), rho_est.out.collect())
}