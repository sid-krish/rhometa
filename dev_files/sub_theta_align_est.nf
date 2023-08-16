#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

// Params
// align_reads
params.fq = "/shared/homes/11849395/rhometaTESTPHB/*{1,2}.fastq.gz" // e.g. with quotes "*{1,2}.fq" for paired end
params.fa = "./contigs/*.fna"

// rho_est
params.seed = 123 // used for samtools subsamping and final bootstrap algorithm
params.filename_prefix = "none"
params.tract_len = 1000
params.lookup_grid = "0,0.25,25,0.5,50,1,100" // The range of rho values used to generate lookup tables
params.window_size = 1000 // For single end this is the read size, for paired end this is the max insert length (1000bp is a practical upper limit)
params.depth_range = "3,10" // min_depth, max_depth
params.lookup_tables = "/shared/homes/11849395/Lookup_tables/lk_multi_stp_0-100_depth_300/Lookup_tables/"

// align_reads and rho_est
params.single_end = false

// Import modules
include { FILTER_FASTQ_SINGLE_END;
        FILTER_FASTQ_PAIRED_END;
        BWA_MEM_SINGLE_END;
        BWA_MEM_PAIRED_END;
        FILTER_BAM;
        SAMTOOLS_COVERAGE
        } from './align_reads'


include {
        FILENAME_PREFIX as TE_FILENAME_PREFIX;
        SORT_BAM as TE_SORT_BAM;
        MAKE_PILEUP as TE_MAKE_PILEUP;
        SUBSAMPLE as TE_SUBSAMPLE;
        FREEBAYES as TE_FREEBAYES;
        // freebayes returns two channels, we just need the first
        THETA_ESTIMATE
        } from './sub_theta_est' addParams(output_dir: "Theta_Est_Output") // output_dir is only called internally within module, not called here


// include { FILENAME_PREFIX as RE_FILENAME_PREFIX;
//         SORT_BAM as RE_SORT_BAM;
//         MAKE_PILEUP as RE_MAKE_PILEUP;
//         SUBSAMPLE as RE_SUBSAMPLE;
//         FREEBAYES as RE_FREEBAYES;
//         PAIRWISE_TABLE_SINGLE_END as RE_PAIRWISE_TABLE_SINGLE_END;
//         PAIRWISE_TABLE_PAIRED_END as RE_PAIRWISE_TABLE_PAIRED_END;
//         RHO_ESTIMATE;
//         RESULTS_PLOT as RE_RESULTS_PLOT
//         } from './rho_est' addParams(output_dir: "Rho_Est_Output") // output_dir is only called internally within module, not called here


workflow align_reads {
    main:
        if (params.single_end == true) {
            // Channels
            fastqs_channel = Channel.fromPath(params.fq)
            reference_genome_channel = Channel.fromPath(params.fa)
            combined_inputs = fastqs_channel.combine(reference_genome_channel)

            // Process execution
            FILTER_FASTQ_SINGLE_END(combined_inputs)

            BWA_MEM_SINGLE_END(FILTER_FASTQ_SINGLE_END.out[0])

            FILTER_BAM(BWA_MEM_SINGLE_END.out)

            SAMTOOLS_COVERAGE(FILTER_BAM.out[0])
    }

        else if (params.single_end == false) {
            // Channels
            fastqs_channel = Channel.fromFilePairs(params.fq)
            reference_genome_channel = Channel.fromPath(params.fa)
            combined_inputs = fastqs_channel.combine(reference_genome_channel)

            // Process execution
            FILTER_FASTQ_PAIRED_END(combined_inputs)

            BWA_MEM_PAIRED_END(FILTER_FASTQ_PAIRED_END.out[0])

            FILTER_BAM(BWA_MEM_PAIRED_END.out)

            SAMTOOLS_COVERAGE(FILTER_BAM.out[0])
    }

    emit:
        FILTER_BAM.out[0]

}


workflow theta_est {
    take:
        theta_est_input

    main:
        TE_FILENAME_PREFIX(theta_est_input,
                        params.filename_prefix)

        TE_SORT_BAM(TE_FILENAME_PREFIX.out)

        TE_MAKE_PILEUP(TE_SORT_BAM.out)
        
        TE_SUBSAMPLE(TE_MAKE_PILEUP.out,
                    params.seed, 
                    params.depth_range)

        TE_FREEBAYES(TE_SUBSAMPLE.out[0])

        // freebayes returns two channels, we just need the first
        THETA_ESTIMATE(TE_FREEBAYES.out[0])

}


workflow rho_est {
    // Channels
    downsampled_lookup_tables = Channel.fromPath( "${params.lookup_tables}/lk_downsampled_*.csv", checkIfExists: true ).collect()

    take:
        rho_est_input

    main:
        RE_FILENAME_PREFIX(rho_est_input, 
                        params.filename_prefix, 
                        params.seed)
        
        RE_SORT_BAM(RE_FILENAME_PREFIX.out)
        
        RE_MAKE_PILEUP(RE_SORT_BAM.out)
        
        RE_SUBSAMPLE(RE_MAKE_PILEUP.out, 
                params.depth_range)
        
        RE_FREEBAYES(RE_SUBSAMPLE.out[0]) // freebayes returns two channels, we just need the first

        if (params.single_end == true) {
            RE_PAIRWISE_TABLE_SINGLE_END(RE_FREEBAYES.out[0], 
                        params.single_end, 
                        params.window_size)

            RHO_ESTIMATE(RE_PAIRWISE_TABLE_SINGLE_END.out, 
                                downsampled_lookup_tables, 
                                params.tract_len, 
                                params.depth_range,
                                params.lookup_grid)
        }
        
        else if (params.single_end == false) {
            RE_PAIRWISE_TABLE_PAIRED_END(RE_FREEBAYES.out[0], 
                        params.single_end, 
                        params.window_size)

            RHO_ESTIMATE(RE_PAIRWISE_TABLE_PAIRED_END.out, 
                                downsampled_lookup_tables, 
                                params.tract_len, 
                                params.depth_range,
                                params.lookup_grid)
        }

        RE_RESULTS_PLOT(RHO_ESTIMATE.out)
}


workflow {
    align_reads()

    theta_est(align_reads.out)

    // rho_est(align_reads.out)
}