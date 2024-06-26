#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


def helpMessage() {

    log.info"""
    Usage:
    nextflow run lookup_table_gen.nf [options]

    Downsample Only:
    nextflow run lookup_table_gen.nf --lk_table [str] --lookup_grid [str] --lk_table_max_depth [int]

    Help:
    nextflow run lookup_table_gen.nf --help

    Options:
    --lookup_grid [str], default:["101,100"], ["num_rh,max_rh"] The grid of rho values used to generate lookup tables for using the ldpop algorithm.
                                                   ldpop help: The grid has num_rh uniformly spaced points from 0 to max_rh, inclusive. (((Alternatively, to create 
                                                   a non-uniform grid, use r0,step0,r1,step1,r2,...rK. This creates a grid {r0,r0+step0,r0+2*step0,...,r1,r1+step1,...,rK}
                                                   similar to ldhelmet. Note that non-uniform grid is incompatible with vanilla ldhat.)))
    --lk_table_max_depth [int], default:[85], The max depth to generate lookup tables for
    --lk_table [str], Provide lookup table to run downsample step only
    --theta [float], default:[0.01], Population mutation rate per site, can be estimated value from theta_est.nf or a different value
    --output_dir [str], default:[Lookup_tables], Directory to save results in

    """.stripIndent()

}


process LDPOP_TABLE_GEN {
    publishDir params.output_dir, mode: 'copy', saveAs: {filename -> "ldpop_table_gen/${filename}"}
    
    maxForks 1

    label 'LDPOP_TABLE_GEN'

    input:
        val lk_table_max_depth
        val theta
        val lookup_grid

    output:
        path "lookup_table.txt", emit: lookup_table_txt

    script:
    // There are other parameters that can be adjusted, I've left them out for the time being
    // also they mention twice muation and recom rate, for the mutation and recom parameters which I am unsure how to interpret
    """
    ldtable.py --cores $task.cpus -n ${lk_table_max_depth} -th ${theta} -rh ${lookup_grid} --approx > lookup_table.txt
    """
}


process DOWNSAMPLE_LOOKUP_TABLE {
    publishDir params.output_dir, mode: "copy"

    label 'DOWNSAMPLE_LOOKUP_TABLE'

    input:
        path lookup_table
        val lookup_grid
        each downsample_val

    output:
        path "lk_downsampled_${downsample_val}.csv"

    script:
    """
    downsample_lk_table.py ${lookup_table} ${lookup_grid} ${downsample_val}
    """
}

// process TABLE_STORE {
//     publishDir params.output_dir, mode: 'copy'

//     input:
//         path downsampled_table
//         val theta
//         val lookup_grid
//         val store_name

//     output:
//         path "${store_name}"

//     script:
//     """
//     m_hdf5.py --table-format ${table_fmt} --init -t ${theta} -r ${lookup_grid} . ${store_name}
//     """
// }

workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // Params
    params.help = false
    params.lk_table = 'none'
    params.theta = 0.01 // Theta can be based on estimate or as desired
    params.lookup_grid = "101,100" // The grid of rho values used to generate lookup tables for using the ldpop algorithm.
    params.lk_table_max_depth = 85
    params.output_dir = 'Lookup_tables'

    depth_range = Channel.of(3 .. params.lk_table_max_depth) // 3 to max_depth val

    // Input verification
    if (params.help) {
        // Show help message from helpMessage() function
        // params.help = false by default
        helpMessage()
        exit 0
    }

    // Process execution
    if (params.lk_table == 'none') {
        LDPOP_TABLE_GEN(params.lk_table_max_depth, params.theta, params.lookup_grid)
        DOWNSAMPLE_LOOKUP_TABLE(LDPOP_TABLE_GEN.out.lookup_table_txt, params.lookup_grid, depth_range)
    }

    else {
        lk_table_file = Channel.fromPath(params.lk_table)
        DOWNSAMPLE_LOOKUP_TABLE(lk_table_file, params.lookup_grid, depth_range)
    }
}
