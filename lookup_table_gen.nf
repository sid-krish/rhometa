#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


def helpMessage() {

    log.info"""
    Description:

    Usage:
    nextflow run lookup_table_gen.nf [options]

    Downsample Only:
    nextflow run lookup_table_gen.nf --lk_table [str] --ldpop_rho_range [str] --lk_table_max_depth [int]

    Help:
    nextflow run lookup_table_gen.nf --help

    Options:
    --ldpop_rho_range [str], default:["0,0.01,1,1,100"], ["num_rh,max_rh"] The grid of rho values used to generate lookup tables for using the ldpop algorithm.
                                                        ldpop help: The grid has num_rh uniformly spaced points from 0 to max_rh, inclusive. (((Alternatively, to create 
                                                        a non-uniform grid, use r0,step0,r1,step1,r2,...rK. This creates a grid {r0,r0+step0,r0+2*step0,...,r1,r1+step1,...,rK}
                                                        similar to ldhelmet. Note that non-uniform grid is incompatible with vanilla ldhat.)))
    --lk_table_max_depth [int], default:[100], The max depth to generate lookup tables for
    --lk_table [str], Provide lookup table to run downsample step only
    --theta [int], default:[0.01], Population mutation rate, can be estimated value from theta_est.nf or a different value

    """.stripIndent()

}


process LDPOP_TABLE_GEN {
    publishDir "Lookup_tables", mode: "copy"
    
    maxForks 1

    input:
        val lk_table_max_depth
        val theta
        val ldpop_rho_range

    output:
        path "lookup_table.txt", emit: lookup_table_txt

    script:
    // There are other parameters that can be adjusted, I've left them out for the time being
    // also they mention twice muation and recom rate, for the mutation and recom parameters which I am unsure how to interpret
    """
    ldtable.py --cores $task.cpus -n ${lk_table_max_depth} -th ${theta} -rh ${ldpop_rho_range} --approx > lookup_table.txt
    """
}


process DOWNSAMPLE_LOOKUP_TABLE {
    publishDir "Lookup_tables", mode: "copy"

    input:
        path lookup_table
        val ldpop_rho_range
        each downsample_val

    output:
        path "lk_downsampled_${downsample_val}.csv"

    script:
    """
    m_downsample_lk_table.py ${lookup_table} ${ldpop_rho_range} ${downsample_val}
    """
}

workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // Params
    params.help = false
    params.lk_table = 'none'
    params.theta = 0.01 // Theta can be based on estimate or as desired
    params.ldpop_rho_range = "0,0.01,1,1,100"
    params.lk_table_max_depth = 200

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
        LDPOP_TABLE_GEN(params.lk_table_max_depth, params.theta, params.ldpop_rho_range)

        DOWNSAMPLE_LOOKUP_TABLE(LDPOP_TABLE_GEN.out.lookup_table_txt, params.ldpop_rho_range, depth_range)
    }

    else {
        lk_table_file = Channel.fromPath(params.lk_table)
        DOWNSAMPLE_LOOKUP_TABLE(lk_table_file, params.ldpop_rho_range, depth_range)
    }

    

}