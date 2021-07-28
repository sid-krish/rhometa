#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


def helpMessage() {

    log.info"""
    Description:

    Usage:
    nextflow run lookup_table_gen.nf [options]

    Help:
    nextflow run lookup_table_gen.nf --help

    Options:
    --ldpop_rho_range [int,int], default:[101,100], The range of rho values used to generate lookup tables
    --num_cores [int], default:[4], The max number of cores the pipeline should use
    --lk_table_max_depth [int], default:[100], The max depth to generate lookup tables for
    --theta [int], default:[0.01], Population mutation rate, can be estimated value from theta_est.nf or a different value

    """.stripIndent()

}


process LDPOP_TABLE_GEN {
    publishDir "Lookup_tables", mode: "copy"
    
    maxForks 1

    cpus {num_cores}
    memory { 50.GB * task.attempt }
    time { 10.hours * task.attempt }
    scratch true

    errorStrategy 'retry'
    maxRetries 5

    input:
        val num_cores
        val lk_table_max_depth
        val theta
        val ldpop_rho_range

    output:
        path "lookup_table.txt", emit: lookup_table_txt

    script:
    // There are other parameters that can be adjusted, I've left them out for the time being
    // also they mention twice muation and recom rate, for the mutation and recom parameters which I am unsure how to interpret
    """
    ldtable.py --cores ${num_cores} -n ${lk_table_max_depth} -th ${theta} -rh ${ldpop_rho_range} --approx > lookup_table.txt
    """
}


process DOWNSAMPLE_LOOKUP_TABLE {
    publishDir "Lookup_tables", mode: "copy"

    maxForks 1

    cpus {num_cores}
    memory { 10.GB * task.attempt }
    time { 1.h * task.attempt }
    scratch true

    errorStrategy 'retry'
    maxRetries 5

    input:
        val num_cores
        path lookup_table
        val ldpop_rho_range
        val lk_table_max_depth

    output:
        path "lk_downsampled_*.csv"

    script:
    """
    m_downsample_lk_tables.py ${lookup_table} ${ldpop_rho_range} ${lk_table_max_depth} ${num_cores}
    """
}

workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // Params
    params.help = false
    params.theta = 0.01 // Theta can be based on estimate or as desired
    params.ldpop_rho_range = "101,100"
    params.lk_table_max_depth = 10

    params.num_cores = 4

    // Input verification
    if (params.help) {
        // Show help message from helpMessage() function
        // params.help = false by default
        helpMessage()
        exit 0
    }

    // Process execution
    LDPOP_TABLE_GEN(params.num_cores, params.lk_table_max_depth, params.theta, params.ldpop_rho_range)

    DOWNSAMPLE_LOOKUP_TABLE(params.num_cores, LDPOP_TABLE_GEN.out.lookup_table_txt, params.ldpop_rho_range, params.lk_table_max_depth)

}