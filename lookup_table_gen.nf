#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process LDPOP_TABLE_GEN {
    publishDir "Lookup_tables", mode: "copy"
    
    maxForks 1

    cpus {num_cores}

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

    // Results will be output to Lookup_tables folder in current dir

    params.num_cores = 4
   
    params.theta = 0.01 // Theta can be based on estimate or as desired
    params.ldpop_rho_range = "101,100"

    params.lk_table_max_depth = 100

    LDPOP_TABLE_GEN(params.num_cores, params.lk_table_max_depth, params.theta, params.ldpop_rho_range)

    DOWNSAMPLE_LOOKUP_TABLE(params.num_cores, LDPOP_TABLE_GEN.out.lookup_table_txt, params.ldpop_rho_range, params.lk_table_max_depth)

}