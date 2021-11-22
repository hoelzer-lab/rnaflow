process buscoGetDB {
    label 'basic_tools'
    label 'smallTask'

    //if (params.cloudProcess) { publishDir "${params.permanentCacheDir}/databases/busco/${params.busco_db}", mode: 'copy', pattern: "${params.busco_db}.tar.gz" }
    //else { storeDir "${params.permanentCacheDir}/databases/busco/${params.busco_db}" }  

  output:
    file("${params.busco_db}.tar.gz")

  script:
    """
    touch ${params.busco_db}.tar.gz # so the process wont break just now
    #wget http://busco.ezlab.org/v2/datasets/${params.busco_db}.tar.gz 
    """
}

/*
putting the database name into the channel for busco later on
*/
