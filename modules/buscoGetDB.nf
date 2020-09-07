process buscoGetDB {
    label 'busco'
    label 'smallTask'

    if (params.cloudProcess) { publishDir "${params.permanentCacheDir}/databases/busco/${params.busco}", mode: 'copy', pattern: "${params.busco}.tar.gz" }
    else { storeDir "${params.permanentCacheDir}/databases/busco/${params.busco}" }  

  output:
    file("${params.busco}.tar.gz")

  script:
    """
    wget http://busco.ezlab.org/v2/datasets/${params.busco}.tar.gz 
    """
}

/*
putting the database name into the channel for busco later on
*/
