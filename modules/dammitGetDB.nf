process dammitGetDB {
    label 'basic_tools'
    //if (!params.cloudProcess) { label 'smallTask' }

    errorStrategy 'retry'
    maxRetries 2

    if (params.cloudProcess) { publishDir "${params.permanentCacheDir}/databases/dammit/${params.busco_db}", mode: 'copy', pattern: "dbs.tar.gz" }
    else { storeDir "${params.permanentCacheDir}/databases/dammit/${params.busco_db}/" }

  output:
    path("dbs.tar.gz")

  script:
    """
    wget https://zenodo.org/api/files/773602b2-e239-4489-8f97-cebb99f36a81/dbs.tar.gz
    """
    
}