process buscoGetDB {
    label 'basic_tools'
    label 'smallTask'

    errorStrategy 'retry'
    maxRetries 2

    if (params.cloudProcess) { publishDir "${params.permanentCacheDir}/databases/busco/${params.busco_db}", mode: 'copy', pattern: "${params.busco_db}_odb10.tar.gz" }
    else { storeDir "${params.permanentCacheDir}/databases/busco/${params.busco_db}" }

  output:
    path("${params.busco_db}_odb10.tar.gz")

  script:
    """
    [ ! -f file_versions.tsv ] && wget --no-check-certificate https://busco-data.ezlab.org/v5/data/file_versions.tsv
    date=\$(awk 'BEGIN{FS="\\t"} \$1 == "${params.busco_db}_odb10" {print \$2}' file_versions.tsv)
    wget --no-check-certificate https://busco-data.ezlab.org/v5/data/lineages/${params.busco_db}_odb10.\${date}.tar.gz
    mv ${params.busco_db}_odb10.\${date}.tar.gz ${params.busco_db}_odb10.tar.gz
    """
    
}