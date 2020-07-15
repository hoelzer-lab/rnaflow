process dammitGetDB {
    label 'dammitDB'
    label 'smallTask'

    if (params.full) {
      if (params.cloudProcess) { publishDir "${params.permanentCacheDir}/databases/dammit/full/${params.busco}", mode: 'copy', pattern: "dbs.tar.gz" }
      else { storeDir "${params.permanentCacheDir}/databases/dammit/full/${params.busco}/" }  
    }
    else {
      if (params.cloudProcess) { publishDir "${params.permanentCacheDir}/databases/dammit/${params.busco}", mode: 'copy', pattern: "dbs.tar.gz" }
      else { storeDir "${params.permanentCacheDir}/databases/dammit/${params.busco}/" }  
    }

  input:
    path(busco_db)

  output:
    path("dbs.tar.gz")

  script:
    if (params.full)
    """
    BUSCO=\$(echo ${params.busco} | awk 'BEGIN{FS="_"};{print \$1}')
    dammit databases --install --database-dir \${PWD}/dbs --busco-group \${BUSCO} --full
    # the busco download fails so use the busco db we anyway already downloaded
    tar -zxvf ${busco_db} -C dbs/busco2db/
    # in addition, download metadata from uniprot/swissprot
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
    gunzip uniprot_sprot.dat.gz
    awk '{if(\$1=="ID" || \$1=="AC" || \$1=="DE" || \$1=="OS" || \$1=="OC" || \$2=="GO;"){print \$0}}' uniprot_sprot.dat > dbs/uniprot_sprot_reduced.dat
    rm uniprot_sprot.dat
    tar zcvf dbs.tar.gz dbs
    rm -rf dbs
    """
  else
    """
    BUSCO=\$(echo ${params.busco} | awk 'BEGIN{FS="_"};{print \$1}')
    dammit databases --install --database-dir \${PWD}/dbs --busco-group \${BUSCO}
    # the busco download fails so use the busco db we anyway already downloaded
    tar -zxvf ${busco_db} -C dbs/busco2db/
    # in addition, download metadata from uniprot/swissprot
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
    gunzip uniprot_sprot.dat.gz
    awk '{if(\$1=="ID" || \$1=="AC" || \$1=="DE" || \$1=="OS" || \$1=="OC" || \$2=="GO;"){print \$0}}' uniprot_sprot.dat > dbs/uniprot_sprot_reduced.dat
    rm uniprot_sprot.dat
    tar zcvf dbs.tar.gz dbs
    rm -rf dbs
    """
}

