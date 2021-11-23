process dammit {
    label 'dammit'

    if ( params.softlink_results ) { publishDir "${params.output}/${params.rnaseq_annotation_dir}/dammit/${params.uniref90_dir}", pattern: "*" }
    else { publishDir "${params.output}/${params.rnaseq_annotation_dir}/dammit/${params.uniref90_dir}", mode: 'copy', pattern: "*" }

  input:
    path(transcriptome_assembly)
    path(busco_db)
    path(dammit_db)
    val(tool)

  output:
    tuple path("${tool}", type: 'dir'), path('uniprot_sprot_reduced.dat')

  script:
    """
    #untar dammitDB and buscoDB and substitute old buscoDB with odb10 file
    tar -zxvf ${dammit_db}
    rm -rf dbs/busco2db/*
    mkdir -p dbs/busco2db/${params.busco_db}_odb10
    tar -zxvf ${busco_db} -C dbs/busco2db/${params.busco_db}_odb10

    #create .done file to bypass dammit download check
    touch dbs/busco2db/download_and_untar:busco2db-${params.busco_db}.done

    #actual annotation
    dammit annotate ${transcriptome_assembly} --database-dir \${PWD}/dbs --busco-group ${params.busco_db} -n dammit -o ${tool} --n_threads ${task.cpus}
    cp dbs/uniprot_sprot_reduced.dat .
    """
  }

//    #BUSCO=\$(echo ${busco_db} | awk 'BEGIN{FS="_"};{print \$1}')
