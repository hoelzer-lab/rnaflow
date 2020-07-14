process dammit {
    //label 'dammitDB'
    container 'nanozoo/dammit:1.2--b47259e'
    publishDir "${params.output}/${params.annotation_dir}/dammit", mode: 'copy', pattern: "${tool}"

  input:
    tuple val(name), path(transcriptome_assembly)
    path(dbs)
    val(tool)

  output:
    tuple val(name), path("${tool}", type: 'dir'), path('uniprot_sprot_reduced.dat')

  script:
    """
    tar zxvf ${dbs}
    BUSCO=\$(echo ${params.busco} | awk 'BEGIN{FS="_"};{print \$1}')
    dammit annotate ${transcriptome_assembly} --database-dir \${PWD}/dbs --busco-group \${BUSCO} -n ${name} -o ${tool} --n_threads ${task.cpus} #--full 
    cp dbs/uniprot_sprot_reduced.dat .
    rm -rf dbs
    """
  }

