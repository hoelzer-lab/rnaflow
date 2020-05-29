/************************************************************************
* FastQC
************************************************************************/
process fastqc {
    label 'fastqc'
    label 'smallTask'

    input:
    tuple val(name), path(reads) 

    output:
    tuple val(name), path("*.tar.gz"), emit: tar

    script:
    """
    mkdir -p fastqc_\$(basename \$PWD)
    fastqc -t ${task.cpus} ${reads} -o fastqc_\$(basename \$PWD)
    tar zcvf fastqc_\$(basename \$PWD).tar.gz fastqc_\$(basename \$PWD)
    rm -rf fastqc_\$(basename \$PWD)
    """
}
