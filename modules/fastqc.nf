/************************************************************************
* FastQC
************************************************************************/
process fastqc {
    label 'fastqc'
    if (!params.workflow.contains('node')) { label 'smallTask' }
    tag "$meta.sample"
 
    input:
    tuple val(meta), path(reads) 

    output:
    path("*_fastqc.zip", emit: zip)

    script:
    """
    fastqc --noextract -t ${task.cpus} ${reads}
    """
}
