/************************************************************************
* FastQC
************************************************************************/
process fastqc {
    label 'fastqc'
    if (!workflow.profile.contains('node')) { label 'smallTask' }
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
