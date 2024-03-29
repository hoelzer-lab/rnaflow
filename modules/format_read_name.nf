process format_read_name {
    label 'basic_tools'
    tag "$meta.sample"
    if (!workflow.profile.contains('node')) { label 'smallTask' }

    input:
    tuple val(meta), path(reads) 

    output:
    tuple val(meta), path("${meta.sample}*.raw.fastq.gz") 

    script:
    """
    mv *fastq.gz ${meta.sample}.raw.fastq.gz
    """
}