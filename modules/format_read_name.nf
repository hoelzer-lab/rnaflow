process format_read_name {
    label 'format_read_name'
    tag "$meta.sample"
    label 'smallTask'

    input:
    tuple val(meta), path(reads) 

    output:
    tuple val(meta), path("${meta.sample}*.raw.fastq.gz") 

    script:
    """
    mv *fastq.gz ${meta.sample}.raw.fastq.gz
    """
}