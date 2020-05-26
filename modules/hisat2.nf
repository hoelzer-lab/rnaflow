/************************************************************************
* HISAT2 INDEX
************************************************************************/
process hisat2index {
    label 'hisat2'
    input:
    path(reference)

    output:
    tuple path(reference), path("${reference.baseName}*.ht2")

    script:
    """
    hisat2-build -p ${task.cpus} ${reference} ${reference.baseName}
    """
}


/************************************************************************
* HISAT2
************************************************************************/
process hisat2 {
    label 'hisat2'
    memory { task.memory * task.attempt }

    if (params.cloudProcess) { publishDir "${params.output}/${params.hisat2_dir}", mode: 'copy', pattern: "*.sorted.bam" }
    else { publishDir "${params.output}/${params.hisat2_dir}", pattern: "*.sorted.bam" }

    input:
    tuple val(sample_name), path(reads)
    tuple path(reference), path(index)

    output:
    tuple val(sample_name), path("${sample_name}.sorted.bam"), emit: sample_bam 
    path "${sample_name}_summary.log", emit: log

    script:
    if (params.mode == 'single') {
    """
    hisat2 -x ${reference.baseName} -U ${reads[0]} -p ${task.cpus} --new-summary --summary-file ${sample_name}_summary.log | samtools view -bS | samtools sort -o ${sample_name}.sorted.bam -T tmp --threads ${task.cpus}
    """
    }
    else {
    """
    hisat2 -x ${reference.baseName} -1 ${reads[0]} -2 ${reads[1]} -p ${task.cpus} --new-summary --summary-file ${sample_name}_summary.log | samtools view -bS | samtools sort -o ${sample_name}.sorted.bam -T tmp --threads ${task.cpus}
    """
    } 
}


