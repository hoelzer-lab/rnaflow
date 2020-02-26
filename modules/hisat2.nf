/************************************************************************
* HISAT2 INDEX
************************************************************************/
process hisat2index {
    label 'hisat2'

    if (params.cloudProcess) { publishDir "${params.cloudDatabase}/genomes/${params.species}", mode: 'copy', pattern: "${params.species}*.ht2" }
    else { storeDir "nextflow-autodownload-databases/genomes/${params.species}" }  

    input:
    file(reference)

    output:
    tuple file(reference), file("${params.species}*.ht2")

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
    publishDir "${params.output}/${params.hisat2_dir}", mode: 'copy', pattern: "${sample_name}.sorted.bam"

    input:
    tuple val(sample_name), file(reads)
    tuple file(reference), file(index)

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


