/************************************************************************
* MINIMAP2 INDEX
************************************************************************/
process minimap2index {
    label 'minimap2'
    input:
    path(reference)

    output:
    tuple path(reference), path("${reference.baseName}.mmi")

    script:
    """
    minimap2 -x map-ont -d ${reference.baseName}.mmi ${reference}
    """
}


/************************************************************************
* MINIMAP2
************************************************************************/
process minimap2 {
    label 'minimap2'
    tag "$meta.sample"
    
    if ( params.softlink_results ) { publishDir "${params.output}/${params.minimap2_dir}", pattern: "*.sorted.bam" }
    else { publishDir "${params.output}/${params.minimap2_dir}", mode: 'copy', pattern: "*.sorted.bam" }

    input:
    tuple val(meta), path(reads)
    tuple path(reference), path(index)
    val(additionalParams)

    output:
    tuple val(meta), path("${meta.sample}.sorted.bam"), emit: sample_bam 
    path "${meta.sample}_summary.log", emit: log

    script:
    """
    minimap2 -a -t ${task.cpus} ${index} ${reads[0]} 2> ${meta.sample}_summary.log | samtools view -bS | samtools sort -o ${meta.sample}.sorted.bam -T tmp --threads ${task.cpus}    
    """
}
    

process index_bam_minimap2 {
    label 'minimap2'
    label 'smallTask'    
    
    if ( params.softlink_results ) { publishDir "${params.output}/${params.minimap2_dir}", pattern: "*.bai" }
    else { publishDir "${params.output}/${params.minimap2_dir}", mode: 'copy', pattern: "*.bai" }

    input:
    tuple val(sample_name), path(bam_file)

    output:
    path("${bam_file}.bai")

    script:
    """
    samtools index ${bam_file}
    """
}
