/************************************************************************
* TRIMMING
*
* TODO: pimp the trimming command for adapters, sliding-window QC, ...
************************************************************************/
process fastp {
    label 'fastp'
    tag "$meta.sample"

    if ( params.softlink_results ) { publishDir "${params.output}/${params.fastp_dir}", pattern: "*.trimmed.fastq.gz" }
    else { publishDir "${params.output}/${params.fastp_dir}", mode: 'copy', pattern: "*.trimmed.fastq.gz" }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.sample}*.trimmed.fastq.gz"), emit: sample_trimmed
    path "${meta.sample}_fastp.json", emit: json_report

    script:
    if ( !meta.paired_end ) {
    """
    fastp -i ${reads[0]} -o ${meta.sample}.trimmed.fastq.gz --thread ${task.cpus} --json ${meta.sample}_fastp.json ${params.fastp_additional_params}
    """
    }
    else {
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${meta.sample}.R1.trimmed.fastq.gz -O ${meta.sample}.R2.trimmed.fastq.gz --thread ${task.cpus} --json ${meta.sample}_fastp.json ${params.fastp_additional_params}
    """
    }
}
