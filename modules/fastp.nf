/************************************************************************
* TRIMMING
*
* TODO: pimp the trimming command for adapters, sliding-window QC, ...
************************************************************************/
process fastp {
    label 'fastp'
    
    publishDir "${params.output}/${params.fastp_dir}", mode: 'copy', pattern: "*.trimmed.fastq.gz"

    input:
    tuple val(name), path(reads)

    output:
    tuple val(name), path("${name}*.trimmed.fastq.gz"), emit: sample_trimmed
    path "${name}_fastp.json", emit: json_report

    script:
    if (params.mode == 'single') {
    """
    fastp -i ${reads[0]} -o ${name}.trimmed.fastq.gz --thread ${task.cpus} --json ${name}_fastp.json ${params.fastp_additional_params}
    """
    }
    else {
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${name}.R1.trimmed.fastq.gz -O ${name}.R2.trimmed.fastq.gz --thread ${task.cpus} --json ${name}_fastp.json ${params.fastp_additional_params}
    """
    }
}
