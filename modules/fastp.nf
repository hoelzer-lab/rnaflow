/************************************************************************
* TRIMMING
*
* TODO: pimp the trimming command for adapters, sliding-window QC, ...
************************************************************************/
process fastp {
    label 'fastp'
    
    if (params.cloudProcess) { publishDir "${params.output}/${params.fastp_dir}", mode: 'copy', pattern: "*.trimmed.fastq.gz" }
    else { publishDir "${params.output}/${params.fastp_dir}", pattern: "*.trimmed.fastq.gz" }

    input:
    tuple val(name), path(reads)
    val(additionalParams)

    output:
    tuple val(name), path("${name}*.trimmed.fastq.gz"), emit: sample_trimmed
    path "${name}_fastp.json", emit: json_report

    script:
    if (params.mode == 'single') {
    """
    fastp -i ${reads[0]} -o ${name}.trimmed.fastq.gz -n 5 --thread ${task.cpus} --json ${name}_fastp.json -z 6 ${additionalParams}
    """
    }
    else {
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${name}.R1.trimmed.fastq.gz -O ${name}.R2.trimmed.fastq.gz -n 5 --thread ${task.cpus} --json ${name}_fastp.json -z 6 ${additionalParams}
    """
    }
}
