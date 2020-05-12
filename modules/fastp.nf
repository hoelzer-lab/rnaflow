/************************************************************************
* TRIMMING
*
* TODO: pimp the trimming command for adapters, sliding-window QC, ...
************************************************************************/
process fastp {
    label 'fastp'
    
    if (params.cloudProcess) { publishDir "${params.output}/${params.fastp_dir}", mode: 'copy', pattern: "*.trimmed.fastq" }
    else { publishDir "${params.output}/${params.fastp_dir}", pattern: "*.trimmed.fastq" }

    input:
    tuple val(name), path(reads)

    output:
    tuple val(name), path("${name}*.trimmed.fastq"), emit: sample_trimmed
    path "${name}_fastp.json", emit: json_report

    script:
    if (params.mode == 'single') {
    """
    fastp -i ${reads[0]} -o ${name}.trimmed.fastq -n 5 --thread ${task.cpus} --json ${name}_fastp.json
    """
    }
    else {
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${name}.R1.trimmed.fastq -O ${name}.R2.trimmed.fastq -n 5 --thread ${task.cpus} --json ${name}_fastp.json
    """
    }
}
