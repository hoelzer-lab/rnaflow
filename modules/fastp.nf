/************************************************************************
* TRIMMING
*
* TODO: pimp the trimming command for adapters, sliding-window QC, ...
************************************************************************/
process fastp {
    label 'fastp'

    publishDir "${params.output}/${params.fastp_dir}", mode: 'copy', pattern: "${name}*.trimmed.fastq"

    input:
    tuple val(name), file(reads)

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
