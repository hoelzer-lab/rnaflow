/************************************************************************
* TRIMMING
*
* TODO: pimp the trimming command for adapters, sliding-window QC, ...
************************************************************************/
process fastp {
  conda 'envs/fastp.yaml'
  publishDir "${params.output}/${params.dir}", mode: 'copy', pattern: "${name}*.trimmed.fastq"

  input:
  tuple val(name), file(reads)

  output:
  tuple val(name), file("${name}*.trimmed.fastq")

  script:
  if (params.mode == 'single') {
  """
  fastp -i ${reads[0]} -o ${name}.trimmed.fastq -n 5 --thread ${params.cores}  
  """
  }
  else {
  """
  fastp -i ${reads[0]} -I ${reads[1]} -o ${name}.R1.trimmed.fastq -O ${name}.R2.trimmed.fastq -n 5 --thread ${params.cores}  
  """         
  }
}
