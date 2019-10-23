/************************************************************************
* HISAT2 INDEX
************************************************************************/
process hisat2index {
  conda 'envs/hisat2.yaml'
  if (params.cloudProcess) { publishDir "${params.cloudDatabase}/genomes/${params.reference}", mode: 'copy', pattern: "${params.reference}*.ht2" }
  else { storeDir "nextflow-autodownload-databases/genomes/${params.reference}" }  

  input:
  file(reference)

  output:
  tuple file(reference), file("${params.reference}*.ht2")

  script:
  """
  hisat2-build -p ${params.cores} ${reference} ${reference.baseName}
  """
}