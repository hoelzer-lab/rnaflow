/**************************************************
* DESEQ2
*
* ToDo: 
***************************************************/
process deseq2 {
  publishDir "${params.output}/${params.dir}", mode: 'copy', pattern: "*"

  input:
  tuple val(name), file(counts)
  file(annotation)
  file(ensembl2id)

  //output:
  //tuple val(name), file(counts) into final_ch 

  script:
  """
  echo ${name}.counts.formated
  """
}

/*
project_dir <- "results/" 
samples <- c("sample1","sample2","sample3","sample4","sample5","sample6") 
conditions <- c("mock","mock","mock","treated","treated","treated") 
col.labels <- c("mock_rep1","mock_rep2","mock_rep3","treated_rep1","treated_rep2","treated_rep3") 
levels <- c("mock","treated") 
comparisons <- c("mock:treated")
ensembl2genes <- "${ensembl2id}"
species <- "${params.}"
patients <- c()
*/

}