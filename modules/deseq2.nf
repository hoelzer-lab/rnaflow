/**************************************************
* DESEQ2
*
* ToDo: 
* 1) here main parts of the R script included in scripts/deseq2.R should be executed
* 2) usage of a R conda environment with the libraries should work, I did not test it fully
***************************************************/
process deseq2 {
    conda 'envs/deseq2.yaml'
    publishDir "${params.output}/${params.dir}", mode: 'copy', pattern: "*"

    input:
    tuple val(name), file(counts)
    file(ensembl2id)

    //output:
    //tuple val(name), file(counts)

    script:
    """
    echo ${name}.counts.formated
    """
}

/*
* this is how the input for the R script in scripts/ should be formated
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
