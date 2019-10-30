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
    val(name)
    file(counts)
    file(ensembl2id)

    output:
    stdout()
    //tuple val(name), file(counts)
    
    // exec:
    // println "$name"
    // println "$counts"
    // println "${params.species}"

    script:
    """
    R CMD BATCH '--args c(${params.output}/${params.dir}/") c(${counts}) c("treat","mock","treat","mock","mock","treat") c(${name}) c("treat","mock") c("mock:treat") c("${ensembl2id}") c("${params.species}") c()' scripts/deseq2.R 
    """
}
// R CMD BATCH '--args c("test/") c("results/03-Counting/featurecounts/treat_rep1.counts.formated","results/03-Counting/featurecounts/mock_rep1.counts.formated","results/03-Counting/featurecounts/treat_rep3.counts.formated","results/03-Counting/featurecounts/mock_rep2.counts.formated","results/03-Counting/featurecounts/mock_rep3.counts.formated","results/03-Counting/featurecounts/treat_rep2.counts.formated") c("treat","mock","treat","mock","mock","treat") c("treat_rep1","mock_rep1","treat_rep3","mock_rep2","mock_rep3","treat_rep2") c("treat","mock") c("mock:treat") c("results/04-Annotation/eco.id2ensembl") c("eco") c()' scripts/deseq2.R 
/*
project_dir <- "${params.output}/${params.dir}/" 
samples <- "$counts"
conditions <-
col.labels <- "$name"
levels <-
comparisons <-
ensembl2genes <- "${ensembl2id}"
species <- "${params.species}"
patients <- c()
*/

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
