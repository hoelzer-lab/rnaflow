/**************************************************
* DESEQ2
***************************************************/
process deseq2 {
    conda 'envs/deseq2.yaml'
    publishDir "${params.output}/${params.dir}", mode: 'copy', pattern: "*"

    input:
    val(name)
    file(counts)
    file(ensembl2id)
    file(annotation_genes)
    file(script)
    file(script_refactor_reportingtools_table)
    file(script_improve_deseq_table)

    output:
    path("*")
    
    script:

    count_files = counts.collect { "\"${it}\"" }.join(",")
    col_labels = name.collect { "\"${it}\"" }.join(",")
    conditions = name.collect { "${it}".tokenize('_')[0] }.collect { "\"${it}\"" }.join(",")
    levels = name.collect { "${it}".tokenize('_')[0] }.collect { "\"${it}\"" }.toSet().join(",")
    comparisons = "\"" + name.collect { "${it}".tokenize('_')[0] }.toSet().join(":") + "\""

    """
    R CMD BATCH --no-save --no-restore '--args c(".") c(${count_files}) c(${conditions}) c(${col_labels}) c(${levels}) c(${comparisons}) c("${ensembl2id}") c("${annotation_genes}") c("${params.species}") c()' ${script}
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
