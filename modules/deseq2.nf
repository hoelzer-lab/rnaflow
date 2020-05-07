/**************************************************
* DESEQ2
***************************************************/
process deseq2 {
    label 'deseq2'
    publishDir "${params.output}/${params.deseq2_dir}", mode: 'copy', pattern: "*"

    input:
    path(fc_counts_formated)
    val(condition)
    val(col_labels)
    val(comparisons)
    path(ensembl2id)
    path(annotation_genes)
    val(patients)
    path(script)
    path(script_refactor_reportingtools_table)
    path(script_improve_deseq_table)

    output:
    path("*")
    
    script:

    sample_files = fc_counts_formated.collect { "\"${it}\"" }.join(",")
    col_labels = col_labels.collect { "\"${it}\"" }.join(",")
    conditions = condition.collect { "\"${it}\"" }.join(",")
    levels = condition.collect { "\"${it}\"" }.toSet().join(",")
    if ( patients.toSet().size() == 1 && ! patients.toSet()[0] ) {
        // patients is a list with empty emlements ([, , , ])
        patients = ''
    } else {
        patients = patients.collect { "\"${it}\"" }.join(",")
    }

    """
    R CMD BATCH --no-save --no-restore '--args c(".") c(${sample_files}) c(${conditions}) c(${col_labels}) c(${levels}) c(${comparisons}) c("${ensembl2id}") c("${annotation_genes}") c(${patients})' ${script}
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
