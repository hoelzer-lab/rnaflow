/**************************************************
* DESEQ2
***************************************************/
process deseq2 {
    label 'deseq2'
    time '3h'

    //errorStrategy 'retry'
    //maxRetries 1

    if ( params.softlink_results ) { publishDir "${params.output}/${params.deseq2_dir}", pattern: "*" }
    else { publishDir "${params.output}/${params.deseq2_dir}", mode: 'copy', pattern: "*" }

    input:
    path(regionReport_config)
    path(fc_counts_formated)
    val(condition)
    val(col_labels)
    val(comparisons)
    path(ensembl2id)
    path(annotation_genes)
    val(sources)
    val(species)
    path(script)
    val(id_type)
    path(script_refactor_reportingtools_table)
    path(script_improve_deseq_table)

    output:
        path("*")
        path "*_vs_*/results/deseq2_*_filtered_padj_0.05.csv", glob: true , emit: resFold05

    
    script:

    sample_files = fc_counts_formated.collect { "\"${it}\"" }.join(",")
    col_labels_in = col_labels.collect { "\"${it}\"" }.join(",")
    conditions = condition.collect { "\"${it}\"" }.join(",")
    levels = condition.collect { "\"${it}\"" }.toSet().join(",")
    if ( sources.toSet().size() == 1 && ! sources.toSet()[0] ) {
        // sources is a list with empty emlements ([, , , ])
        sources_in = ''
    } else {
        sources_in = sources.collect { "\"${it}\"" }.join(",")
    }
    """
    R CMD BATCH --no-save --no-restore '--args c(".") c(${sample_files}) c(${conditions}) c(${col_labels_in}) c(${levels}) c(${comparisons}) c("${ensembl2id}") c("${annotation_genes}") c(${sources_in}) c("${species}") c("${regionReport_config}") c(${task.cpus}) c("${id_type}")' ${script}
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
sources <- c()
*/
