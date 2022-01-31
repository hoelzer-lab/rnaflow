process piano {
    label 'deseq2'
    tag "$comparison"

    publishDir "${params.output}/${params.deseq2_dir}/${comparison}/downstream_analysis", mode: 'copy', pattern: "piano"

    errorStrategy 'retry'
    maxRetries 1 

    input:
    path(piano_script)
    each(resFold05)
    val(species)
    val(id_type)
    path(script_improve_deseq_table)

    output:
    path("piano")
    
    script:
    comparison = resFold05.toString().findAll("[a-zA-Z]*_vs_[a-zA-Z]*")[0]

    """
    R CMD BATCH --no-save --no-restore '--args c(".") c("${resFold05}") c("${species}") c("${id_type}") c("${task.cpus}")' ${piano_script}
    """
}
