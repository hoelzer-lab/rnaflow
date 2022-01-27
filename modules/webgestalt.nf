process webgestalt {
    label 'deseq2'

    errorStrategy 'retry'
    maxRetries 1

    input:
    path(webgestalt_script)
    path(resFold05)
    val(species)
    val(id_type)
    path(script_improve_deseq_table)
    val(l1)
    val(l2)

    output:
    path("*vs*/downstream_analysis/piano", glob: True)
    
    script:
    """
    R CMD BATCH --no-save --no-restore '--args c(".") c(${resFold05}) c(${species}) c(${id_type}) c(${l1}) c(${l2})' ${webgestalt_script}
    """
}
