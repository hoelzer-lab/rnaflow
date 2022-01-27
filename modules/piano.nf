process piano {
    label 'deseq2'

    errorStrategy 'retry'
    maxRetries 1

    input:
    path(piano_script)
    path(resFold05)
    val(species)
    val(id_type)
    path(script_improve_deseq_table)

    output:
    path("downstream_analysis/piano")
    
    script:
    
    """
    R CMD BATCH --no-save --no-restore '--args c(".") c(${resFold05}) c(${species}) c(${id_type}) c(${task.cpus}) ${piano_script}
    """
}
