process webgestalt {
    label 'deseq2'
    tag "$comparison"

    publishDir "${params.output}/${params.deseq2_dir}/${comparison}/downstream_analysis", mode: 'copy', pattern: "WebGestalt"

    //errorStrategy 'retry'
    //maxRetries 1 

    input:
    path(webgestalt_script)
    each(path(resFold05))
    val(species)
    val(id_type)

    output:
    path("WebGestalt")
    
    script:
    comparison = resFold05.toString().findAll("[a-zA-Z]*_vs_[a-zA-Z]*")[0]
    """
    R CMD BATCH --no-save --no-restore '--args c(".") c("${resFold05}") c("${species}") c("${id_type}")' ${webgestalt_script}
    """
}
