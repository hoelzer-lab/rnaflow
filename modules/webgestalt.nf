process webgestalt {
    label 'webgestalt'
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
    path("WebGestalt") optional true
    
    script:
    comparison = resFold05.toString().split("deseq2_")[1].split("_filtered")[0]//findAll("[a-zA-Z]*_vs_[a-zA-Z]*")[0]
    """
    # set env var for zip command for zipping in webgestalt since sometimes its not set in conda envs somehow and crashes the report making
    export R_ZIPCMD=\$(which zip)
    R CMD BATCH --no-save --no-restore '--args c(".") c("${resFold05}") c("${species}") c("${id_type}")' ${webgestalt_script}
    """
}
