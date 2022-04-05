process extract_tar_bz2 {
    label 'basic_tools'
    if (!params.cloudProcess) { label 'smallTask' }

    input:
    path(tar_bz2)

    output:
    path("*")

    script:
    """
    tar zxvf ${tar_bz2}
    """
}