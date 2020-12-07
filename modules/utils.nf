process extract_tar_bz2 {
    input:
    path(tar_bz2)

    output:
    path("*")

    script:
    """
    tar zxvf ${tar_bz2}
    """
}