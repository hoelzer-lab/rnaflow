/************************************************************************
* MultiQC
************************************************************************/
process multiqc {
    conda 'envs/multiqc.yaml'
    publishDir "${params.output}/${params.dir}", mode: 'copy'

    input:
    path(fastp)
    path(sortmerna)
    path(hisat2)
    path(featurecounts)

    output:
    path "*multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc .    
    """
}
