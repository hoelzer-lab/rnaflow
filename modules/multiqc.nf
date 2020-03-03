/************************************************************************
* MultiQC
************************************************************************/
process multiqc {
    label 'multiqc'
    publishDir "${params.output}/${params.multiqc_dir}", mode: 'copy'

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
