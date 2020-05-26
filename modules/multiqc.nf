/************************************************************************
* MultiQC
************************************************************************/
process multiqc {
    label 'multiqc'
    label 'smallTask'

    if (params.cloudProcess) { publishDir "${params.output}/${params.multiqc_dir}", mode: 'copy' }
    else { publishDir "${params.output}/${params.multiqc_dir}" }

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
