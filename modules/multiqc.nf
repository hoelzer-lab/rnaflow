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
    path(fastqc)

    output:
    path "*multiqc_report.html"
    path "multiqc_data"

    script:
    """
    for i in *.tar.gz; do tar zxvf \$i; done 
    cp fastqc_*/* .
    multiqc .    
    """
}
