/************************************************************************
* MultiQC
************************************************************************/
process multiqc {
    label 'multiqc'
    label 'smallTask'

    if (params.cloudProcess) { publishDir "${params.output}/${params.multiqc_dir}", mode: 'copy' }
    else { publishDir "${params.output}/${params.multiqc_dir}" }

    input:
    path(config)
    path(sample_names)
    path(fastp)
    path(sortmerna)
    path(hisat2)
    path(featurecounts)
    path(fastqcPre)
    path(fastqcPost)

    output:
    path "*multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc . -c ${config} --sample-names ${sample_names}
    """
}

process multiqc_sample_names {
    label 'smallTask'

    input:
    val(list)

    output:
    path('multiqc_sample_names.tsv')

    script:
    def tbl = ''
    if (params.mode == 'single') {
        for( int i=0; i<list.size()-1; i=i+2 ){
            tbl += "${list[i+1].baseName.tokenize('.')[0]}\t${list[i]}\n"
        }
    }
    else {
        for( int i=0; i<list.size()-2; i=i+3 ){
            tbl += "${list[i+1].baseName.tokenize('.')[0]}\t${list[i]}.R1\n"
            tbl += "${list[i+2].baseName.tokenize('.')[0]}\t${list[i]}.R2\n"
        }
    }
    """
    echo "Read File Names\tSample Names" > multiqc_sample_names.tsv
    echo "${tbl}" >> multiqc_sample_names.tsv
    """
}