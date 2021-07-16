/************************************************************************
* COUNTING WITH FEATURECOUNTS
* TODO: check that we do not miss a gene due to the 'sed 1d' remove of the first two lines
************************************************************************/
process featurecounts {
    label 'subread'

    if ( params.softlink_results ) { publishDir "${params.output}/${params.featurecounts_dir}", pattern: "*.tsv" }
    else { publishDir "${params.output}/${params.featurecounts_dir}", mode: 'copy', pattern: "*.tsv" }

    input:
    tuple val(name), path(bam)
    path(annotation)
    val(additionalParams)

    output:
    tuple val(name), path("${name}.counts.tsv"), emit: counts // [mock_rep1, /home/hoelzer/git/nanozoo/wf_gene_expression/work/9e/7fb58903c9e4163d526ef749c0d088/mock_rep1.tsv]
    path "${name}.counts.tsv.summary", emit: log

    script:
    if (params.mode == 'single') {
        if (params.nanopore == true) {
        """
        featureCounts -L -T ${task.cpus} -s ${params.strand} -a ${annotation} -o ${name}.counts.tsv ${additionalParams} ${bam}
        """
        } else {
            """
            featureCounts -T ${task.cpus} -s ${params.strand} -a ${annotation} -o ${name}.counts.tsv  ${additionalParams} ${bam}
            """
        }
    }
    else {
    """
    featureCounts -pBP -T ${task.cpus} -s ${params.strand} -a ${annotation} -o ${name}.counts.tsv ${additionalParams} ${bam}
    """
    }
}
