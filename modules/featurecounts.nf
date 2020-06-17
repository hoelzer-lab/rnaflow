/************************************************************************
* COUNTING WITH FEATURECOUNTS
* TODO: check that we do not miss a gene due to the 'sed 1d' remove of the first two lines
************************************************************************/
process featurecounts {
    label 'subread'

    if (params.cloudProcess) { publishDir "${params.output}/${params.featurecounts_dir}", mode: 'copy', pattern: "*.tsv" }
    else { publishDir "${params.output}/${params.featurecounts_dir}", pattern: "*.tsv" }

    input:
    tuple val(name), path(bam)
    path(annotation)

    output:
    tuple val(name), path("${name}.counts.tsv"), emit: counts // [mock_rep1, /home/hoelzer/git/nanozoo/wf_gene_expression/work/9e/7fb58903c9e4163d526ef749c0d088/mock_rep1.tsv]
    path "${name}.counts.tsv.summary", emit: log

    script:
    if (params.mode == 'single') {
    """
    featureCounts -T ${task.cpus} -s ${params.strand} -a ${annotation} -o ${name}.counts.tsv ${bam}
    """
    }
    else {
    """
    featureCounts -pBP -T ${task.cpus} -s ${params.strand} -a ${annotation} -o ${name}.counts.tsv ${bam}
    """
    }
}
