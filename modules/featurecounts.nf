/************************************************************************
* COUNTING WITH FEATURECOUNTS
* TODO: check that we do not miss a gene due to the 'sed 1d' remove of the first two lines
************************************************************************/
process featurecounts {
    label 'subread'
    tag "$meta.sample"

    if ( params.softlink_results ) { publishDir "${params.output}/${params.featurecounts_dir}", pattern: "*.tsv" }
    else { publishDir "${params.output}/${params.featurecounts_dir}", mode: 'copy', pattern: "*.tsv" }

    input:
    tuple val(meta), path(bam)
    path(annotation)
    val(additionalParams)

    output:
    tuple val(meta.sample), path("${meta.sample}.counts.tsv"), emit: counts // [mock_rep1, /home/hoelzer/git/nanozoo/wf_gene_expression/work/9e/7fb58903c9e4163d526ef749c0d088/mock_rep1.tsv]
    path "${meta.sample}.counts.tsv.summary", emit: log

    script:
    if ( !meta.paired ) {
        if (params.nanopore == true) {
        """
        featureCounts -L -T ${task.cpus} -s ${meta.strandedness} -a ${annotation} -o ${meta.sample}.counts.tsv ${additionalParams} ${bam}
        """
        } else {
            """
            featureCounts -T ${task.cpus} -s ${meta.strandedness} -a ${annotation} -o ${meta.sample}.counts.tsv  ${additionalParams} ${bam}
            """
        }
    }
    else {
    """
    featureCounts -pBP -T ${task.cpus} -s ${meta.strandedness} -a ${annotation} -o ${meta.sample}.counts.tsv ${additionalParams} ${bam}
    """
    }
}
