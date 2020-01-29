/**************************************************
* PREPARE ANNOTATION FOR LATER INPUT AND USAGE
***************************************************/
process prepare_annotation_gene_rows {
    label 'python3'
    publishDir "${params.output}/${params.dir}", mode: 'copy', pattern: "${annotation.baseName}.gene.gtf"

    input: 
    file(annotation)

    output:
    file("${annotation.baseName}.gene.gtf")

    script:
    """
    awk '{if(\$3=="gene"){print \$0}}' ${annotation} > ${annotation.baseName}.gene.gtf
    """
}
