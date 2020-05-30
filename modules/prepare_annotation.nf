/**************************************************
* PREPARE ANNOTATION FOR LATER INPUT AND USAGE
***************************************************/
process format_annotation {
    label 'python3'
    label 'smallTask'

    if (params.cloudProcess) { publishDir "${params.output}/${params.annotation_dir}", mode: 'copy', pattern: "*.id2ensembl" }
    else { publishDir "${params.output}/${params.annotation_dir}", pattern: "*.id2ensembl" }

    input: 
    path(annotation)

    output:
    path("${annotation.baseName}.id2ensembl")

    shell:
    '''
    #!/usr/bin/env python3
    with open("!{annotation}", 'r') as gtf, open("!{annotation.baseName}.id2ensembl", 'a') as out:
        for line in gtf:
            if not line.startswith('#'):
                split_line = line.split('\\t')
                if split_line[2] == 'gene' or split_line[2] == 'pseudogene':
                    desc = split_line[8]
                    gene_id = line.split('gene_id')[1].split(';')[0].replace('"', '').strip()
                    if 'gene_name' in desc:
                        gene_name = desc.split('gene_name')[1].split(';')[0].replace('"','').strip()
                    else:
                        gene_name = gene_id
                    if 'gene_biotype' in desc:
                        gene_biotype = desc.split('gene_biotype')[1].split(';')[0].replace('"','').strip()
                    else:
                        gene_biotype = 'NA'
                    out.write('\\t'.join([gene_id, gene_name, gene_biotype]) + '\\n')
    '''
}

/**************************************************
* PREPARE ANNOTATION FOR LATER INPUT AND USAGE
***************************************************/
process format_annotation_gene_rows {
    label 'python3'
    label 'smallTask'
    publishDir "${params.output}/${params.annotation_dir}", mode: 'copy', pattern: "${annotation.baseName}.gene.gtf"

    input: 
    path(annotation)

    output:
    path("${annotation.baseName}.gene.gtf")

    script:
    """
    awk '{if(\$3=="gene" || \$3=="pseudogene"){print \$0}}' ${annotation} > ${annotation.baseName}.gene.gtf
    """
}
