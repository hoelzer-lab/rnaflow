/**************************************************
* PREPARE ANNOTATION FOR LATER INPUT AND USAGE
***************************************************/
process format_annotation {
    label 'python3'
    label 'smallTask'

    if ( params.softlink_results ) { publishDir "${params.output}/${params.annotation_dir}", pattern: "*.id2details" }
    else { publishDir "${params.output}/${params.annotation_dir}", mode: 'copy', pattern: "*.id2details" }

    input: 
    path(annotation)
    val(gtf_attr_type)
    val(gtf_feature_type_of_attr_type)

    output:
    path("${annotation.baseName}.id2details")

    shell:
    '''
    #!/usr/bin/env python3
    with open("!{annotation}", 'r') as gtf, open("!{annotation.baseName}.id2details", 'a') as out:
        for line in gtf:
            if not line.startswith('#'):
                split_line = line.split('\\t')
                if split_line[2] == '!{gtf_feature_type_of_attr_type}' or split_line[2] == 'pseudogene':
                    desc = split_line[8]
                    target_id = line.split('!{gtf_attr_type}')[1].split(';')[0].replace('"', '').strip()
                    if 'gene_name' in desc:
                        gene_name = desc.split('gene_name')[1].split(';')[0].replace('"','').strip()
                    else:
                        gene_name = target_id
                    if 'gene_biotype' in desc:
                        gene_biotype = desc.split('gene_biotype')[1].split(';')[0].replace('"','').strip()
                    else:
                        gene_biotype = 'NA'
                    out.write('\\t'.join([target_id, gene_name, gene_biotype]) + '\\n')
    '''
}

/**************************************************
* PREPARE ANNOTATION FOR LATER INPUT AND USAGE
***************************************************/
process format_annotation_gene_rows {
    label 'smallTask'

    if ( params.softlink_results ) { publishDir "${params.output}/${params.annotation_dir}", pattern: "*.gtf" }
    else { publishDir "${params.output}/${params.annotation_dir}", mode: 'copy', pattern: "*.gtf" }

    input: 
    path(annotation)
    val(gtf_feature_type_of_attr_type)

    output:
    path("${annotation.baseName}.${gtf_feature_type_of_attr_type}.gtf")

    script:
    """
    awk '{if(\$3=="${gtf_feature_type_of_attr_type}" || \$3=="pseudogene"){print \$0}}' ${annotation} > ${annotation.baseName}.${gtf_feature_type_of_attr_type}.gtf
    """
}
