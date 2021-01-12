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
    import sys
    target_id_set = set()
    with open("!{annotation}", 'r') as gtf, open("!{annotation.baseName}.id2details", 'a') as out:
        for line in gtf:
            if not line.startswith('#'):
                split_line = line.split('\\t')
                if split_line[2] == '!{gtf_feature_type_of_attr_type}' or split_line[2] == 'pseudogene':
                    if '!{gtf_attr_type}' not in line:
                        sys.exit(f"ERROR: No '!{gtf_attr_type}' found in a '!{gtf_feature_type_of_attr_type}'-type line:\\n{line}\\nCheck your annotation file or adapt the parameter `--featurecounts_additional_params`")
                    target_id = line.split('!{gtf_attr_type}')[1].split(';')[0].replace('"', '').strip()
                    if target_id not in target_id_set:
                        desc = split_line[8]
                        chr = split_line[0]
                        start = split_line[3]
                        stop = split_line[4]
                        strand = split_line[6]
                        if '!{gtf_attr_type}' == 'transcript_id' and 'transcript_name' in desc:
                            target_name = desc.split('transcript_name')[1].split(';')[0].replace('"','').strip()
                        else:
                            if 'gene_name' in desc:
                                target_name = desc.split('gene_name')[1].split(';')[0].replace('"','').strip()
                            else:
                                target_name = target_id
                        if target_name == 'NA':
                            target_name = target_id
                        if '!{gtf_attr_type}' == 'transcript_id' and 'transcript_biotype' in desc:
                            target_biotype = desc.split('transcript_biotype')[1].split(';')[0].replace('"','').strip()
                        else:
                            if 'gene_biotype' in desc:
                                target_biotype = desc.split('gene_biotype')[1].split(';')[0].replace('"','').strip()
                            else:
                                target_biotype = 'NA'
                        target_id_set.add(target_id)
                        out.write('\\t'.join([target_id, target_name, target_biotype, chr, start, stop, strand, desc.rstrip()]) + '\\n')
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
