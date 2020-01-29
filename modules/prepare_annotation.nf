/**************************************************
* PREPARE ANNOTATION FOR LATER INPUT AND USAGE
***************************************************/
process prepare_annotation {
    label 'python3'
    publishDir "${params.output}/${params.dir}", mode: 'copy', pattern: "${annotation.baseName}.id2ensembl"

    input: 
    file(annotation)

    output:
    file("${annotation.baseName}.id2ensembl")

    shell:
    '''
#!/usr/bin/env python3
with open("!{annotation}", 'r') as gtf, open("!{annotation.baseName}.id2ensembl", 'a') as out:
  for line in gtf:
    if not line.startswith('#'):
      split_line = line.split('\\t')
      if split_line[2] == 'gene':
        desc = split_line[8]
        gene_id = line.split('gene_id')[1].split(';')[0].replace('"', '').strip()
        if 'gene_name' in line:
          gene_name = desc.split('gene_name')[1].split(';')[0].replace('"','').strip()
        else:
          gene_name = gene_id
        gene_biotype = desc.split('gene_biotype')[1].split(';')[0].replace('"','').strip()
        out.write('\\t'.join([gene_id, gene_name, gene_biotype]) + '\\n')
'''
}
