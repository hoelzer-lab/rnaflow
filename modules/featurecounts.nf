/************************************************************************
* COUNTING WITH FEATURECOUNTS
* TODO: check that we do not miss a gene due to the 'sed 1d' remove of the first two lines
************************************************************************/
process featurecounts {
  conda 'envs/subread.yaml'
  publishDir "${params.output}/${params.dir}", mode: 'copy', pattern: "${name}.counts*"

  input:
  tuple val(name), file(bam)
  file(annotation)

  output:
  tuple val(name), file("${name}.counts.formated"), emit: formated_counts // [mock_rep1, /home/hoelzer/git/nanozoo/wf_gene_expression/work/9e/7fb58903c9e4163d526ef749c0d088/mock_rep1.counts.formated]
  path "${name}.counts.summary", emit: log

  shell:
  if (params.mode == 'single') {
  '''
  featureCounts -T !{params.cores} -s !{params.strand} -a !{annotation} -o !{name}.counts !{bam}
  awk '{print $1"\t"$7}' !{name}.counts | sed 1d | sed 1d > !{name}.counts.formated
  '''
  }
  else {
  '''
  featureCounts -pBP -T !{params.cores} -s !{params.strand} -a !{annotation} -o !{name}.counts !{bam}
  awk '{print $1"\t"$7}' !{name}.counts | sed 1d | sed 1d > !{name}.counts.formated
  '''
  }
}
