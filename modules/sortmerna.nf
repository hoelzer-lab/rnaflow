/************************************************************************
* SORTMERNA
*
* Remove rRNA reads
************************************************************************/
process sortmerna {
  conda 'envs/sortmerna.yaml'
  publishDir "${params.output}/${params.dir}", mode: 'copy', pattern: "${name}.other.fastq"
  publishDir "${params.output}/${params.dir}", mode: 'copy', pattern: "${name}.aligned.log"

  input:
  tuple val(name), file(reads)
  file(db)

  output:
  tuple val(name), file("${name}*.other.fastq"), emit: no_rna_fastq
  path "${name}.aligned.log", emit: log

  script:
  if (params.mode == 'single') {
  """
  sortmerna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./rRNA_databases/silva-bac-16s-id90:./rRNA_databases/silva-bac-23s-id98.fasta,./rRNA_databases/silva-bac-23s-id98:./rRNA_databases/silva-arc-16s-id95.fasta,./rRNA_databases/silva-arc-16s-id95:./rRNA_databases/silva-arc-23s-id98.fasta,./rRNA_databases/silva-arc-23s-id98:./rRNA_databases/silva-euk-18s-id95.fasta,./rRNA_databases/silva-euk-18s-id95:./rRNA_databases/silva-euk-28s-id98.fasta,./rRNA_databases/silva-euk-28s-id98:./rRNA_databases/rfam-5s-database-id98.fasta,./rRNA_databases/rfam-5s-database-id98:./rRNA_databases/rfam-5.8s-database-id98.fasta,./rRNA_databases/rfam-5.8s-database-id98 \
--reads ${reads[0]} \
--aligned ${name}.aligned \
--other ${name}.other \
--sam --fastx --log --blast 1 --num_alignments 1 -v \
-a ${params.cores}
  """
  }
  else {
  """
  merge-paired-reads.sh ${reads[0]} ${reads[1]} ${name}.merged.fastq
  sortmerna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./rRNA_databases/silva-bac-16s-id90:./rRNA_databases/silva-bac-23s-id98.fasta,./rRNA_databases/silva-bac-23s-id98:./rRNA_databases/silva-arc-16s-id95.fasta,./rRNA_databases/silva-arc-16s-id95:./rRNA_databases/silva-arc-23s-id98.fasta,./rRNA_databases/silva-arc-23s-id98:./rRNA_databases/silva-euk-18s-id95.fasta,./rRNA_databases/silva-euk-18s-id95:./rRNA_databases/silva-euk-28s-id98.fasta,./rRNA_databases/silva-euk-28s-id98:./rRNA_databases/rfam-5s-database-id98.fasta,./rRNA_databases/rfam-5s-database-id98:./rRNA_databases/rfam-5.8s-database-id98.fasta,./rRNA_databases/rfam-5.8s-database-id98 \
--reads ${name}.merged.fastq \
--paired_in --aligned ${name}.aligned \
--other ${name}.other_merged \
--sam --fastx --log --blast 1 --num_alignments 1 -v \
-a ${params.cores}
  unmerge-paired-reads.sh ${name}.other_merged.fastq ${name}.R1.other.fastq ${name}.R2.other.fastq
  """
  }
}