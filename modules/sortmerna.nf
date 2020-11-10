/************************************************************************
* SORTMERNA
*
* Remove rRNA reads
************************************************************************/
process sortmerna {
    label 'sortmerna'
    
    publishDir "${params.output}/${params.sortmerna_dir}", mode: 'copy', pattern: "*.other.fastq.gz"

    input:
    tuple val(name), path(reads)
    path(db)

    output:
    tuple val(name), path("${name}*.other.fastq.gz"), emit: no_rna_fastq
    path "${name}.aligned.log", emit: log

    script:
    if (params.mode == 'single') {
    assert "$reads".split('.gz').size() == 1 : "suffix should be .gz"
    uncompr_reads = "$reads".split('.gz')[0]
    """
    unpigz -f -p ${task.cpus} ${reads[0]}
    sortmerna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./rRNA_databases/silva-bac-16s-id90:./rRNA_databases/silva-bac-23s-id98.fasta,./rRNA_databases/silva-bac-23s-id98:./rRNA_databases/silva-arc-16s-id95.fasta,./rRNA_databases/silva-arc-16s-id95:./rRNA_databases/silva-arc-23s-id98.fasta,./rRNA_databases/silva-arc-23s-id98:./rRNA_databases/silva-euk-18s-id95.fasta,./rRNA_databases/silva-euk-18s-id95:./rRNA_databases/silva-euk-28s-id98.fasta,./rRNA_databases/silva-euk-28s-id98:./rRNA_databases/rfam-5s-database-id98.fasta,./rRNA_databases/rfam-5s-database-id98:./rRNA_databases/rfam-5.8s-database-id98.fasta,./rRNA_databases/rfam-5.8s-database-id98 \
    --reads ${uncompr_reads} \
    --aligned ${name}.aligned \
    --other ${name}.other \
    --fastx --log --num_alignments 1 -v \
    -a ${task.cpus}
    pigz -p ${task.cpus} ${name}.other.fastq
    rm ${name}.aligned.fastq ${uncompr_reads}
    """
    }
    else {
    split_R1 = "$reads".split(' ')[0].split('.gz')
    split_R2 = "$reads".split(' ')[1].split('.gz')
    assert split_R1.size() == 1 : "suffix should be .gz" 
    assert split_R2.size() == 1 : "suffix should be .gz"
    uncompr_reads_R1 = split_R1[0]
    uncompr_reads_R2 = split_R2[0]
    """
    unpigz -f -p ${task.cpus} ${reads[0]}
    unpigz -f -p ${task.cpus} ${reads[1]}
    merge-paired-reads.sh ${uncompr_reads_R1} ${uncompr_reads_R2} ${name}.merged.fastq
    sortmerna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./rRNA_databases/silva-bac-16s-id90:./rRNA_databases/silva-bac-23s-id98.fasta,./rRNA_databases/silva-bac-23s-id98:./rRNA_databases/silva-arc-16s-id95.fasta,./rRNA_databases/silva-arc-16s-id95:./rRNA_databases/silva-arc-23s-id98.fasta,./rRNA_databases/silva-arc-23s-id98:./rRNA_databases/silva-euk-18s-id95.fasta,./rRNA_databases/silva-euk-18s-id95:./rRNA_databases/silva-euk-28s-id98.fasta,./rRNA_databases/silva-euk-28s-id98:./rRNA_databases/rfam-5s-database-id98.fasta,./rRNA_databases/rfam-5s-database-id98:./rRNA_databases/rfam-5.8s-database-id98.fasta,./rRNA_databases/rfam-5.8s-database-id98 \
    --reads ${name}.merged.fastq \
    --paired_in --aligned ${name}.aligned \
    --other ${name}.other_merged \
    --fastx --log --num_alignments 1 -v \
    -a ${task.cpus}
    unmerge-paired-reads.sh ${name}.other_merged.fastq ${name}.R1.other.fastq ${name}.R2.other.fastq
    pigz -p ${task.cpus} ${name}.R1.other.fastq
    pigz -p ${task.cpus} ${name}.R2.other.fastq
    rm ${name}.merged.fastq ${name}.aligned.fastq ${name}.other_merged.fastq ${uncompr_reads_R1} ${uncompr_reads_R2}
    """
    }
}