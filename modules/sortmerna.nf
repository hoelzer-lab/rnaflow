/************************************************************************
* SORTMERNA
*
* Remove rRNA reads
************************************************************************/
process sortmerna {
    label 'sortmerna'
    tag "$meta.sample"

    if ( params.nanopore ) { errorStrategy 'ignore' }

    if ( params.softlink_results ) { publishDir "${params.output}/${params.sortmerna_dir}", pattern: "*.other.fastq.gz" }
    else { publishDir "${params.output}/${params.sortmerna_dir}", mode: 'copy', pattern: "*.other.fastq.gz" }

    input:
    tuple val(meta), path(reads)
    path(db)

    output:
    tuple val(meta), path("${meta.sample}*.other.fastq.gz"), emit: no_rna_fastq
    path "${meta.sample}.aligned.log", emit: log

    script:
    if (!meta.paired_end) {
    assert "$reads".split('.gz').size() == 1 : "suffix should be .gz"
    uncompr_reads = "$reads".split('.gz')[0]
    """
	unpigz -f -p ${task.cpus} ${reads[0]}
	tar -zxvf ${db}

	# Database Options (sorted by type and how fast they are)
	# smr_v4.3_fast_db.fasta / smr_v4.3_default_db.fasta
	# smr_v4.3_sensitive_db_rfam_seeds.fasta / smr_v4.3_sensitive_db.fasta 
	
	sortmerna --ref smr_v4.3_default_db.fasta \
	    --reads ${uncompr_reads} \
	    --aligned ${meta.sample}.aligned --other ${meta.sample}.other \
	    --fastx \
            --num_alignments 1 \
	    -threads ${task.cpus}

	pigz -p ${task.cpus} ${meta.sample}.other.fq

	# Name aligment to a previous schema
	mv ${meta.sample}.other.fq.gz ${meta.sample}.other.fastq.gz

	# Cleanup
	rm ${meta.sample}.aligned.fq ${uncompr_reads} *fasta #dbs
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
	tar -zxvf ${db}

	# Database Options (sorted by type and how fast they are)
	# smr_v4.3_fast_db.fasta / smr_v4.3_default_db.fasta
	# smr_v4.3_sensitive_db_rfam_seeds.fasta / smr_v4.3_sensitive_db.fasta 

	sortmerna --ref smr_v4.3_default_db.fasta \
	    --reads ${uncompr_reads_R1} --reads ${uncompr_reads_R2} \
	    --paired_in --out2 --aligned ${meta.sample}.aligned --other ${meta.sample}.other_merged \
	    --fastx \
	    --num_alignments 1 \
	    -threads ${task.cpus}
	
	pigz -p ${task.cpus} ${meta.sample}.other_merged_fwd.fq
	pigz -p ${task.cpus} ${meta.sample}.other_merged_rev.fq

	# Name aligment to a previous schema
	mv ${meta.sample}.other_merged_fwd.fq.gz ${meta.sample}.R1.other.fastq.gz
	mv ${meta.sample}.other_merged_rev.fq.gz ${meta.sample}.R2.other.fastq.gz

	# Cleanup
	rm ${meta.sample}.aligned_fwd.fq ${meta.sample}.aligned_rev.fq ${uncompr_reads_R1} ${uncompr_reads_R2} *fasta #dbs
    """
    }
}
