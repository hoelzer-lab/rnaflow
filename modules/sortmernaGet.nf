/*******************************************
* DOWNLOAD SORTMERNA DATABASES
********************************************/
process sortmernaGet {
    label 'sortmerna'
    label 'smallTask'

    if (params.cloudProcess) { publishDir "${params.permanentCacheDir}/databases/", mode: 'copy', pattern: "sortmerna/data/rRNA_databases" }
    else { storeDir "${params.permanentCacheDir}/databases/" }  

    output:
    path("sortmerna/data/rRNA_databases")

    script:
    """
    git clone https://github.com/biocore/sortmerna.git
    indexdb_rna --ref ./sortmerna/data/rRNA_databases/silva-bac-16s-id90.fasta,./sortmerna/data/rRNA_databases/silva-bac-16s-id90:./sortmerna/data/rRNA_databases/silva-bac-23s-id98.fasta,./sortmerna/data/rRNA_databases/silva-bac-23s-id98:./sortmerna/data/rRNA_databases/silva-arc-16s-id95.fasta,./sortmerna/data/rRNA_databases/silva-arc-16s-id95:./sortmerna/data/rRNA_databases/silva-arc-23s-id98.fasta,./sortmerna/data/rRNA_databases/silva-arc-23s-id98:./sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta,./sortmerna/data/rRNA_databases/silva-euk-18s-id95:./sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta,./sortmerna/data/rRNA_databases/silva-euk-28s-id98:./sortmerna/data/rRNA_databases/rfam-5s-database-id98.fasta,./sortmerna/data/rRNA_databases/rfam-5s-database-id98:./sortmerna/data/rRNA_databases/rfam-5.8s-database-id98.fasta,./sortmerna/data/rRNA_databases/rfam-5.8s-database-id98 -v
    """
    }



