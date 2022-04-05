/*******************************************
* DOWNLOAD SORTMERNA DATABASES
********************************************/
process sortmernaGet {
    label 'sortmerna'
    if (!params.cloudProcess) { label 'smallTask' }

    if (params.cloudProcess) { publishDir "${params.permanentCacheDir}/databases/sortmerna/", mode: 'copy', pattern: "rRNA_databases.tar.gz" }
    else { storeDir "${params.permanentCacheDir}/databases/sortmerna/" }  

    output:
    path("rRNA_databases.tar.gz")

    script:
    """
    git clone https://github.com/biocore/sortmerna.git
    indexdb_rna --ref ./sortmerna/data/rRNA_databases/silva-bac-16s-id90.fasta,./sortmerna/data/rRNA_databases/silva-bac-16s-id90:./sortmerna/data/rRNA_databases/silva-bac-23s-id98.fasta,./sortmerna/data/rRNA_databases/silva-bac-23s-id98:./sortmerna/data/rRNA_databases/silva-arc-16s-id95.fasta,./sortmerna/data/rRNA_databases/silva-arc-16s-id95:./sortmerna/data/rRNA_databases/silva-arc-23s-id98.fasta,./sortmerna/data/rRNA_databases/silva-arc-23s-id98:./sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta,./sortmerna/data/rRNA_databases/silva-euk-18s-id95:./sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta,./sortmerna/data/rRNA_databases/silva-euk-28s-id98:./sortmerna/data/rRNA_databases/rfam-5s-database-id98.fasta,./sortmerna/data/rRNA_databases/rfam-5s-database-id98:./sortmerna/data/rRNA_databases/rfam-5.8s-database-id98.fasta,./sortmerna/data/rRNA_databases/rfam-5.8s-database-id98 -v
    mv sortmerna/data/rRNA_databases . 
    tar zcvf rRNA_databases.tar.gz rRNA_databases
    rm -r rRNA_databases
    rm -r sortmerna
    """
    }



