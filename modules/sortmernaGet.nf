/*******************************************
* DOWNLOAD SORTMERNA DATABASES
********************************************/
process sortmernaGet {
    label 'sortmerna'
    if (!workflow.profile.contains('node')) { label 'smallTask' }

    if (params.cloudProcess) { publishDir "${params.permanentCacheDir}/databases/sortmerna/", mode: 'copy', pattern: "rRNA_databases.tar.gz" }
    else { storeDir "${params.permanentCacheDir}/databases/sortmerna/" }  

    output:
    path("database.tar.gz")

    script:
    """
	# Database Files (sorted by type and how fast they are)
	# smr_v4.3_fast_db.fasta / smr_v4.3_default_db.fasta
	# smr_v4.3_sensitive_db_rfam_seeds.fasta / smr_v4.3_sensitive_db.fasta 
    # See for details: https://github.com/sortmerna/sortmerna?tab=readme-ov-file#databases

    wget https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
    """
    }



