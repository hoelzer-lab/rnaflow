/*******************************************
* DOWNLOAD SORTMERNA DATABASES
********************************************/
process sortmernaGET {
    //conda 'envs/hisat2.yaml'
    if (params.cloudProcess) { publishDir "${params.cloudDatabase}/databases/XXXXXXXX", mode: 'copy', pattern: "XXXXXXXX" }
    else { storeDir "nextflow-autodownload-databases/databases/XXXXXXXX" }  
    output:
        file("XXXXXXXX")
    script:
        """
        wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
        gunzip -f Homo_sapiens.GRCh38.98.gtf.gz
        mv Homo_sapiens.GRCh38.98.gtf ${params.annotation}.gtf
        """
}



