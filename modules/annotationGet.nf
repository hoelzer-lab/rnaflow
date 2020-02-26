/*******************************************
* DOWNLOAD ANNOTATION
********************************************/
process annotationGet {
    label 'python3'
    
    if (params.cloudProcess) { publishDir "${params.cloudDatabase}/annotations/${params.species}", mode: 'copy', pattern: "${params.species}.gtf" }
    else { storeDir "nextflow-autodownload-databases/annotations/${params.species}" }  
    
    output:
    file("${params.species}.gtf")
    
    script:
    if (params.species == 'hsa') {
    """
    wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
    gunzip -f Homo_sapiens.GRCh38.98.gtf.gz
    mv Homo_sapiens.GRCh38.98.gtf ${params.species}.gtf
    """
    }
    if (params.species == 'eco') {
    """
    wget ftp://ftp.ensemblgenomes.org/pub/release-45/bacteria//gtf/bacteria_90_collection/escherichia_coli_k_12/Escherichia_coli_k_12.ASM80076v1.45.gtf.gz
    gunzip -f Escherichia_coli_k_12.ASM80076v1.45.gtf.gz
    mv Escherichia_coli_k_12.ASM80076v1.45.gtf ${params.species}.gtf
    """
    }
}



