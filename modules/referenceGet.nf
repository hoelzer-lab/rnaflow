/*******************************************
* DOWNLOAD A REFERENCE GENOME
********************************************/
process referenceGet {
    label 'python3'
    if (params.cloudProcess) { publishDir "${params.cloudDatabase}/genomes/${params.species}", mode: 'copy', pattern: "${params.species}.fa" }
    else { storeDir "nextflow-autodownload-databases/genomes/${params.species}" }  
    
    output:
    file("${params.species}.fa")
    
    script:
    if (params.species == 'hsa') {
    """
    wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    gunzip -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 
    mv Homo_sapiens.GRCh38.dna.primary_assembly.fa ${params.species}.fa
    """
    }
    if (params.species == 'eco') {
    """
    wget ftp://ftp.ensemblgenomes.org/pub/release-45/bacteria//fasta/bacteria_90_collection/escherichia_coli_k_12/dna/Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa.gz
    gunzip -f Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa.gz
    mv Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa ${params.species}.fa
    """
    }
}



