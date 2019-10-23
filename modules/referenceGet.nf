/*******************************************
* DOWNLOAD A REFERENCE GENOME
********************************************/
process referenceGet {
    conda 'envs/hisat2.yaml'
    if (params.cloudProcess) { publishDir "${params.cloudDatabase}/genomes/${params.reference}", mode: 'copy', pattern: "${params.reference}.fa" }
    else { storeDir "nextflow-autodownload-databases/genomes/${params.reference}" }  
    output:
        file("${params.reference}.fa")
    script:
        if (params.reference == 'hsa') {
        """
        wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        gunzip -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 
        mv Homo_sapiens.GRCh38.dna.primary_assembly.fa ${params.reference}.fa
        """
        }
        if (params.reference == 'eco') {
        """
        wget ftp://ftp.ensemblgenomes.org/pub/release-45/bacteria//fasta/bacteria_90_collection/escherichia_coli_k_12/dna/Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa.gz
        gunzip -f Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa.gz
        mv Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa ${params.reference}.fa
        """
        }
}



