process get_reduced_genome_test {
    label 'basic_tools'
    label 'smallTask'
    
    if (params.cloudProcess) { publishDir "${params.permanentCacheDir}/genomes/", mode: 'copy', pattern: "*.fa" }
    else { storeDir "${params.permanentCacheDir}/genomes/" }  
    
    input:
    val(species)

    output:
    path("${species}_small.fa")

    script:
    if (species == 'hsa') 
        """
        wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        gunzip -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 

        NSEQS=3
        awk "/^>/ {n++} n>\$NSEQS {exit} {print}" Homo_sapiens.GRCh38.dna.primary_assembly.fa > ${species}_small.fa

        rm Homo_sapiens.GRCh38.dna.primary_assembly.fa
        """
    else
        error "Invalid species: ${species} for test profile."
}

process reduce_genome_test{
    label 'basic_tools'
    label 'smallTask'

    if (params.cloudProcess) { publishDir "${params.permanentCacheDir}/genomes/", mode: 'copy', pattern: "*.fa" }
    else { storeDir "${params.permanentCacheDir}/genomes/" }  

    input:
    path(complete_genome)

    output:
    path("${complete_genome.baseName}_small.fa")

    script:
    """
    NSEQS=3
    awk "/^>/ {n++} n>\$NSEQS {exit} {print}" ${complete_genome} > ${complete_genome.baseName}_small.fa
    """

}

process get_reduced_annotation_test {
    label 'basic_tools'
    label 'smallTask'

    if (params.cloudProcess) { publishDir "${params.permanentCacheDir}/annotations/", mode: 'copy', pattern: "*.gtf" }
    else { storeDir "${params.permanentCacheDir}/annotations/" }  
    
    input:
    val(species)

    output:
    path("${species}_small.gtf")

    script:
    if (species == 'hsa') 
        """
        wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
        gunzip -f Homo_sapiens.GRCh38.98.gtf.gz

        grep -P '^1[01]?\\s' Homo_sapiens.GRCh38.98.gtf > ${species}_small.gtf

        rm Homo_sapiens.GRCh38.98.gtf
        """
    else
    error "Invalid species: ${species} for test profile."

}

process reduce_annotation_test{
    label 'basic_tools'
    label 'smallTask'

    if (params.cloudProcess) { publishDir "${params.permanentCacheDir}/annotations/", mode: 'copy', pattern: "*.gtf" }
    else { storeDir "${params.permanentCacheDir}/annotations/" }  

    input:
    path(complete_annotation)

    output:
    path("${complete_annotation.baseName}_small.gtf")

    script:
    """
    grep -P '^1[01]?\\s' ${complete_annotation} > ${complete_annotation.baseName}_small.gtf
    """

}