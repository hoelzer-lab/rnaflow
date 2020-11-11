/*******************************************
* DOWNLOAD A REFERENCE GENOME
********************************************/
process referenceGet {
    label 'python3'
    label 'smallTask'
    
    if (params.cloudProcess) { publishDir "${params.permanentCacheDir}/genomes/", mode: 'copy', pattern: "*.fa" }
    else { storeDir "${params.permanentCacheDir}/genomes/" }  
    
    input:
    val(species)

    output:
    path("${species}.fa")
    
    script:
    if (species == 'hsa')
      """
      wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
      gunzip -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz 
      mv Homo_sapiens.GRCh38.dna.primary_assembly.fa ${species}.fa
      """
    else if (species == 'mmu')
      """
      wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
      gunzip -f Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
      mv Mus_musculus.GRCm38.dna.primary_assembly.fa ${species}.fa
      """
    else if (species == 'eco')
      """
      wget ftp://ftp.ensemblgenomes.org/pub/release-45/bacteria//fasta/bacteria_90_collection/escherichia_coli_k_12/dna/Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa.gz
      gunzip -f Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa.gz
      mv Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa ${species}.fa
      """
    else if (species == 'mau')
      """
      wget ftp://ftp.ensembl.org/pub/release-100/fasta/mesocricetus_auratus/dna/Mesocricetus_auratus.MesAur1.0.dna.toplevel.fa.gz
      gunzip Mesocricetus_auratus.MesAur1.0.dna.toplevel.fa.gz
      mv Mesocricetus_auratus.MesAur1.0.dna.toplevel.fa ${species}.fa
      """
    else
      error "Invalid species: ${species}"
}

process concat_genome {
  label 'smallTask'
  
  input:
  path '*'

  output:
  path 'reference.fa'
  
  script:
  """
  cat * > reference.fa
  """
}