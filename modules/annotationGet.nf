/*******************************************
* DOWNLOAD ANNOTATION
********************************************/
process annotationGet {
    label 'basic_tools'
    if (!params.workflow.contains('node')) { label 'smallTask' }
    
    if (params.cloudProcess) { publishDir "${params.permanentCacheDir}/annotations/", mode: 'copy', pattern: "*.gtf" }
    else { storeDir "${params.permanentCacheDir}/annotations/" }  
    
    input:
    val(species)

    output:
    path("${species}.gtf")
    
    script:
    if (species == 'hsa') 
      """
      wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
      gunzip -f Homo_sapiens.GRCh38.98.gtf.gz
      mv Homo_sapiens.GRCh38.98.gtf ${species}.gtf
      """
    else if (species == 'mmu')
      """
      wget ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
      gunzip -f Mus_musculus.GRCm38.99.gtf.gz
      mv Mus_musculus.GRCm38.99.gtf ${species}.gtf
      """
    else if (species == 'eco')
      """
      wget ftp://ftp.ensemblgenomes.org/pub/release-45/bacteria//gtf/bacteria_90_collection/escherichia_coli_k_12/Escherichia_coli_k_12.ASM80076v1.45.gtf.gz
      gunzip -f Escherichia_coli_k_12.ASM80076v1.45.gtf.gz
      mv Escherichia_coli_k_12.ASM80076v1.45.gtf ${species}.gtf
      """
    else if (species == 'mau') 
      """
      wget ftp://ftp.ensembl.org/pub/release-100/gtf/mesocricetus_auratus/Mesocricetus_auratus.MesAur1.0.100.gtf.gz
      gunzip -f Mesocricetus_auratus.MesAur1.0.100.gtf.gz
      mv Mesocricetus_auratus.MesAur1.0.100.gtf ${species}.gtf
      """
    else
      error "Invalid species: ${species}"
}

process concat_annotation {
  label 'basic_tools'
  if (!params.workflow.contains('node')) { label 'smallTask' }

  input:
  path annotation

  output:
  path 'annotation.gtf'
  
  script:
  """
  cat ${annotation} > annotation.gtf
  """
}