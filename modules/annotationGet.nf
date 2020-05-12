/*******************************************
* DOWNLOAD ANNOTATION
********************************************/
process annotationGet {
    label 'python3'
    
    if (params.cloudProcess) { publishDir "${params.cloudDatabase}/annotations/", mode: 'copy', pattern: "*.gtf" }
    else { storeDir "nextflow-autodownload-databases/annotations/" }  
    
    input:
    val(species)

    output:
    path("${species}.gtf")
    
    script:
    if (species == 'hsa') {
    """
    wget ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
    gunzip -f Homo_sapiens.GRCh38.98.gtf.gz
    mv Homo_sapiens.GRCh38.98.gtf ${species}.gtf
    """
    }
    if (species == 'eco') {
    """
    wget ftp://ftp.ensemblgenomes.org/pub/release-45/bacteria//gtf/bacteria_90_collection/escherichia_coli_k_12/Escherichia_coli_k_12.ASM80076v1.45.gtf.gz
    gunzip -f Escherichia_coli_k_12.ASM80076v1.45.gtf.gz
    mv Escherichia_coli_k_12.ASM80076v1.45.gtf ${species}.gtf
    """
    }
}

process concat_annotation {

  input:
  path '*'

  output:
  path 'annotation.gtf'
  
  script:
  """
  cat * > annotation.gtf
  """
}