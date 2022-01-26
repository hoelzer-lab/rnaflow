library("piano")
library("biomaRt")
library("openxlsx")
library("snowfall")

# read arguments and evaluate
args <- commandArgs(TRUE)

project_dir <- eval( parse(text=args[1]) )[1] 
gene_file <- eval( parse(text=args[2]) )
species <- eval( parse(text=args[3]) )
id_type <- eval( parse(text=args[4]) )
cpus <- eval( parse(text=args[5]) )

resFold05 <- read.csv(file = gene_file, row.names = 1)


##########################################
## BiomaRt object
##########################################
try.biomart <- try(
  if (species == 'mmu'){
    biomart.ensembl <- useMart('ensembl', dataset='mmusculus_gene_ensembl')
  } else if (species == 'hsa') {
    biomart.ensembl <- useMart('ensembl', dataset='hsapiens_gene_ensembl')
  } else if (species == 'mau') {
    biomart.ensembl <- useMart('ensembl', dataset='mauratus_gene_ensembl')
  } else {
    biomart.ensembl <- NA
    print('SKIPPING: BiomaRt. Species not accasible with BiomaRt.')
  }
)
if (class(try.biomart) == "try-error") {
  biomart.ensembl <- NA
  print('SKIPPING: BiomaRt. BiomaRt is not accessible.')
}


##########################################
## Functions
##########################################
write.table.to.file <- function(as.data.frame.object, output.path, output.name, id2name, row.names=TRUE, col.names=TRUE) {
  output.file.basename <- paste0(output.path, "/", output.name)
  write.table(as.data.frame.object, file=paste0(output.file.basename, ".csv"), sep = ",", row.names=row.names, col.names=col.names)
  if( is.na(col.names) ){
    write.xlsx(as.data.frame.object, file=paste0(output.file.basename, ".xlsx"), row.names=row.names, col.names=TRUE, asTable=TRUE)
  } else {
    write.xlsx(as.data.frame.object, file=paste0(output.file.basename, ".xlsx"), row.names=row.names, col.names=col.names, asTable=TRUE)
  }

  if ( !missing(id2name)) {
    output.file.basename.extended <- paste0(output.path, "/", output.name, "_extended")
    ## add real gene names and biotypes to the csv files
    system(paste("./improve_deseq_table.rb", paste0(output.file.basename.extended, ".csv" ), paste0(output.file.basename, ".csv"), id2name, sep=" "), wait=TRUE)
    write.xlsx(read.csv(paste0(output.file.basename.extended, ".csv" )), paste0(output.file.basename.extended, ".xlsx" ), asTable=TRUE)
  }
}


piano <- function(out.dir, resFold, mapGO, cpus) {
  mapGO <- mapGO[mapGO[,2]!="",]
  write.table.to.file(mapGO, out.dir, "ENSG_GOterm", row.names = FALSE)
  
  myGsc <- loadGSC(mapGO)
  
  myPval <- resFold$padj
  names(myPval) <- rownames(resFold)
  myFC <- resFold$log2FoldChange
  names(myFC) <- rownames(resFold)
  
  if (cpus >= 10) {
    piano_cpus = 10
  } else {
    piano_cpus = 1
  }
  
  gene.set.min <- 20
  gene.set.max <- 'inf' # 9999999999999
  gsaRes1 <- runGSA(myFC, geneSetStat="maxmean", gsc=myGsc,
                    gsSizeLim=c(gene.set.min,gene.set.max), ncpus=piano_cpus)
  gsaRes2 <- runGSA(myFC, geneSetStat="gsea", gsc=myGsc,
                    gsSizeLim=c(gene.set.min,gene.set.max), ncpus=piano_cpus)
  gsaRes3 <- runGSA(myFC, geneSetStat="fgsea", gsc=myGsc,
                    gsSizeLim=c(gene.set.min,gene.set.max), ncpus=piano_cpus)
  gsaRes4 <- runGSA(myFC, geneSetStat="page", gsc=myGsc,
                    gsSizeLim=c(gene.set.min,gene.set.max), ncpus=piano_cpus)
  gsaRes5 <- runGSA(myPval, myFC, geneSetStat="fisher", gsc=myGsc,
                    gsSizeLim=c(gene.set.min,gene.set.max), ncpus=piano_cpus)
  gsaRes6 <- runGSA(myPval, myFC, geneSetStat="stouffer", gsc=myGsc,
                    gsSizeLim=c(gene.set.min,gene.set.max), ncpus=piano_cpus)
  gsaRes7 <- runGSA(myPval, myFC, geneSetStat="reporter", gsc=myGsc,
                    gsSizeLim=c(gene.set.min,gene.set.max), ncpus=piano_cpus)
  gsaRes8 <- runGSA(myPval, myFC, geneSetStat="tailStrength", gsc=myGsc,
                    gsSizeLim=c(gene.set.min,gene.set.max), ncpus=piano_cpus)
  
  resList <- list(gsaRes1,gsaRes2,gsaRes3,gsaRes4,gsaRes5,gsaRes6,gsaRes7,gsaRes8)
  names(resList) <- c("maxmean", "gsea", "fgsea", "page", "fisher", "stouffer", "reporter", "tailStrength")
  
  try.piano <- try( {
    pdf(paste(out.dir,"/consensus_heatmap.pdf",sep=""), width = 10, height = 10)
    ch <- consensusHeatmap(resList,cutoff=10,method="mean",colorkey=FALSE,cellnote="consensusScore",ncharLabel = 120) ## medianPvalue or consensusScore or nGenes
    dev.off()
    svg(paste(out.dir,"/consensus_heatmap.svg",sep=""), width = 10, height = 10)
    ch <- consensusHeatmap(resList,cutoff=10,method="mean",colorkey=FALSE,cellnote="consensusScore",ncharLabel = 120) ## medianPvalue or consensusScore
    dev.off()
    
    downregulated_paths <- as.data.frame(ch$pMat[,1][ch$pMat[,1] < 0.05])
    upregulated_paths <- as.data.frame(ch$pMat[,5][ch$pMat[,5] < 0.05])
    
    write.table.to.file(downregulated_paths, out.dir, "paths_sigdown", col.names=FALSE)
    write.table.to.file(upregulated_paths, out.dir, "paths_sigup", col.names=FALSE)
  }
  )
  if (class(try.piano) == "try-error") {
    print('SKIPPING: piano consensusHeatmap.')
  }
}

#####################
## Piano
#####################
if ( ! is.na(biomart.ensembl) ) {
  dir.create(file.path('downstream_analysis/piano'), showWarnings = FALSE, recursive = TRUE)
  if (any(grepl(id_type, listAttributes(biomart.ensembl)$name, fixed=TRUE))){
    results.gene <- getBM(attributes =  c(id_type, "name_1006"), filters = id_type, values = rownames(resFold05), mart=biomart.ensembl)
    if ( length(row.names(results.gene)) > 0 ) {
      try.piano <- try( 
        piano(paste('downstream_analysis', 'piano', sep='/'), resFold05, results.gene, cpus)
      )
      if (class(try.piano) == "try-error") {
        print ('SKIPPING: Piano. Some error occurred.')
      }
    } else {
      print(paste('SKIPPING: Piano. No matching feature IDs with type', id_type, 'found.'))
    }
  } else {
    print(paste('SKIPPING: Piano. Feature ID type', id_type, 'not supported by biomaRt.'))
  }
}