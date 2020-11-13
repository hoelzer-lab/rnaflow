library("DESeq2")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("ReportingTools")
library("pheatmap")
library("biomaRt")
library("svglite")
library("piano")
library("apeglm")
library("EnhancedVolcano")
library("regionReport")
library("stringr")
library("WebGestaltR")
library("snowfall")
library("openxlsx")


#####################################################################################
## FUNCTIONS
#####################################################################################
build.project.structure <- function(out.dir) {
  dir.create(file.path(out.dir), showWarnings = FALSE)
  # Build necessary project structure
  dir.create(file.path(out.dir, '/results'), showWarnings = FALSE)
  dir.create(file.path(out.dir, '/reports'), showWarnings = FALSE)
  dir.create(file.path(out.dir, '/input'), showWarnings = FALSE)
  for (plot.type in c('volcano', 'PCA', 'heatmaps', 'MA', 'sample2sample')) {
    dir.create(file.path(out.dir, paste0('/plots/', plot.type)), showWarnings = FALSE, recursive = TRUE)
  }
}

write.table.to.file <- function(as.data.frame.object, output.path, output.name, ensembl2genes, row.names=TRUE, col.names=TRUE) {
  output.file.basename <- paste0(output.path, "/", output.name)
  write.table(as.data.frame.object, file=paste0(output.file.basename, ".csv"), sep = ",", row.names=row.names, col.names=col.names)
  if( is.na(col.names) ){
    write.xlsx(as.data.frame.object, file=paste0(output.file.basename, ".xlsx"), row.names=row.names, col.names=TRUE, asTable=TRUE)
  } else {
    write.xlsx(as.data.frame.object, file=paste0(output.file.basename, ".xlsx"), row.names=row.names, col.names=col.names, asTable=TRUE)
  }

  if ( !missing(ensembl2genes)) {
    output.file.basename.extended <- paste0(output.path, "/", output.name, "_extended")
    ## add real gene names and biotypes to the csv files
    system(paste("./improve_deseq_table.rb", paste0(output.file.basename.extended, ".csv" ), paste0(output.file.basename, ".csv"), ensembl2genes, sep=" "), wait=TRUE)
    write.xlsx(read.csv(paste0(output.file.basename.extended, ".csv" )), paste0(output.file.basename.extended, ".xlsx" ), asTable=TRUE)
  }
}

plot.sample2sample <- function(out.dir, col.labels, trsf_data, trsf_type, colors) {
  ## get sample-to-sample distances
  sampleDists <- dist(t(assay(trsf_data)))
  sampleDistMatrix <- as.matrix(sampleDists)
  ## add names
  rownames(sampleDistMatrix) <- with(colData(trsf_data), col.labels)
  colnames(sampleDistMatrix) <- with(colData(trsf_data), col.labels)


  pdf(paste(out.dir, paste0("sample2sample_", trsf_type, ".pdf"), sep="/"))
  pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, color = colors)
  dev.off()
  svg(paste(out.dir, paste0("sample2sample_", trsf_type, ".svg"), sep="/"))
  pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, color = colors)
  dev.off()

  # sample2sample heatmap with color key and histogram
  # hc <- hclust(sampleDists)
  # heatmap.2(sampleDistMatrix, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col = colors, margin=c(13, 13))
}

plot.ma <- function(output.dir, deseq2.res, alpha) {
  pdf(paste(output.dir, paste0("MA_alpha", alpha, ".pdf"), sep="/"))
  plotMA(deseq2.res, alpha = alpha, main = paste('MA plot with alpha =', alpha))
  dev.off()
  svg(paste(output.dir, paste0("MA_alpha", alpha, ".svg"), sep="/"))
  plotMA(deseq2.res, alpha = alpha, main = paste('MA plot with alpha =', alpha))
  dev.off()
}

reportingTools.html <- function(out.dir, dds, deseq2.result, pvalueCutoff, condition1, condition2, annotation_genes, make.plots=TRUE) {
  # Exporting results to HTML and CSV
  if (pvalueCutoff == 1.1){
    shortName <- 'RNAseq_analysis_with_DESeq2_full'
    title <- paste0('RNA-seq analysis of differential expression using DESeq2, no P value cutoff')
  } else {
    shortName <- paste0('RNAseq_analysis_with_DESeq2_p', pvalueCutoff)
    title <- paste0('RNA-seq analysis of differential expression using DESeq2, P value cutoff ', pvalueCutoff)
  }
  if (make.plots == FALSE) {
    dir.create(file.path(paste0(out.sub, '/reports/figures', shortName)), showWarnings = FALSE)
    for ( id in rownames(deseq2.result[ !is.na(deseq2.result$padj) & deseq2.result$padj < pvalueCutoff, ]) ) {
      system(paste0('cp ', out.dir, '/reports/figuresRNAseq_analysis_with_DESeq2_full/boxplot.', id, '.pdf ', out.sub, '/reports/figures', shortName))
      system(paste0('cp ', out.dir, '/reports/figuresRNAseq_analysis_with_DESeq2_full/mini.', id, '.png ', out.sub, '/reports/figures', shortName))
    }
  }
  des2Report <- HTMLReport(shortName=shortName, title=title, basePath=out.dir, reportDirectory="reports/")
  publish(dds, des2Report, pvalueCutoff=pvalueCutoff, annotation.db=NULL, factor=colData(dds)$condition, reportDir=out.dir, n=length(row.names(deseq2.result)), contrast=c("condition",condition1,condition2), make.plots=make.plots)
  finish(des2Report)
  system(paste('./refactor_reportingtools_table.rb', paste0(out.dir, '/reports/', shortName,'.html'), annotation_genes, 'add_plots', sep=" "))
}

plot.pca <- function(out.dir, col.labels, trsf_data, trsf_type, ntop) {
  # calculate the variance for each gene
  rv <- rowVars(assay(trsf_data))
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  # Extract the data 
  X <- t(assay(trsf_data)[select,]) # Transpose this as our read count matrix as R has dimensions as columns and not as rows (thanks, R!!!)
  
  # Using R's internal function for improved speed (and accuracy as they use SDV)
  # Caution: R will not consider all eigenvectors (there are thousands of genes)
  # Theoretically, we need to calculate ALL of them (we then obtain PC1, PC2, ... PCm with m dimensions = genes)
  # But R will truncate it to PC1, PC2, ... PCn with n data points (if n < m), which is fast.
  # Don't let this confuse you
  pca <- prcomp(X, center = TRUE, scale = FALSE) # default: center = TRUE, scale = FALSE
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  intgroup <- c("condition")
  if (!all(intgroup %in% names(colData(trsf_data)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(trsf_data)[, intgroup, drop=FALSE])

  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(trsf_data)[[intgroup]]
  }

  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=col.labels)

  ggplot(data=d, aes_string(x="PC1", y="PC2", colour="condition")) +
    geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    ggtitle(paste("PC1 vs PC2: top ", ntop, " variable genes")) +
    coord_fixed() +
    ggsave(paste(out.dir, paste0("PCA_simple_", trsf_type, "_top", ntop, ".pdf"), sep="/"))
    ggsave(paste(out.dir, paste0("PCA_simple_", trsf_type, "_top", ntop, ".svg"), sep="/"))
}

plot.heatmap.most_var <- function(out.dir, dds, trsf_data, trsf_type, ntop, samples.info=df.samples.info, genes.info=df.gene.anno) {
  select <- order(rowVars(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:ntop]
  selected.ensembl.ids <- row.names(trsf_data[select,])

  file <- paste(out.dir, paste0("heatmap_count_matrix_", trsf_type, "_mostVar", ntop, "_row-scaled.pdf"), sep="/")
  pheatmap(assay(trsf_data)[select,], cluster_cols = FALSE, cluster_rows = TRUE,
        scale = "row", border_color = NA,
        labels_row = as.character(genes.info[selected.ensembl.ids,]$gene_type),
        annotation_col=samples.info[ , !(colnames(samples.info) == 'columns'), drop=FALSE],
        labels_col = as.character(samples.info[colnames(trsf_data),]$columns),
        height = 12, width = 8, file = file)
}

plot.heatmap.top_counts <- function(out.dir, dds, trsf_data, trsf_type, ntop, samples.info=df.samples.info, genes.info=df.gene.anno) {
  select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:ntop]
  selected.ensembl.ids <- row.names(counts(dds,normalized=TRUE)[select,])

  file <- paste(out.dir, paste0("heatmap_count_matrix_", trsf_type, "_top", ntop, "_row-scaled.pdf"), sep="/")
  pheatmap(assay(trsf_data)[select,], cluster_cols = FALSE, cluster_rows = TRUE,
          scale = "row", border_color = NA,
          labels_row = as.character(genes.info[selected.ensembl.ids,]$gene_type),
          annotation_col=samples.info[ , !(colnames(samples.info) == 'columns'), drop=FALSE],
          labels_col = as.character(samples.info[colnames(trsf_data),]$columns),
          height = 12, width = 8, file = file)
}

plot.heatmap.top_fc <- function(out.dir, resFold, trsf_data, trsf_type, ntop, samples.info=df.samples.info, genes.info=df.gene.anno) {
  selected.ensembl.ids <- row.names(resFold[order(resFold$log2FoldChange, decreasing=TRUE), ])[1:ntop]
  
  file <- paste(out.dir, paste0("heatmap_count_matrix_", trsf_type, "_top", ntop, "log2FC_row-scaled.pdf"), sep="/")
  pheatmap(assay(trsf_data)[selected.ensembl.ids,], cluster_cols = FALSE, cluster_rows = TRUE, 
          scale = "row", border_color = NA, 
          labels_row = as.character(genes.info[selected.ensembl.ids,]$gene_type),
          annotation_col=samples.info[ , !(colnames(samples.info) == 'columns'), drop=FALSE],
          labels_col = as.character(samples.info[colnames(trsf_data),]$columns),
          height = 12, width = 8, file = file)
  }

piano <- function(out.dir, resFold, mapGO, cpus) {
  mapGO <- mapGO[mapGO[,2]!="",]
  write.table.to.file(mapGO, out.dir, "ENSG_GOterm", row.names = FALSE)

  myGsc <- loadGSC(mapGO)

  myPval <- resFold$padj
  names(myPval) <- rownames(resFold)
  myFC <- resFold$log2FoldChange
  names(myFC) <- rownames(resFold)

  gene.set.min <- 20
  gene.set.max <- 'inf' # 9999999999999
  gsaRes1 <- runGSA(myFC, geneSetStat="maxmean", gsc=myGsc,
                  gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
  gsaRes2 <- runGSA(myFC, geneSetStat="gsea", gsc=myGsc,
                  gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
  gsaRes3 <- runGSA(myFC, geneSetStat="fgsea", gsc=myGsc,
                  gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
  gsaRes4 <- runGSA(myFC, geneSetStat="page", gsc=myGsc,
                  gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
  gsaRes5 <- runGSA(myPval, myFC, geneSetStat="fisher", gsc=myGsc,
                  gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
  gsaRes6 <- runGSA(myPval, myFC, geneSetStat="stouffer", gsc=myGsc,
                  gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
  gsaRes7 <- runGSA(myPval, myFC, geneSetStat="reporter", gsc=myGsc,
                  gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
  gsaRes8 <- runGSA(myPval, myFC, geneSetStat="tailStrength", gsc=myGsc,
                  gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)

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

  # for (i in 1:length(resList)){
  #   svg(paste(out.dir, paste0(names(resList)[i], '.svg'), sep='/'), width = 10, height = 10)
  #   networkPlot(resList[[i]], class="non")
  #   dev.off()
  # }
}

##################### TODO
# plot.ma.go <- function(out, deseq2.res, ma.size, results.gene, go.terms, trsf_data, trsf_type) {
#   ## We can also make an MA-plot for the results table in which we raised
#   ## the log2 fold change threshold (Figure below). We can label individual
#   ## points on the MA-plot as well. Here we use the with R function to plot
#   ## a circle and text for a selected row of the results object. Within the
#   ## with function, only the baseMean and log2FoldChange values for the
#   ## selected rows of res are used.
#   ##-----------------------------
#   for (go.term.ma in go.terms) {
#     #go.term.ma <- "GO:0009615"
#     pdf(paste(out,"statistics/ma_",trsf_type,"_", gsub(":", "", go.term.ma), ".pdf",sep=""))
#     plotMA(deseq2.res, main=paste("DESeq2, ", go.term.ma, sep=''), ylim=ma.size)
#     results.gene.GO.ma <- grep(go.term.ma, results.gene$go_id, fixed=TRUE)  ### e.g. GO:0002376, immune system process in mice
#     trsf_data.go.ma <- rownames(assay(trsf_data)[results.gene[results.gene.GO.ma,]$ensembl_gene_id,]) # get the ensembl ids corresponding to this go term
#     for (gene in trsf_data.go.ma) {
#       index = which(ensembl.ids == gene)
#       gene.name <- toString(gene.ids[index])
#       with(deseq2.res[gene, ], {
#         if (gene %in% rownames(resFold05)) {
#           points(baseMean, log2FoldChange, col="dodgerblue", cex=0.8, lwd=2, bg="dodgerblue")
#           text(baseMean, log2FoldChange, gene.name, pos=2, col="dodgerblue")
#         }
#       })
#     }
#     dev.off()
#   }
# }
#####################################################################################
## END FUNCTIONS
#####################################################################################


#####################################################################################
## MAIN 
#####################################################################################

##########################################
## Preparation
##########################################

#####################
## Parse arguments

args <- commandArgs(TRUE) # Read the arguments passed from the command line and assigns them to a vector of characters

## Parse the arguments (in characters) and evaluate them
project_dir <- eval( parse(text=args[1]) )[1] 
samples <- eval( parse(text=args[2]) )
conditions <- eval( parse(text=args[3]) )
col.labels <- eval( parse(text=args[4]) )
levels <- eval( parse(text=args[5]) )
comparisons <- eval( parse(text=args[6]) )
ensembl2genes <- eval( parse(text=args[7]) )[1]
annotation_genes <- eval( parse(text=args[8]) )[1]
patients <- eval( parse(text=args[9]) )
species <- eval( parse(text=args[10]) )
regionReport_config  <- eval( parse(text=args[11]) )[1]
regionReport_config <- normalizePath(regionReport_config) # regionReport needs the absolute path
cpus <- eval( parse(text=args[12]) )
#go.terms <- c()
#go.terms <- eval( parse(text=args[12]) ) # c("GO:004563","GO:0011231",...)

#####################
## Read in ensembl ids, gene names and biotypes from a tab seperated table
df.gene.anno <- as.data.frame( read.table(ensembl2genes, header=FALSE, sep="\t") )
rownames(df.gene.anno) <- df.gene.anno$V1
df.gene.anno$V1 <- NULL
colnames(df.gene.anno) <- c('gene_symbol', 'biotype')
df.gene.anno$gene_type <- paste(df.gene.anno$gene_symbol, str_replace(df.gene.anno$biotype, '_', ' '), sep=', ')

#####################
## Build project structure
out <- paste(project_dir,'/',sep='') # deseq2 dir is created by nextflow in the results dir ()
dir.create(file.path(out), showWarnings = FALSE)
dir.create(file.path(out, 'plots'), showWarnings = FALSE)
dir.create(file.path(out, 'plots/PCA'), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out, 'plots/heatmaps'), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out, 'data/input'), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out, 'data/counts'), showWarnings = FALSE, recursive = TRUE)

##########################################
## DESeq2 stuff
##########################################

#####################
## Create input object
if (length(patients) > 0) {
    sampleTable <- data.frame(sampleName = samples, fileName = samples, condition = conditions, type = col.labels, patients = patients, design = paste(patients, conditions, sep = ':'))
    ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design= ~ patients + condition)
} else {
    sampleTable <- data.frame(sampleName = samples, fileName = samples, condition = conditions, type = col.labels, design = conditions)
    ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design= ~ condition)
}
## adjust DESeq2 object to avoid comparision in alphabetical order(!!!1!1):
## order factor() in every level according to input files
ddsHTSeq$condition <- factor(ddsHTSeq$condition, levels=levels)
ddsHTSeq$type <- factor(ddsHTSeq$type, levels=col.labels)

#####################
## Create DESeqDataSet object by running DESeq2 on the input object
dds <- DESeq(ddsHTSeq)

##########################################
## Write input
##########################################

## raw input
if (length(patients) > 0) {
  df.samples.info <- data.frame(samples = samples, columns = col.labels, conditions = conditions, patients = patients)
} else {
  df.samples.info <- data.frame(samples = samples, columns = col.labels, conditions = conditions)
}
rownames(df.samples.info) <- df.samples.info$samples
df.samples.info$samples <- NULL
write.table.to.file(df.samples.info, paste0(out, "/data/input"), "input", col.names=NA)

## DESeq2 input
input.summary <- paste(out, "/data/input/", "DESeq2_input_summary.txt", sep="/") 
cat("Count input object:\n", file=input.summary, append=TRUE)
sink(input.summary, append=TRUE)
print(ddsHTSeq)
sink()
cat("\n\nCondition of count input object:\n", file=input.summary, append=TRUE)
sink(input.summary, append=TRUE)
print(ddsHTSeq$condition)
sink()
cat("\n\nDESeqDataSet object:\n", file=input.summary, append=TRUE)
sink(input.summary, append=TRUE)
print(dds)
sink()

##########################################
## Normalization and transformation
##########################################

#####################
## normalize counts
norm.counts <- counts(dds, normalized=T)

#####################
## write normalized counts and size factors
write.table.to.file(as.data.frame(norm.counts), paste0(out, "/data/counts"), "normalized_counts", col.names=NA)
write.table.to.file(as.data.frame(dds$sizeFactor), paste0(out, "/data/counts"), "sizeFactors", col.names=NA)

#####################
## transform counts
rld <- rlog(dds, blind=FALSE) # regularized log transformation
try.vst <- try(
  vsd <- vst(dds, blind=FALSE) # variance stabilizing transformation (VST)
)
if (class(try.vst) == "try-error") {
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
}
ntd <- normTransform(dds) # log2(n + 1) transformation

#####################
## collect transformed counts for easy iterating
transformed.counts = vector(mode="list", length=3)
names(transformed.counts) = c("vsd", "rld", "ntd")
transformed.counts[[1]] <- vsd; transformed.counts[[2]] <- rld; transformed.counts[[3]] <- ntd

#####################
## write transformed counts 
for (i in 1:length(transformed.counts)) {
  write.table.to.file(as.data.frame(assay(transformed.counts[[i]])), paste0(out, "/data/counts"), paste0("transformed_counts_", names(transformed.counts)[[i]]), col.names=NA)
}

##########################################
## Visualisation
##########################################

#####################
## PCA
for (i in 1:length(transformed.counts)) { 
  for (ntop in c(500, 100, 50)){
    plot.pca(paste(out, "plots/PCA/", sep="/"), col.labels, transformed.counts[[i]], names(transformed.counts)[[i]], ntop)
  }
}

#####################
## Heatmaps on counts
for (i in 1:length(transformed.counts)) { 
  for (ntop in c(50, 100)){
    plot.heatmap.top_counts(paste(out, "plots/heatmaps/", sep="/"), dds, transformed.counts[[i]], names(transformed.counts)[[i]], ntop)
  }
}

#####################
## Heatmaps on counts, most variable transformed genes
for (i in 1:length(transformed.counts)) { 
  for (ntop in c(50, 100)){
    plot.heatmap.most_var(paste(out, "plots/heatmaps/", sep="/"), dds, transformed.counts[[i]], names(transformed.counts)[[i]], ntop)
  }
}

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

#####################################################################################
## PERFORM PAIRWISE COMPARISONS
#####################################################################################
for (comparison in comparisons) {

  l1 <- strsplit(comparison, ':')[[1]][1]
  l2 <- strsplit(comparison, ':')[[1]][2]

  out.sub <- paste(out, l1, '_vs_', l2, '/', sep='')
  build.project.structure(out.sub)

  name <- paste("deseq2_",l1,"_",l2,sep="")

  ##########################################
  ## Adjust data, count data, levles and valiables to current pairwise comparison
  ##########################################

  dds.sub <- dds[ , dds$condition %in% c(l1, l2) ]
  dds.sub$condition <- droplevels(dds.sub$condition)
  dds.sub$type <- droplevels(dds.sub$type)

  ## transformed count data
  rld.sub <- rld[ , rld$condition %in% c(l1, l2) ]
  rld.sub$condition <- droplevels(rld.sub$condition)
  rld.sub$type <- droplevels(rld.sub$type)

  vsd.sub <- vsd[ , vsd$condition %in% c(l1, l2) ]
  vsd.sub$condition <- droplevels(vsd.sub$condition)
  vsd.sub$type <- droplevels(vsd.sub$type)

  ntd.sub <- ntd[ , ntd$condition %in% c(l1, l2) ]
  ntd.sub$condition <- droplevels(ntd.sub$condition)
  ntd.sub$type <- droplevels(ntd.sub$type)

  transformed.counts.sub = vector(mode="list", length=3)
  names(transformed.counts.sub) = c("vsd", "rld", "ntd")
  transformed.counts.sub[[1]] <- vsd.sub; transformed.counts.sub[[2]] <- rld.sub; transformed.counts.sub[[3]] <- ntd.sub

  ## adjust variabels
  conditions.sub <- c()
  col.labels.sub <- c()
  samples.sub <- c()
  levels.sub <- c(l1, l2)
  for (pos in which(conditions == l1)) {
    conditions.sub <- c(conditions.sub, conditions[pos])
    col.labels.sub <- c(col.labels.sub, col.labels[pos])
    samples.sub <- c(samples.sub, samples[pos])
  }
  for (pos in which(conditions == l2)) {
    conditions.sub <- c(conditions.sub, conditions[pos])
    col.labels.sub <- c(col.labels.sub, col.labels[pos])
    samples.sub <- c(samples.sub, samples[pos])
  }

  ##########################################
  ## Perform the pairwise comparison 
  ## code inspiration: https://github.com/acidgenomics/DESeqAnalysis/blob/master/R/apeglmResults-methods.R#L95
  ##########################################

  factor <- "condition"
  numerator <- l2
  denominator <- l1
  group <- colData(dds)[[factor]]
  group <- relevel(x = group, ref = denominator)
  colData(dds)[[factor]] <- group
  dds <- DESeq(dds) # nbinomWaldTest() via DESeq(), but the dispersion does not have to be estimated again | was not done before
  resultsNames <- resultsNames(dds)
  coef <- match(
    x = paste(factor, numerator, "vs", denominator, sep = "_"),
    table = resultsNames )
  deseq2.res <- lfcShrink(
        dds = dds,
        type = "apeglm",
        coef = coef
  )

  ##########################################
  ## Order and filter output
  ##########################################

  ## ordered by smallest adjusted p value
  resOrdered <- deseq2.res[order(deseq2.res$padj), ]

  ## filter NA values in fc and padj
  resNA = deseq2.res[ !is.na(deseq2.res$log2FoldChange) , ]
  resNA = resNA[ !is.na(resNA$padj) , ]

  ## filter 0 baseMean
  resBaseMean = resNA[ resNA$baseMean > 0.0 , ]

  ## resFold is now sorted by abs(foldchange) and all NA entries are removed as well as all zero baseMean values
  resFold <<-resBaseMean[rev(order(abs(resBaseMean$log2FoldChange))),]

  ## filter for specific adjusted P value
  resFold05 <<- resFold[ resFold$padj < 0.05 , ]
  resFold01 <<- resFold[ resFold$padj < 0.01 , ]

  ##########################################
  ## Write input and output
  ##########################################

  #####################
  ## input
  df.samples.info.sub <- df.samples.info[samples.sub,]
  write.table.to.file(df.samples.info.sub,  paste0(out.sub, "/input"), "input", col.names=NA)

  #####################
  ## DESeq2 results
  out.sub.output.dir <- paste0(out.sub, "/results/")

  write.table.to.file(as.data.frame(resOrdered), out.sub.output.dir, paste(name, "full", sep="_"), ensembl2genes, col.names=NA)
  write.table.to.file(as.data.frame(resFold), out.sub.output.dir, paste(name, "filtered_NA", sep="_"), ensembl2genes, col.names=NA)
  write.table.to.file(as.data.frame(resFold05), out.sub.output.dir, paste(name, "filtered_padj_0.05", sep="_"), ensembl2genes, col.names=NA)
  write.table.to.file(as.data.frame(resFold01), out.sub.output.dir, paste(name, "filtered_padj_0.01", sep="_"), ensembl2genes, col.names=NA)

  #####################
  ## DESeq2 results summary
  summary <- paste(out.sub.output.dir,"summary.txt",sep="/")
  sink(summary)
  summary(deseq2.res)
  sink()
  cat("#deseq2.res$padj < 0.1:\nFALSE\tTRUE\n", file=summary, append=TRUE)
  cat(table(deseq2.res$padj < 0.1), file=summary, append=TRUE)
  cat("\n\n", file=summary, append=TRUE)
  cat("#deseq2.res$padj < 0.05:\nFALSE\tTRUE\n", file=summary, append=TRUE)
  cat(table(deseq2.res$padj < 0.05), file=summary, append=TRUE)
  cat("\n\n", file=summary, append=TRUE)
  cat("#deseq2.res$padj < 0.01:\nFALSE\tTRUE\n", file=summary, append=TRUE)
  cat(table(deseq2.res$padj < 0.01), file=summary, append=TRUE)
  cat("\n", file=summary, append=TRUE)
  
  ##########################################
  ## Plots
  ##########################################
  
  #####################
  ## Volcano plot
  deseq2.res.anno <- merge(as.data.frame(deseq2.res), df.gene.anno, by=0)
  rownames(deseq2.res.anno) <- deseq2.res.anno$Row.names
  volcano = EnhancedVolcano(deseq2.res.anno, lab = deseq2.res.anno$gene_symbol, x = 'log2FoldChange', y = 'padj', 
    legendLabels = c('NS', expression(Log[2]~FC), "adj. p-value", expression(adj.~p-value~and~log[2]~FC)))
  volcano + 
    ggsave(paste(out.sub,"/plots/volcano/volcano.svg", sep='/')) +
    ggsave(paste(out.sub,"/plots/volcano/volcano.pdf", sep='/'))

  #####################
  ## MA plots
  plot.ma(paste(out.sub, "/plots/MA/", sep="/"), deseq2.res, metadata(deseq2.res)$alpha)
  plot.ma(paste(out.sub, "/plots/MA/", sep="/"), deseq2.res, metadata(deseq2.res)$alpha / 2)

  #####################
  ## Heatmaps of sample2sample distances
  for (i in 1:length(transformed.counts.sub)) {
    plot.sample2sample(paste(out.sub, "/plots/sample2sample/", sep="/"), col.labels.sub, 
      transformed.counts.sub[[i]], names(transformed.counts.sub)[[i]], colorRampPalette( rev(brewer.pal(9, "Blues")) )(255))
  }

  #####################
  ## Heatmaps on counts, most variable transformed genes
  for (i in 1:length(transformed.counts.sub)) { 
    for (ntop in c(50, 100)){
      plot.heatmap.most_var(paste(out.sub, "plots/heatmaps/", sep="/"), dds.sub, transformed.counts.sub[[i]], names(transformed.counts.sub)[[i]], ntop)
    }
  }

  #####################
  ## Heatmaps on counts, top count genes
  for (i in 1:length(transformed.counts.sub)) { 
    for (ntop in c(50, 100)){
      plot.heatmap.top_counts(paste(out.sub, "plots/heatmaps/", sep="/"), dds.sub, transformed.counts.sub[[i]], names(transformed.counts.sub)[[i]], ntop)
    }
  }
  
  #####################
  ## Heatmaps on counts, top FC genes
  for (i in 1:length(transformed.counts.sub)) {
    for (ntop in c(50, 100)){
      plot.heatmap.top_fc(paste(out.sub, "plots/heatmaps/", sep="/"), resFold, transformed.counts.sub[[i]], names(transformed.counts.sub)[[i]], ntop)
    }
  }

  #####################
  ## PCA
  for (i in 1:length(transformed.counts.sub)) {
    for (ntop in c(500, 100, 50)) {
      plot.pca(paste(out.sub, "/plots/PCA/", sep="/"), col.labels.sub, transformed.counts.sub[[i]], names(transformed.counts.sub)[[i]], ntop)
    }
  }

  ##########################################
  ## Further analysis
  ##########################################

  #####################
  ## MA plots with genes colored by GO terms
  # go.terms <- eval( c("GO:004563","GO:0011231") )
  # if ( ! is.na(biomart.ensembl)  && length(go.terms) > 0) {
  #   results.gene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "go_id","name_1006"), filters = "ensembl_gene_id", values = rownames(deseq2.res), mart=biomart.ensembl)
  #   #   plot.ma.go(out.sub, deseq2.res, ma.size, results.gene, go.terms, transformed.counts.sub[[i]], names(transformed.counts.sub)[[i]])
  # }

  #####################
  ## Piano
  if ( ! is.na(biomart.ensembl) ) {
    dir.create(file.path(out.sub, '/downstream_analysis/piano'), showWarnings = FALSE, recursive = TRUE)
    results.gene <- getBM(attributes =  c("ensembl_gene_id", "name_1006"), filters = "ensembl_gene_id", values = rownames(resFold05), mart=biomart.ensembl)
    piano(paste(out.sub, 'downstream_analysis', 'piano', sep='/'), resFold05, results.gene, cpus)
  }

  #####################
  ## Webgestalt
  if ( species == 'hsa' ){
    organism <- "hsapiens"
  } else if (species == 'mmu') {
    organism <- "mmusculus"
  } else {
    organism <- NA
  }
  if (! is.na(organism)) {
    dir.create(file.path(out.sub, '/downstream_analysis/WebGestalt'), showWarnings = FALSE, recursive = TRUE)
    interestGene <- as.data.frame(resFold05)[, 'log2FoldChange', drop=FALSE]
    interestGene$id <- rownames(interestGene)
    rownames(interestGene) <- NULL
    colnames(interestGene) <- NULL
    interestGene <- interestGene[c(2,1)]
    webgestalt.out.dir <- paste(out.sub, "downstream_analysis", "WebGestalt", sep='/')
    try.webgestalt <- try(
      for (enrDB in c("geneontology_Biological_Process_noRedundant", "pathway_KEGG")){
        enrichResult <-WebGestaltR(enrichMethod="GSEA", organism=organism, enrichDatabase=enrDB, interestGene=interestGene, interestGeneType="ensembl_gene_id", collapseMethod="mean", minNum=10, maxNum=500, fdrMethod="BH", sigMethod="fdr", fdrThr=0.01, topThr=10, perNum=1000, isOutput=TRUE, outputDirectory=webgestalt.out.dir, projectName=paste0(l1, '_vs_', l2))
      }
    )
    if (class(try.webgestalt) == "try-error") {
      print('SKIPPING: WebGestaltR. The number of annotated IDs for all functional categories are not from 10 to 500 for the GSEA enrichment method.')
    }
  }

  ##########################################
  ## Reports
  ##########################################

  #####################
  ## regionReport report
  ## needs knitr version 1.29. Some bug seems not to be fixed in 1.30 Anaconda version
  ## set output
  report.project.name <- paste(l1, "vs", l2, sep=" ")
  report.dir <- paste(out.sub, "reports", sep="/")
  report.output <- paste0('DESeq2_results_exploration')

  ## create html
  report_html <- DESeq2Report(dds, project = report.project.name,
    intgroup = c('condition', 'type'), res = deseq2.res, template = regionReport_config,
    outdir = report.dir, output = report.output, theme = theme_bw())

  ## and also pfd
  try( report_pdf <- DESeq2Report(dds, project = report.project.name,
    intgroup = c('condition', 'type'), res = deseq2.res, template = regionReport_config,
    outdir = report.dir, output = report.output, theme = theme_bw(),
    output_format = 'pdf_document', device = 'pdf') )


  #####################
  ## ReportingTools
  reportingTools.html(out.sub, dds, deseq2.res, 1.1, l2, l1, annotation_genes)
  if (length(rownames(resFold05)) > 0) {
    reportingTools.html(out.sub, dds, deseq2.res, 0.05, l2, l1, annotation_genes, make.plots=FALSE)
  }
  if (length(rownames(resFold01)) > 0) {
    reportingTools.html(out.sub, dds, deseq2.res, 0.01, l2, l1, annotation_genes, make.plots=FALSE)
  }

}
#####################################################################################
## END PAIRWISE COMPARISONS
#####################################################################################

