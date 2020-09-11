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

#####################################################################################
## FUNCTIONS
#####################################################################################
build.project.structure <- function(out) {
  dir.create(file.path(out), showWarnings = FALSE)
  # Build necessary project structure
  dir.create(file.path(out, '/statistics'), showWarnings = FALSE)
  dir.create(file.path(out, '/heatmaps'), showWarnings = FALSE)
  dir.create(file.path(out, '/results'), showWarnings = FALSE)
  dir.create(file.path(out, '/input'), showWarnings = FALSE)
  for (plot.type in c('volcano', 'PCA', 'heatmaps', 'MA', 'sample2sample')) {
    dir.create(file.path(out, paste0('/plots/', plot.type)), showWarnings = FALSE, recursive = TRUE)
  }
}

write.table.to.file <- function(as.data.frame.object, output.path, output.name, ensembl2genes) {
  output.file.basename <- paste0(output.path, "/", output.name)
  write.csv(as.data.frame.object, file=paste0(output.file.basename, ".csv"))
  system(paste("./csv_to_excel.py", paste0(output.file.basename, ".csv"), paste0(output.file.basename, ".xlsx"), sep=" "))

  if ( !missing(ensembl2genes)) {
    output.file.basename.extended <- paste0(output.path, "/", output.name, "_extended")
    ## add real gene names and biotypes to the csv files
    system(paste("./improve_deseq_table.rb", paste0(output.file.basename.extended, ".csv" ), paste0(output.file.basename, ".csv"), ensembl2genes, sep=" "), wait=TRUE)
    system(paste("./csv_to_excel.py", paste0(output.file.basename.extended, ".csv" ), paste0(output.path, "/", output.name, "_extended", ".xlsx"), sep=" "))
  }
}

plot.sample2sample <- function(out, col.labels, trsf_data, trsf_type, colors) {
  ## get sample-to-sample distances
  sampleDists <- dist(t(assay(trsf_data)))
  sampleDistMatrix <- as.matrix(sampleDists)
  ## add names
  rownames(sampleDistMatrix) <- with(colData(trsf_data), col.labels)
  colnames(sampleDistMatrix) <- with(colData(trsf_data), col.labels)


  pdf(paste(out, paste0("sample2sample_", trsf_type, ".pdf"), sep="/"))
  pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, color = colors)
  dev.off()
  svg(paste(out, paste0("sample2sample_", trsf_type, ".svg"), sep="/"))
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

reportingTools.html <- function(out, dds, deseq2.result, pvalueCutoff, condition1,condition2, annotation_genes) {
  # Exporting results to HTML and CSV
  if (pvalueCutoff == 1.1){
    shortName <- 'RNAseq_analysis_with_DESeq2_full'
    title <- paste0('RNA-seq analysis of differential expression using DESeq2, no P value cutoff')
  } else {
    shortName <- paste0('RNAseq_analysis_with_DESeq2_p', pvalueCutoff)
    title <- paste0('RNA-seq analysis of differential expression using DESeq2, P value cutoff ', pvalueCutoff)
  }
  des2Report <- HTMLReport(shortName = shortName, title = title, basePath = out, reportDirectory = "reports/")
  publish(dds, des2Report, pvalueCutoff=pvalueCutoff, annotation.db=NULL, factor = colData(dds)$condition, reportDir=out, n = length(row.names(deseq2.result)), contrast = c("condition",condition1,condition2), make.plots = TRUE)
  finish(des2Report)
  system(paste('./refactor_reportingtools_table.rb', paste0(out, '/reports/', shortName,'.html'), annotation_genes, 'add_plots', sep=" "))
}


##################### HEATMAPS TODO
plot.heat.countmatrix <- function(out, dds, col.labels, count, trsf_data, trsf_type) {
  # Plot a heat map of the count matrix, top basemeans
  select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:count]
  selected.ensembl.ids <- row.names(counts(dds,normalized=TRUE)[select,])
  
  ## read in gene and sample names for replacement
  row_names = c()
  for (gene in selected.ensembl.ids) {
    index = which(ensembl.ids == gene)
    gene_name <- toString(gene.ids[index])
    biotype <- toString(biotype.ids[index])
    row_names <- c(row_names, paste(gene_name, biotype, sep=", "))
  }
  
  ### RAW
  #svg(paste(out,"heatmaps/heatmap_count_matrix_raw.svg",sep=""))
  #heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol, Rowv = TRUE, Colv = TRUE, scale="none", dendrogram="both", trace="none", margin=c(10,6), labCol=col.labels)
  #dev.off()
  
  ### LOG
  #svg(paste(out,"heatmaps/heatmap_count_matrix_log.svg",sep=""))
  #heatmap.2(assay(vsd)[select,], col = hmcol, Rowv = TRUE, Colv = TRUE, scale="none", dendrogram="both", trace="none", margin=c(10, 6), labCol=col.labels)
  #dev.off()
  
  ### LOG STABILIZED
  file <- paste(out,"heatmaps/heatmap_count_matrix_stabilized_",trsf_type,".pdf",sep="")
  pheatmap(assay(trsf_data)[select,], cluster_cols = FALSE, cluster_rows = TRUE, 
           labels_row = row_names, labels_col = col.labels, scale = "row", border_color = NA, 
           height = 12, width = 8, file = file)
  
  # ##### ROW SCALED HEATMAPS OF THE COUNT MATRIX
  # svg(paste(out,"heatmap_count_matrix_raw_rowscaled.svg",sep=""))
  # heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="row", dendrogram="none", trace="none", margin=c(10,6), labCol=col.labels)
  # dev.off()
  # svg(paste(out,"heatmap_count_matrix_log_rowscaled.svg",sep=""))
  # heatmap.2(assay(rld)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="row", dendrogram="none", trace="none", margin=c(10, 6), labCol=col.labels)
  # dev.off()
  # svg(paste(out,"heatmap_count_matrix_stabilized_rowscaled.svg",sep=""))
  # heatmap.2(assay(vsd)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="row", dendrogram="none", trace="none", margin=c(10, 6), labCol=col.labels)
  # dev.off()
}

plot.heat.fc <- function(out, deseq2.res, resFold, dds, col.labels, count, trsf_data, trsf_type) {
  # heat map of log2 foldchanges
  gene_names <- row.names(deseq2.res) #all gene names
  fc_gene_names <- row.names(resFold)[1:count] # for the pattern search for the '30' gene IDs with highest foldchange
  
  fc_gene_names_refac <- numeric()
  for(i in fc_gene_names){
    x <- (paste("^",i,"$", sep=""))
    print(x)
    fc_gene_names_refac <- c(fc_gene_names_refac, x)
  }
  print(fc_gene_names_refac)
  
  pattern <- paste(fc_gene_names_refac, collapse = '|')
  
  isFOLD = grepl( pattern,  gene_names)
  
  select <- order(isFOLD, rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:count]
  
  ## rename ensmbl ids with gene names and biotype
  selected.ensembl.ids <- row.names(counts(dds,normalized=TRUE)[select,])
  
  row_names = c()
  for (gene in selected.ensembl.ids) {
    index = which(ensembl.ids == gene)
    gene_name <- toString(gene.ids[index])
    biotype <- toString(biotype.ids[index])
    row_names <- c(row_names, paste(gene_name, biotype, sep=", "))
  }
  #svg(paste(out,"/heatmaps/heatmap_foldchange_stabilized.svg",sep=""))
  #heatmap.2(assay(vsd)[select,], col = hmcol, Rowv = TRUE, Colv = TRUE, scale="none", 
  #          dendrogram="both", trace="none", margin=c(10, 6), labCol=col.labels, labRow=row_names)
  #dev.off()
  
  file <- paste(out,"heatmaps/heatmap_foldchange_stabilized_",trsf_type,".pdf",sep="")
  pheatmap(assay(trsf_data)[select,], cluster_cols = FALSE, cluster_rows = TRUE, 
           labels_row = row_names, labels_col = col.labels, scale = "row", border_color = NA, 
           height = 12, width = 8, file = file)
  }

##################### PCA TODO
# plot.pca.highest.variance <- function(out, vsd, Pvars, ntops, comparison) {
#   ###############
#   ## Since PCA can be slightly problematic with high dimensional data,
#   ## we first select only the 500 genes showing the highest
#   ## variance.
  
#   point_size = 3
#   point_stroke = 1
#   shape_default = 21
  
#   for (ntop in ntops) {
    
#     select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, length(Pvars)))]
    
#     PCA <- prcomp(t(assay(vsd)[select, ]), scale = F)
#     percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
    
#     dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
#                         PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
#                         sampleNO = colData(vsd)$type,
#                         condition = colData(vsd)$condition)
    
#     rownames(dataGG) = dataGG$sampleNO
    
#     cond1 <- strsplit(strsplit(comparison,":")[[1]][1],'_')[[1]][1] 
#     cond2 <- strsplit(strsplit(comparison,":")[[1]][2],'_')[[1]][1] 
#     time1 <- strsplit(strsplit(comparison,":")[[1]][1],'_')[[1]][2] 
#     time2 <- strsplit(strsplit(comparison,":")[[1]][2],'_')[[1]][2]
#     cond1 <- paste(cond1,time1,sep="_")
#     cond2 <- paste(cond2,time2,sep="_")
#     dataGG$condition <- c(cond1, cond1, cond1, cond2, cond2, cond2)
#     #dataGG$timepoint <- c(time1, time1, time1, time2, time2, time2)
#     dataGG$replicate <- c('N1','N2','N3','N1','N2','N3')
    
# #      ggplot(dataGG, aes(PC1, PC2, colour=condition, fill=timepoint, shape=replicate)) +
#     ggplot(dataGG, aes(PC1, PC2, colour=condition, shape=replicate)) +
#       geom_point(size=point_size, stroke=point_stroke) +
#         #scale_fill_manual(values = my_fillings, breaks = my_fillings_order) + 
#         #scale_shape_manual(values = my_shapes, breaks = my_shapes_order) + 
#         #scale_colour_manual(values = my_colours, breaks = my_infections_order) + 
#         xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#         ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#         ggtitle(paste("PC1 vs PC2: top", ntop, "variable genes")) +
#         guides(
#           shape = guide_legend(order = 2),
#           colour = guide_legend(order = 1, override.aes = list(shape=shape_default))
#         ) +
#         ggsave(paste(out,"statistics/pca_top",ntop,".svg",sep="")) + 
#         ggsave(paste(out,"statistics/pca_top",ntop,".pdf",sep=""))
#   }
# }

# plot.pca <- function(out, col.labels, trsf_data, trsf_type) {
#   # Plot certain Principal Component Analyses of 500 (default) most variable genes
  
#   data <- plotPCA(trsf_data, intgroup=c("condition", "type"), returnData=TRUE) 
#   percentVar <- round(100 * attr(data, "percentVar"))
   
#   ggplot(data, aes(PC1, PC2, color=condition)) +
#     geom_point(size=3) +
#     xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#     ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#     ggtitle(paste("PC1 vs PC2: 500 genes")) + # 500 default
#     ggsave(paste(out,"statistics/pca_simple_",trsf_type,".svg",sep=""))
#     ggsave(paste(out,"statistics/pca_simple_",trsf_type,".svg",sep=""))
  
#   ggplot(data, aes(PC1, PC2, color=condition, shape=col.labels)) +
#     scale_shape_manual(values=1:length(col.labels)) +
#     geom_point(size=3) +
#     xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#     ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#     coord_fixed() + 
#     theme(legend.box = "horizontal") +
#     ggtitle(paste("PC1 vs PC2: 500 genes"))  +
#     ggsave(paste(out,"statistics/pca_ggsave_bigger_fixed_",trsf_type,".svg",sep=""), width=10, height=10)
# }

plot.pca <- function(out, col.labels, trsf_data, trsf_type, ntop) {
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
    ggsave(paste(out, paste0("PCA_simple_", trsf_type, "_top", ntop, ".pdf"), sep="/"))
    ggsave(paste(out, paste0("PCA_simple_", trsf_type, "_top", ntop, ".svg"), sep="/"))

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

# piano <- function(out, resBaseMean, resFold, ensembl) {
#   piano.out <- paste(out,'/piano',sep='')
#   dir.create(piano.out, showWarnings = FALSE)

#   resFold <-resBaseMean[rev(order(abs(resBaseMean$log2FoldChange))),]
#   resSig <- resFold[resFold$padj <= 0.05,]
#   length(rownames(resSig))

#   mapGO <- getBM(attributes = c("ensembl_gene_id","go_id"), # name_1006
#                filters = "ensembl_gene_id",
#                values = rownames(resSig),###resOrdered
#                mart = ensembl)

#   mapGO <- mapGO[mapGO[,2]!="",]
#   write.csv(mapGO, file = paste(piano.out,"/ENSG_GOterm.csv",sep=''), quote = FALSE, row.names = FALSE)

#   head(mapGO)

#   ### filter for go terms that are in biological processes and 20 <= # of genes in go terms < 300
#   #load("/mnt/dessertlocal/mono_pmn_hg_fungi_vit_a_d_deseq/globalData_go_human_resource.RData") # loads gene2descrtiption and gene2go
#   #head(mapGO_new)
#   #mapGO_new <- merge(mapGO, gene2go, by.x = 'go_id', by.y = 'go_term')
#   #mapGO_new$go_id <- NULL
#   #mapGO_new$gene <- NULL
#   mapGO_new <- mapGO

#   myGsc <- loadGSC(mapGO_new)
#   myTval <- resSig$stat
#   names(myTval) <- rownames(resSig)

#   myPval <- resSig$padj
#   names(myPval) <- rownames(resSig)

#   myFC <- resSig$log2FoldChange
#   names(myFC) <- rownames(resSig)

#   perm <- 10000
#   cpus <- 5
#   gene.set.min <- 20
#   gene.set.max <- 100 # 9999999999999
#   gsaRes1 <- runGSA(myTval,geneSetStat="mean", gsc=myGsc, nPerm=perm, 
#                   gsSizeLim=c(gene.set.min,gene.set.max), adjMethod="fdr", ncpus=cpus)
#   gsaRes2 <- runGSA(myTval,geneSetStat="wilcoxon", gsc=myGsc, nPerm=perm, 
#                   gsSizeLim=c(gene.set.min,gene.set.max), adjMethod="fdr", ncpus=cpus)

#   gsaRes3 <- runGSA(myTval,geneSetStat="median",gsc=myGsc,
#                   nPerm=perm,gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
#   gsaRes4 <- runGSA(myTval,geneSetStat="sum",gsc=myGsc,
#                   nPerm=perm,gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
#   #gsaRes5 <- runGSA(myTval,geneSetStat="maxmean",gsc=myGsc,
#   #                nPerm=perm,gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
#   gsaRes6 <- runGSA(myPval,myFC,geneSetStat="fisher",gsc=myGsc,
#                   nPerm=perm,gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
#   gsaRes7 <- runGSA(myPval,myFC,geneSetStat="stouffer",gsc=myGsc,
#                   nPerm=perm,gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
#   gsaRes8 <- runGSA(myPval,myFC,geneSetStat="tailStrength",gsc=myGsc,
#                   nPerm=perm,gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)

#   resList <- list(gsaRes1,gsaRes2,gsaRes3,gsaRes4,gsaRes6,gsaRes7,gsaRes8)
#   #resList <- list(gsaRes1,gsaRes2)
#   names(resList) <- c("mean", "wilcoxon","median","sum","fisher","stouffer","tailStrength")
#   #names(resList) <- c("mean", "wilcoxon")
#   old.par <- par(mar = c(0, 0, 0, 0))
#   par(old.par)

#   if (length(rownames(resSig)) > 150) {
#     pdf(paste(piano.out,"/consensus_heatmap.pdf",sep=""), width = 10, height = 10)
#     ch <- consensusHeatmap(resList,cutoff=10,method="mean",colorkey=FALSE,cellnote="consensusScore",ncharLabel = 120) ## medianPvalue or consensusScore or nGenes
#     dev.off()
#     svg(paste(piano.out,"/consensus_heatmap.svg",sep=""), width = 10, height = 10)
#     ch <- consensusHeatmap(resList,cutoff=10,method="mean",colorkey=FALSE,cellnote="consensusScore",ncharLabel = 120) ## medianPvalue or consensusScore
#     dev.off()
#   }
  
#   downregulated_paths <- ch$pMat[,1][ch$pMat[,1] < 0.05]
#   upregulated_paths <- ch$pMat[,5][ch$pMat[,5] < 0.05]

#   write.table(downregulated_paths, file = paste(piano.out,"/paths_sigdown.csv",sep=''), sep = ";", quote = F, col.names = F, row.names = T)
#   write.table(upregulated_paths, file = paste(piano.out,"/paths_sigup.csv",sep=''), sep = ";", quote = F, col.names = F, row.names = T)
#   system(paste('ruby /mnt/fass2/projects/mh_myotis_rnaseq_weber/scripts/go.rb ',paste(piano.out,"/paths_sigdown.csv",sep=''),sep=""))
#   system(paste('ruby /mnt/fass2/projects/mh_myotis_rnaseq_weber/scripts/go.rb ',paste(piano.out,"/paths_sigup.csv",sep=''),sep=""))
  
#   #new_ch <- data.frame(up = ch$pMat[,1], dn = ch$pMat[,5])
#   #sig_path <- apply(new_ch, 1, min)
#   #new_ch <- new_ch[sig_path < 0.05,]
#   #new_ch_lod <- -log10(new_ch)
#   #new_ch_lod$up[new_ch_lod$up == Inf] <- 4.5
#   #new_ch_lod$dn[new_ch_lod$dn == Inf] <- 4.5

#   #library("pheatmap")
#   #pheatmap(new_ch_lod, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)
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
regionReport_config  <- eval( parse(text=args[10]) )[1]
regionReport_config <- normalizePath(regionReport_config) # regionReport needs the absolute path
#gene.files <- eval( parse(text=args[11]) ) # c("/this/is/file1","/this/is/file2",...) BEST IF THIS DOES NOT HAVE A FILE ENDING LIKE .csv, .txt, ... because used for header and plot titles
#go.terms <- c()
#go.terms <- eval( parse(text=args[12]) ) # c("GO:004563","GO:0011231",...)

#####################
## Read in ensembl ids, gene names and biotypes from a tab seperated table
gene_file <- read.table(ensembl2genes, header=FALSE, sep="\t")
ensembl.ids <- gene_file$V1
gene.ids <- gene_file$V2
biotype.ids <- gene_file$V3

#####################
## Build project structure
out <- paste(project_dir,'/',sep='') # deseq2 dir is created by nextflow in the results dir ()
dir.create(file.path(out), showWarnings = FALSE)
dir.create(file.path(out, 'statistics'), showWarnings = FALSE)
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
df <- data.frame(samples = samples, columns = col.labels, conditions = conditions)
write.table.to.file(df, paste0(out, "/data/input"), "input")

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
write.table.to.file(as.data.frame(norm.counts), paste0(out, "/data/counts"), "normalized_counts")
write.table.to.file(as.data.frame(dds$sizeFactor), paste0(out, "/data/counts"), "sizeFactors")

#####################
## transform counts
rld <- rlog(dds, blind=FALSE)
vsd <- vst(dds, blind=FALSE)

### stuff
# rlogMat <- assay(rld)
# vstMat <- assay(vsd)
#par(mfrow=c(1,3))
#notAllZero <- (rowSums(counts(dds))>0)
## write out the full transformed table, we want to load them later for pathway heatmap in additional script
###end stuff

#####################
## collect transformed counts for easy iterating
transformed.counts = vector(mode="list", length=2)
names(transformed.counts) = c("vsd", "rld")
transformed.counts[[1]] <- vsd; transformed.counts[[2]] <- rld

#####################
## write transformed counts 
for (i in 1:length(transformed.counts)) {
  write.table.to.file(as.data.frame(assay(transformed.counts[[i]])), paste0(out, "/data/counts"), paste0("transformed_counts_", names(transformed.counts)[[i]]))
}

##########################################
## Visualisation
##########################################

## BIOMART OBJECT
# the below code works actually. But we need to know the reference species. But for example not for ecoli because I think that's a different ensembl db
# we can use this for hsapiens, mmusculus, mlucifugus, ... 
#ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host = "apr2020.archive.ensembl.org")
#mart <- useDataset("mmusculus_gene_ensembl", ensembl)

##########################################
## Visualisation
##########################################

#####################
## PCA
for (i in 1:length(transformed.counts)) { 
  for (ntop in c(500, 100, 50)){
    output.pca.dir <- paste(out, "plots/PCA/", sep="/")

    plot.pca(output.pca.dir, col.labels, transformed.counts[[i]], names(transformed.counts)[[i]], ntop)

    if (length(patients) > 0) {
      pcaData <- plotPCA(transformed.counts[[i]], intgroup=c("condition", "type", "patients"), ntop=ntop, returnData=TRUE)
      percentVar <- round(100 * attr(pcaData, "percentVar"))
      ggplot(pcaData, aes(PC1, PC2, color=patients:type, shape=condition)) + 
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed() +
        ggtitle(paste("PC1 vs PC2: top ", ntop, " variable genes")) +
        ggsave(paste(out,"statistics/pca_detailed_top_",names(transformed.counts)[[i]],"_",ntop,".svg",sep="")) +
        ggsave(paste(out,"statistics/pca_detailed_top_",names(transformed.counts)[[i]],"_",ntop,".pdf",sep=""))
    } else{
      pcaData <- plotPCA(transformed.counts[[i]], intgroup=c("condition", "type"), ntop=ntop, returnData=TRUE)
      percentVar <- round(100 * attr(pcaData, "percentVar"))
      ggplot(pcaData, aes(PC1, PC2, color=type, shape=condition)) + 
        geom_point(size=3) + 
        xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed() +
        ggtitle(paste("PC1 vs PC2: top ", ntop, " variable genes")) +
        ggsave(paste(out,"statistics/pca_detailed_top_",names(transformed.counts)[[i]],"_",ntop,".svg",sep="")) +
        ggsave(paste(out,"statistics/pca_detailed_top_",names(transformed.counts)[[i]],"_",ntop,".pdf",sep=""))
    }
  }
}

#####################
## Heatmap

## TODO

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

  transformed.counts.sub = vector(mode="list", length=2)
  names(transformed.counts.sub) = c("vsd", "rld")
  transformed.counts.sub[[1]] <- vsd; transformed.counts.sub[[2]] <- rld

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
  resOrdered <<- deseq2.res[order(deseq2.res$padj),]

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
  df.sub <- data.frame(samples = samples.sub, columns = col.labels.sub, conditions = conditions.sub)
  write.table.to.file(df.sub,  paste0(out.sub, "/input"), "input")

  #####################
  ## DESeq2 results
  out.sub.output.dir <- paste0(out.sub, "/results/")

  write.table.to.file(as.data.frame(resOrdered), out.sub.output.dir, paste(name, "full", sep="_"), ensembl2genes)
  write.table.to.file(as.data.frame(resFold), out.sub.output.dir, paste(name, "filtered_NA", sep="_"), ensembl2genes)
  write.table.to.file(as.data.frame(resFold05), out.sub.output.dir, paste(name, "filtered_padj_0.05", sep="_"), ensembl2genes)
  write.table.to.file(as.data.frame(resFold01), out.sub.output.dir, paste(name, "filtered_padj_0.01", sep="_"), ensembl2genes)

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
  volcano = EnhancedVolcano(deseq2.res, lab = rownames(deseq2.res), x = 'log2FoldChange', y = 'padj', 
    legendLabels = c('NS', expression(Log[2]~FC), "adj. p-value", expression(adj.~p-value~and~log[2]~FC)))
  volcano + 
    ggsave(paste(out.sub,"/plots/volcano/volcano.svg", sep='/')) +
    ggsave(paste(out.sub,"/plots/volcano/volcano.pdf", sep='/'))

  #####################
  ## MA plots
  plot.ma(paste(out.sub, "/plots/MA/", sep="/"), deseq2.res, metadata(deseq2.res)$alpha)
  plot.ma(paste(out.sub, "/plots/MA/", sep="/"), deseq2.res, metadata(deseq2.res)$alpha / 2)

  #####################
  ## Heatmap of sample2sample distances
  for (i in 1:length(transformed.counts.sub)) {
    plot.sample2sample(paste(out.sub, "/plots/sample2sample/", sep="/"), col.labels.sub, 
      transformed.counts.sub[[i]], names(transformed.counts.sub)[[i]], colorRampPalette( rev(brewer.pal(9, "Blues")) )(255))
  }

  #####################
  ## Heatmap: count and fc
  for (i in 1:length(transformed.counts.sub)) {
    ## HEATMAPs
    #TODO PHEATMAP REBUILD!!!
    #TODO AND BUILD HEATMAP BASED ON GENE LIST
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    plot.heat.countmatrix(out.sub, dds.sub, col.labels.sub, 50, transformed.counts.sub[[i]], names(transformed.counts.sub)[[i]])
    plot.heat.fc(out.sub, deseq2.res, resFold, dds.sub, col.labels.sub, 50, transformed.counts.sub[[i]], names(transformed.counts.sub)[[i]])
  }
  #########################
  # ## Heat map of the count matrix
  # #########################

  # select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
  # selected.ensembl.ids <- row.names(counts(dds,normalized=TRUE)[select,])

  # ## read in gene and sample names for replacement
  # row_names = c()
  # for (gene in selected.ensembl.ids) {
  #   index = which(ensembl.ids == gene)
  #   gene_name <- toString(gene.ids[index])
  #   biotype <- toString(biotype.ids[index])
  #   row_names <- c(row_names, paste(gene_name, biotype, sep=", "))
  # }

  # ### LOG STABILIZED
  # for (i in 1:length(transformed.counts)) {
  #   file <- paste(out,"heatmaps/heatmap_count_matrix_row-scaled_",names(transformed.counts)[[i]],".pdf",sep="")
  #   pheatmap(assay(transformed.counts[[i]])[select,], cluster_cols = FALSE, cluster_rows = TRUE,
  #           labels_row = row_names, labels_col = col.labels, scale = "row", border_color = NA,
  #           height = 12, width = 8, file = file)
  # }

  #####################
  ## PCA

  for (i in 1:length(transformed.counts.sub)) {
    for (ntop in c(500, 100, 50)) {
      plot.pca(paste(out.sub, "/plots/PCA/", sep="/"), col.labels.sub, transformed.counts.sub[[i]], names(transformed.counts.sub)[[i]], ntop)
    }
    # below would generate a nice PCA, but we need to generlize this first. And maybe re-think the way we are providing information about replicates, timepoints, patients, ...
    # Pvars.sub <- rowVars(assay(vsd.sub))
    #plot.pca.highest.variance(out.sub, transformed.counts.sub[[i]], Pvars.sub, ntops, comparison)
  }

  ##########################################
  ## Further analysis
  ##########################################
  #data.set <- rownames(deseq2.res)
  #results.gene <- getBM(attributes = c("ensembl_gene_id","external_gene_name","go_id","name_1006"),  filters="ensembl_gene_id",values = data.set, mart=mart)

  ## piano
  #piano(out.sub, resBaseMean, resFold, ensembl)

  # nice extension for later to color genes that belong to specific GO terms
  # if (length(go.terms) > 0) {
  #   plot.ma.go(out.sub, deseq2.res, ma.size, results.gene, go.terms, transformed.counts.sub[[i]], names(transformed.counts.sub)[[i]])
  # }

  ##########################################
  ## Reports
  ##########################################

  #####################
  ## regionReport report

  ## set output
  report.project.name <- paste(l1, "vs", l2, sep=" ")
  report.dir <- paste(out.sub, "reports", sep="/")
  report.output <- paste0('DESeq2_results_exploration')

  ## create html
  report_html <- DESeq2Report(dds, project = report.project.name,
    intgroup = c('condition', 'type'), res = deseq2.res, template = regionReport_config,
    outdir = report.dir, output = report.output, theme = theme_bw())

  ## and also pfd
  report_pdf <- DESeq2Report(dds, project = report.project.name,
    intgroup = c('condition', 'type'), res = deseq2.res, template = regionReport_config,
    outdir = report.dir, output = report.output, theme = theme_bw(),
    output_format = 'pdf_document', device = 'pdf')

  #####################
  ## ReportingTools
  reportingTools.html(out.sub, dds, deseq2.res, 1.1, l1, l2, annotation_genes)
  if (length(rownames(resFold05)) > 0) {
    reportingTools.html(out.sub, dds, deseq2.res, 0.05, l1, l2, annotation_genes)
  }
  if (length(rownames(resFold01)) > 0) {
    reportingTools.html(out.sub, dds, deseq2.res, 0.01, l1, l2, annotation_genes)
  }

}
#####################################################################################
## END PAIRWISE COMPARISONS
#####################################################################################

