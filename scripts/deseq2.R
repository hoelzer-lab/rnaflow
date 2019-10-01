library("DESeq2")
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("ReportingTools")
library("pheatmap")

###############################################################################################
## FUNCTIONS
###############################################################################################
plot.sample2sample <- function(out, dds, rld, col.labels) {
  ## Heat map of the sample-to-sample distances
  
  distsRL <- dist(t(assay(rld)))
  
  mat <- as.matrix(distsRL)
  rownames(mat) <- colnames(mat) <- with(colData(dds), col.labels)
  hc <- hclust(distsRL)
  
  pdf(paste(out,"/heatmaps/heatmap_sample2sample.pdf",sep=""))
  heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col = rev(hmcol), margin=c(13, 13))
  dev.off()
}

report.html <- function(out, dds, deseq2.res, l1, l2, use.contrasts) {
  ## Exporting results to HTML and CSV
  db <- NULL

  des2Report05 <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2_p05', title = 'RNA-seq analysis of differential expression using DESeq2, pvalue cutoff 0.05', basePath = out, reportDirectory = "html/")
  if (use.contrasts) {
    publish(dds, des2Report05, pvalueCutoff=0.05, annotation.db=db, factor = colData(dds)$condition, reportDir=out, n = length(row.names(deseq2.res)), contrast = c("condition",l1,l2))
  } else {
    publish(dds, des2Report05, pvalueCutoff=0.05, annotation.db=db, factor = colData(dds)$condition, reportDir=out, n = length(row.names(deseq2.res)))
  }
  finish(des2Report05)
  system(paste('/mnt/mahlzeitlocal/projects/myotis_rnaseq_weber/scripts/refactor_reportingtools_table.rb ', out, '/html/', 'RNAseq_analysis_with_DESeq2_p05.html', sep=''))
}


plot.heat.countmatrix <- function(out, dds, vsd, col.labels, count) {
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
  #heatmap.2(assay(rld)[select,], col = hmcol, Rowv = TRUE, Colv = TRUE, scale="none", dendrogram="both", trace="none", margin=c(10, 6), labCol=col.labels)
  #dev.off()
  
  ### LOG STABILIZED
  file <- paste(out,"heatmaps/heatmap_count_matrix_stabilized.pdf",sep="")
  pheatmap(assay(vsd)[select,], cluster_cols = FALSE, cluster_rows = TRUE, 
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

plot.heat.fc <- function(out, deseq2.res, resFold, dds, vsd, col.labels, count) {
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
  
  file <- paste(out,"heatmaps/heatmap_foldchange_stabilized.pdf",sep="")
  pheatmap(assay(vsd)[select,], cluster_cols = FALSE, cluster_rows = TRUE, 
           labels_row = row_names, labels_col = col.labels, scale = "row", border_color = NA, 
           height = 12, width = 8, file = file)
  }


plot.pca.highest.variance <- function(out, rld, Pvars, ntops, comparison) {
  ###############
  ## Since PCA can be slightly problematic with high dimensional data,
  ## we first select only the 500 genes showing the highest
  ## variance.
  
  point_size = 3
  point_stroke = 1
  shape_default = 21
  
  for (ntop in ntops) {
    
    select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, length(Pvars)))]
    
    PCA <- prcomp(t(assay(rld)[select, ]), scale = F)
    percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
    
    dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
                        PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
                        sampleNO = colData(rld)$type,
                        condition = colData(rld)$condition)
    
    rownames(dataGG) = dataGG$sampleNO
    
    cond1 <- strsplit(strsplit(comparison,":")[[1]][1],'_')[[1]][1] 
    cond2 <- strsplit(strsplit(comparison,":")[[1]][2],'_')[[1]][1] 
    time1 <- strsplit(strsplit(comparison,":")[[1]][1],'_')[[1]][2] 
    time2 <- strsplit(strsplit(comparison,":")[[1]][2],'_')[[1]][2]
    cond1 <- paste(cond1,time1,sep="_")
    cond2 <- paste(cond2,time2,sep="_")
    dataGG$condition <- c(cond1, cond1, cond1, cond2, cond2, cond2)
    dataGG$timepoint <- c(time1, time1, time1, time2, time2, time2)
    dataGG$replicate <- c('N1','N2','N3','N1','N2','N3')
    
#      ggplot(dataGG, aes(PC1, PC2, colour=condition, fill=timepoint, shape=replicate)) +
    ggplot(dataGG, aes(PC1, PC2, colour=condition, shape=replicate)) +
      geom_point(size=point_size, stroke=point_stroke) +
        #scale_fill_manual(values = my_fillings, breaks = my_fillings_order) + 
        #scale_shape_manual(values = my_shapes, breaks = my_shapes_order) + 
        #scale_colour_manual(values = my_colours, breaks = my_infections_order) + 
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        ggtitle(paste("PC1 vs PC2: top", ntop, "variable genes")) +
        guides(
          shape = guide_legend(order = 2),
          colour = guide_legend(order = 1, override.aes = list(shape=shape_default))
        ) +
        ggsave(paste(out,"statistics/pca_top",ntop,".svg",sep="")) + 
        ggsave(paste(out,"statistics/pca_top",ntop,".pdf",sep=""))
  }
}


plot.pca <- function(out, rld, col.labels, patients) {
  # Plot certain Principal Component Analyses.
  
  head(colData(rld))
  
  ## old plot with less information
  #pdf(paste(out,"statistics/pca_simple.pdf",sep=""))
  #plotPCA(rld, intgroup=c("condition", "type")) #"sizeFactor" worked somehow....
  #dev.off()
  
  #pdf(paste(out,"statistics/pca.pdf",sep=""))
  data <- plotPCA(rld, intgroup=c("condition", "type"), returnData=TRUE) 
  percentVar <- round(100 * attr(data, "percentVar"))
  
  #ggplot(data, aes(PC1, PC2, color=condition, shape=col.labels)) +
  #  scale_shape_manual(values=1:length(col.labels)) +
  #  geom_point(size=3) +
  #  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  #  ylab(paste0("PC2: ",percentVar[2],"% variance"))
  #dev.off()
  
  ggplot(data, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    ggtitle(paste("PC1 vs PC2: ", length(rownames(rld)), " genes")) +
    ggsave(paste(out,"statistics/pca_simple.svg",sep=""))
  
    ggplot(data, aes(PC1, PC2, color=condition, shape=col.labels)) +
    scale_shape_manual(values=1:length(col.labels)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() + 
    theme(legend.box = "horizontal") +
    ggtitle(paste("PC1 vs PC2: ", length(rownames(rld)), " genes")) +
    ggsave(paste(out,"statistics/pca_ggsave_bigger_fixed.svg",sep=""), width=10, height=10)
  
}


build.project.structure <- function(out) {
  dir.create(file.path(out), showWarnings = FALSE)
  # Build necessary project structure.
  dir.create(file.path(out, '/statistics'), showWarnings = FALSE)
  dir.create(file.path(out, '/heatmaps'), showWarnings = FALSE)
  dir.create(file.path(out, '/html'), showWarnings = FALSE)
  dir.create(file.path(out, '/tmp'), showWarnings = FALSE)
}


plot.ma <- function(out, deseq2.res, ma.size, rld) {
  ##############################
  ## MA plot
  ############################## 
  # These plots show the log2 fold changes from the treatment over the
  # mean of normalized counts, i.e. the average of counts normalized by
  # size factors. The left plot shows the “unshrunken” log2 fold changes,
  # while the right plot, produced by the code above, shows the shrinkage
  # of log2 fold changes resulting from the incorporation of zero-centered
  # normal prior. The shrinkage is greater for the log2 fold change
  # estimates from genes with low counts and high dispersion, as can be
  # seen by the narrowing of spread of leftmost points in the right plot.
  ##################
  
  pdf(paste(out,"statistics/ma.pdf",sep=""))
  plotMA(deseq2.res, main="DESeq2", ylim=ma.size)
  dev.off()
}

plot.ma.go <- function(out, deseq2.res, ma.size, rld, results.gene, go.terms) {
  ## We can also make an MA-plot for the results table in which we raised
  ## the log2 fold change threshold (Figure below). We can label individual
  ## points on the MA-plot as well. Here we use the with R function to plot
  ## a circle and text for a selected row of the results object. Within the
  ## with function, only the baseMean and log2FoldChange values for the
  ## selected rows of res are used.
  ##-----------------------------
  for (go.term.ma in go.terms) {
    #go.term.ma <- "GO:0009615"
    pdf(paste(out,"statistics/ma_", gsub(":", "", go.term.ma), ".pdf",sep=""))
    plotMA(deseq2.res, main=paste("DESeq2, ", go.term.ma, sep=''), ylim=ma.size)
    results.gene.GO.ma <- grep(go.term.ma, results.gene$go_id, fixed=TRUE)  ### e.g. GO:0002376, immune system process in mice
    rld.go.ma <- rownames(assay(rld)[results.gene[results.gene.GO.ma,]$ensembl_gene_id,]) # get the ensembl ids corresponding to this go term
    for (gene in rld.go.ma) {
      index = which(ensembl.ids == gene)
      gene.name <- toString(gene.ids[index])
      with(deseq2.res[gene, ], {
        if (gene %in% rownames(resFold05)) {
          points(baseMean, log2FoldChange, col="dodgerblue", cex=0.8, lwd=2, bg="dodgerblue")
          text(baseMean, log2FoldChange, gene.name, pos=2, col="dodgerblue")
        }
      })
    }
    dev.off()
  }
}

piano <- function(out, resBaseMean, resFold, ensembl) {
  piano.out <- paste(out,'/piano',sep='')
  dir.create(piano.out, showWarnings = FALSE)

  resFold <-resBaseMean[rev(order(abs(resBaseMean$log2FoldChange))),]
  resSig <- resFold[resFold$padj <= 0.05,]
  length(rownames(resSig))

  mapGO <- getBM(attributes = c("ensembl_gene_id","go_id"), # name_1006
               filters = "ensembl_gene_id",
               values = rownames(resSig),###resOrdered
               mart = ensembl)

  mapGO <- mapGO[mapGO[,2]!="",]
  write.csv(mapGO, file = paste(piano.out,"/ENSG_GOterm.csv",sep=''), quote = FALSE, row.names = FALSE)

  head(mapGO)

  ### filter for go terms that are in biological processes and 20 <= # of genes in go terms < 300
  #load("/mnt/dessertlocal/mono_pmn_hg_fungi_vit_a_d_deseq/globalData_go_human_resource.RData") # loads gene2descrtiption and gene2go
  #head(mapGO_new)
  #mapGO_new <- merge(mapGO, gene2go, by.x = 'go_id', by.y = 'go_term')
  #mapGO_new$go_id <- NULL
  #mapGO_new$gene <- NULL
  mapGO_new <- mapGO

  myGsc <- loadGSC(mapGO_new)
  myTval <- resSig$stat
  names(myTval) <- rownames(resSig)

  myPval <- resSig$padj
  names(myPval) <- rownames(resSig)

  myFC <- resSig$log2FoldChange
  names(myFC) <- rownames(resSig)

  perm <- 10000
  cpus <- 5
  gene.set.min <- 20
  gene.set.max <- 100 # 9999999999999
  gsaRes1 <- runGSA(myTval,geneSetStat="mean", gsc=myGsc, nPerm=perm, 
                  gsSizeLim=c(gene.set.min,gene.set.max), adjMethod="fdr", ncpus=cpus)
  gsaRes2 <- runGSA(myTval,geneSetStat="wilcoxon", gsc=myGsc, nPerm=perm, 
                  gsSizeLim=c(gene.set.min,gene.set.max), adjMethod="fdr", ncpus=cpus)

  gsaRes3 <- runGSA(myTval,geneSetStat="median",gsc=myGsc,
                  nPerm=perm,gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
  gsaRes4 <- runGSA(myTval,geneSetStat="sum",gsc=myGsc,
                  nPerm=perm,gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
  #gsaRes5 <- runGSA(myTval,geneSetStat="maxmean",gsc=myGsc,
  #                nPerm=perm,gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
  gsaRes6 <- runGSA(myPval,myFC,geneSetStat="fisher",gsc=myGsc,
                  nPerm=perm,gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
  gsaRes7 <- runGSA(myPval,myFC,geneSetStat="stouffer",gsc=myGsc,
                  nPerm=perm,gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)
  gsaRes8 <- runGSA(myPval,myFC,geneSetStat="tailStrength",gsc=myGsc,
                  nPerm=perm,gsSizeLim=c(gene.set.min,gene.set.max), ncpus=cpus)

  resList <- list(gsaRes1,gsaRes2,gsaRes3,gsaRes4,gsaRes6,gsaRes7,gsaRes8)
  #resList <- list(gsaRes1,gsaRes2)
  names(resList) <- c("mean", "wilcoxon","median","sum","fisher","stouffer","tailStrength")
  #names(resList) <- c("mean", "wilcoxon")
  old.par <- par(mar = c(0, 0, 0, 0))
  par(old.par)

  if (length(rownames(resSig)) > 150) {
    pdf(paste(piano.out,"/consensus_heatmap.pdf",sep=""), width = 10, height = 10)
    ch <- consensusHeatmap(resList,cutoff=10,method="mean",colorkey=FALSE,cellnote="consensusScore",ncharLabel = 120) ## medianPvalue or consensusScore or nGenes
    dev.off()
    svg(paste(piano.out,"/consensus_heatmap.svg",sep=""), width = 10, height = 10)
    ch <- consensusHeatmap(resList,cutoff=10,method="mean",colorkey=FALSE,cellnote="consensusScore",ncharLabel = 120) ## medianPvalue or consensusScore
    dev.off()
  }
  
  downregulated_paths <- ch$pMat[,1][ch$pMat[,1] < 0.05]
  upregulated_paths <- ch$pMat[,5][ch$pMat[,5] < 0.05]

  write.table(downregulated_paths, file = paste(piano.out,"/paths_sigdown.csv",sep=''), sep = ";", quote = F, col.names = F, row.names = T)
  write.table(upregulated_paths, file = paste(piano.out,"/paths_sigup.csv",sep=''), sep = ";", quote = F, col.names = F, row.names = T)
  system(paste('ruby /mnt/fass2/projects/mh_myotis_rnaseq_weber/scripts/go.rb ',paste(piano.out,"/paths_sigdown.csv",sep=''),sep=""))
  system(paste('ruby /mnt/fass2/projects/mh_myotis_rnaseq_weber/scripts/go.rb ',paste(piano.out,"/paths_sigup.csv",sep=''),sep=""))
  
  #new_ch <- data.frame(up = ch$pMat[,1], dn = ch$pMat[,5])
  #sig_path <- apply(new_ch, 1, min)
  #new_ch <- new_ch[sig_path < 0.05,]
  #new_ch_lod <- -log10(new_ch)
  #new_ch_lod$up[new_ch_lod$up == Inf] <- 4.5
  #new_ch_lod$dn[new_ch_lod$dn == Inf] <- 4.5

  #library("pheatmap")
  #pheatmap(new_ch_lod, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)
}
#####################################################################################
## END FUNCTIONS
#####################################################################################


#####################################################################################
#####################################################################################
#####################################################################################

### RUN THESE SCRIPT
# R CMD BATCH --no-save --no-restore '--args c("project_dir") c("a","b","c","d","e","f","g") c(1,2,3)' /home/hoelzer/scripts/R/deseq2.R test.out

#R CMD BATCH --no-save --no-restore '--args c("/home/hoelzer/git/nanozoo/wf_gene_expression/results/") c("results/03-Counting/mock_rep1.counts.formated","results/03-Counting/mock_rep2.counts.formated","results/03-Counting/mock_rep3.counts.formated","results/03-Counting/treated_rep1.counts.formated","results/03-Counting/treated_rep2.counts.formated","results/03-Counting/treated_rep3.counts.formated") c("mock","mock","mock","treated","treated","treated") c("mock_rep1","mock_rep2","mock_rep3","treated_rep1","treated_rep2","treated_rep3") c("mock","treated") c("mock:treated") c("data/db/Rattus_norvegicus.Rnor_6.0.91.chr.id2name") c("rno") c()' scripts/deseq2.R

#######################
# This just reads the two arguments passed from the command line
# and assigns them to a vector of characters.
args <- commandArgs(TRUE)
 
# Parse the arguments (in characters) and evaluate them
project_dir <- eval( parse(text=args[1]) )[1] 
samples <- eval( parse(text=args[2]) )
conditions <- eval( parse(text=args[3]) )
col.labels <- eval( parse(text=args[4]) )
levels <- eval( parse(text=args[5]) )
comparisons <- eval( parse(text=args[6]) )
out <- paste(project_dir,'deseq2/',sep='/')
ensembl2genes <- eval( parse(text=args[7]) )[1]
species <- eval( parse(text=args[8]) )[1]

ntops <- c(500)
patients <- eval( parse(text=args[9]) ) # patients is a vector like c("1","1","1","2","2","2"), so if we have samples from the same patient we want to use as replicates, and not all vs all
#gene.files <- eval( parse(text=args[10]) ) # c("/this/is/file1","/this/is/file2",...) BEST IF THIS DOES NOT HAVE A FILE ENDING LIKE .csv, .txt, ... because used for header and plot titles
#go.terms <- eval( parse(text=args[11]) ) # c("GO:004563","GO:0011231",...)

#name <- paste("deseq2_",levels[1],"_",levels[2],sep="")

#####################
## read in ensembl ids and gene names, as well as biotypes
gene_file <- read.table(ensembl2genes, header=FALSE, sep="\t")
ensembl.ids <- gene_file$V1
gene.ids <- gene_file$V2
biotype.ids <- gene_file$V3

#######################
## build project dirs
#######################
dir.create(file.path(out), showWarnings = FALSE)
dir.create(file.path(out, 'statistics'), showWarnings = FALSE)
dir.create(file.path(out, 'heatmaps'), showWarnings = FALSE)
dir.create(file.path(out, 'html'), showWarnings = FALSE)
dir.create(file.path(out, 'tmp'), showWarnings = FALSE)

#######################
## write out the input files
df <- data.frame(samples = samples, columns = col.labels, conditions = conditions)
input.csv <- paste(out,"input.csv",sep="")
write.csv(as.data.frame(df), file=input.csv)

if (length(patients) > 0) {
    sampleTable <- data.frame(sampleName = samples, fileName = samples, condition = conditions, type = col.labels, patients = patients)
    ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = "", design= ~ patients + condition)
} else {
    sampleTable <- data.frame(sampleName = samples, fileName = samples, condition = conditions, type = col.labels)
    ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = "", design= ~ condition)
}

### IMPORTANT STEP< OTHERWIESE DESEQ WILL DO THE COMPARISON IN ALPHABETICAL ORDER!!!!
#ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref=levels[1]) ## this is enough if we only have two conditions, 
#but for more conditions we need factor() to order every level according to input files
ddsHTSeq$condition <- factor(ddsHTSeq$condition, levels=levels)
ddsHTSeq$type <- factor(ddsHTSeq$type, levels=col.labels)
print(ddsHTSeq)
print(ddsHTSeq$condition)

print("DESeq Data Object:")
dds <- DESeq(ddsHTSeq, betaPrior = TRUE)
head(dds)

## write out this table to have the size factors and normalized read
## counts for each gene and sample
norm.counts <- counts(dds, normalized=T)
csv <- paste(out,"normalized_counts.csv",sep="")
write.csv(as.data.frame(norm.counts), file=csv)
csv <- paste(out,"sizeFactors.txt",sep="")
write.csv(as.data.frame(dds$sizeFactor), file=csv)

###################################
## Extracting transformed values
###################################
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay(rld)
vstMat <- assay(vsd)

par(mfrow=c(1,3))
notAllZero <- (rowSums(counts(dds))>0)

## write out the full vsd table, we want to load them later for pathway heatmap in additional script
csv <- paste(out,"/normalized_counts.vsd.csv",sep="")
write.csv(as.data.frame(assay(vsd)), file=csv)

## BIOMART OBJECT
#ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mlucifugus_gene_ensembl", host = "may2017.archive.ensembl.org")
#mart <- useDataset("mlucifugus_gene_ensembl", ensembl)

#simple pca
pdf(paste(out,"statistics/pca_simple.pdf",sep=""))
plotPCA(rld, intgroup=c("condition", "type")) #"sizeFactor" worked somehow....
dev.off()

ntop = 500

## Heat map of the sample-to-sample distances
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
#hmcol <- colorRampPalette(brewer.pal(9, "RdBu"))(100) # this is a red/blue color map

distsRL <- dist(t(assay(rld)))

mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), col.labels)
hc <- hclust(distsRL)

pdf(paste(out,"heatmaps/heatmap_sample2sample.pdf",sep=""))
heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()


#########################
## Heat map of the count matrix
#########################

select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
selected.ensembl.ids <- row.names(counts(dds,normalized=TRUE)[select,])

## read in gene and sample names for replacement
row_names = c()
for (gene in selected.ensembl.ids) {
  index = which(ensembl.ids == gene)
  gene_name <- toString(gene.ids[index])
  biotype <- toString(biotype.ids[index])
  row_names <- c(row_names, paste(gene_name, biotype, sep=", "))
}

### LOG STABILIZED
file <- paste(out,"heatmaps/heatmap_count_matrix_row-scaled.pdf",sep="")
pheatmap(assay(vsd)[select,], cluster_cols = FALSE, cluster_rows = TRUE, 
         labels_row = row_names, labels_col = col.labels, scale = "row", border_color = NA, 
         height = 12, width = 8, file = file)



## REPORT TO HTML
db <- NULL

des2Report.full <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2_full', title = 'RNA-seq analysis of differential expression using DESeq2, pvalue cutoff 1', basePath = out, reportDirectory = "html/")
publish(dds, des2Report.full, pvalueCutoff=1, annotation.db=db, factor = colData(dds)$condition, reportDir=out, n = length(row.names(dds)))
finish(des2Report.full)
#system(paste('/mnt/mahlzeitlocal/projects/myotis_rnaseq_weber/scripts/refactor_reportingtools_table.rb ', out, '/html/', 'RNAseq_analysis_with_DESeq2_full.html', sep=''))
