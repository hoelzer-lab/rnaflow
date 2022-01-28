library("DESeq2")
  
  

# read in args and evaluate
args <- commandArgs(TRUE)

comparison <- eval( parse(text=args[1]) )
  
#####################################################################################
## PERFORM PAIRWISE COMPARISONS
#####################################################################################

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
numerator <- l1
denominator <- l2
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

write.table.to.file(as.data.frame(resOrdered), out.sub.output.dir, paste(name, "full", sep="_"), id2name, col.names=NA)
write.table.to.file(as.data.frame(resFold), out.sub.output.dir, paste(name, "filtered_NA", sep="_"), id2name, col.names=NA)
write.table.to.file(as.data.frame(resFold05), out.sub.output.dir, paste(name, "filtered_padj_0.05", sep="_"), id2name, col.names=NA)
write.table.to.file(as.data.frame(resFold01), out.sub.output.dir, paste(name, "filtered_padj_0.01", sep="_"), id2name, col.names=NA)

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
    plot.heatmap.top_fc(paste(out.sub, "plots/heatmaps/", sep="/"), resFold05, transformed.counts.sub[[i]], names(transformed.counts.sub)[[i]], ntop, 'pcutoff0-05')
    plot.heatmap.top_fc(paste(out.sub, "plots/heatmaps/", sep="/"), resFold01, transformed.counts.sub[[i]], names(transformed.counts.sub)[[i]], ntop, 'pcutoff0-01')
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

#####################################################################################
## END PAIRWISE COMPARISONS
#####################################################################################

