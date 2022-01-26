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
  system(paste('./refactor_reportingtools_table.rb', paste0(out.dir, '/reports/', shortName,'.html'), annotation_genes, 'add_plots', pvalueCutoff, sep=" "))
}


  ## ReportingTools
  reportingTools.html(out.sub, dds, deseq2.res, 1.1, l2, l1, annotation_genes)
  if (length(rownames(resFold05)) > 0) {
    reportingTools.html(out.sub, dds, deseq2.res, 0.05, l2, l1, annotation_genes, make.plots=FALSE)
  }
  if (length(rownames(resFold01)) > 0) {
    reportingTools.html(out.sub, dds, deseq2.res, 0.01, l2, l1, annotation_genes, make.plots=FALSE)
  }