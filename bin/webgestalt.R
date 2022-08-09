library("WebGestaltR")

# read arguments and evaluate
args <- commandArgs(TRUE)

project_dir <- eval( parse(text=args[1]) )[1] 
gene_file <- eval( parse(text=args[2]) )
species <- eval( parse(text=args[3]) )
id_type <- eval( parse(text=args[4]) )

resFold05 <- read.csv(file = gene_file, row.names = 1)


#####################
## Webgestalt
#####################
if ( species == 'hsa' ){
organism <- "hsapiens"
} else if (species == 'mmu') {
organism <- "mmusculus"
} else {
organism <- NA
}
if (! is.na(organism)) {
    dir.create(file.path('WebGestalt'), showWarnings = FALSE, recursive = TRUE)
    interestGene <- as.data.frame(resFold05)[, 'log2FoldChange', drop=FALSE]
    interestGene$id <- rownames(interestGene)
    rownames(interestGene) <- NULL
    colnames(interestGene) <- NULL
    interestGene <- interestGene[c(2,1)]
    webgestalt.out.dir <- paste("WebGestalt", sep='/')
    if (any(grepl(id_type, listIdType(), fixed=TRUE))) {
        try.webgestalt <- try(
        for (enrDB in c("geneontology_Biological_Process_noRedundant", "pathway_KEGG")){
            enrichResult <- WebGestaltR(enrichMethod="GSEA", organism=organism, enrichDatabase=enrDB, interestGene=interestGene, interestGeneType=id_type, collapseMethod="mean", minNum=10, maxNum=500, fdrMethod="BH", sigMethod="fdr", fdrThr=0.01, topThr=10, perNum=1000, isOutput=TRUE, outputDirectory=webgestalt.out.dir, projectName="pairwise_comparison")
        }
        )
        if (class(try.webgestalt) == "try-error") {
        print('SKIPPING: WebGestaltR. The number of annotated IDs for all functional categories are not from 10 to 500 for the GSEA enrichment method.')
        }
    } else {
        print(paste('SKIPPING: WebGestaltR. Feature ID', id_type, 'not supported.'))
    }
}