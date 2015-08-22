library("knitr")
library('clusterProfiler')
library(GO.db)
library("org.Dm.eg.db")

## Parameters
dir.main <- "~/dr-chip-rna-seq"
setwd(dir.main)

dir.deg <- file.path(dir.main, "results/rna-seq/DEG")
dir.enrichment <- file.path(dir.deg, "enrichment_analysis")
dir.create(dir.enrichment, showWarnings = FALSE, recursive = TRUE)

cond1 <- "S2"
cond2 <- "WT"
comparison <- paste(sep="", cond1, "_vs_", cond2)

parameters <- list(
  "dir.deg" = dir.deg,
  "dir.enrichment" =  dir.enrichment,
  "condition1" = cond1,
  "condition1" = cond2,
  "deg.threshold" = 0.05,
  "enrich.threshold" = 0.05,
  "groupGo.level" = 3,
  "barplot_catnb" = 36)


## List of tools for differential expression analysis
deg.tools <- c("edgeR", "DESeq2")
for (deg.tool in deg.tools) {
  parameters[[deg.tool]] <- 
    file.path(
      dir.deg, comparison, 
      paste(sep="", cond1, "_VS_", cond2, "_sickle_se_q20_subread_featurecounts_",deg.tool,".tab"))
}  

## Print out the parameters in the Rmd report
param.table <- t(as.data.frame(parameters))
colnames(param.table) <- "value"
# kable(param.table)


## Create a table for the summary of analyses for the different programs
deg.summary <- data.frame() ## Summary table for differentially expressed genes (DEG)
enrich.summary <- data.frame() ## Summary table for functional enrichment analysis
genesets <- list()
## Iterate over DEG analysis tools
deg.tool <- "edgeR" ## Default tool for testing
for (deg.tool in deg.tools) {

  ## Load the DEG analysis result table
  deg.table <- read.table(parameters[[deg.tool]], row.names=1)
  #dim(deg.table)
  
  ## Select the threshold column according to the DEG tool.
  ## Give similar names to similar columns.
  if (deg.tool == "DESeq2") {
    names(deg.table)[names(deg.table) == "log2FoldChange"] = "log2FC"
  } else if (deg.tool == "edgeR") {
    #threshold.column <- "FDR"
    names(deg.table)[names(deg.table) == "FDR"] = "padj"
    names(deg.table)[names(deg.table) == "logFC"] = "log2FC"
  }
  threshold.column <- "padj" ## Select genes based on adjusted p-values
  deg.table$FC <- 2^deg.table$log2FC

  ## Select genes according to the user-specified threshold
  deg.table$positive <- deg.table[,threshold.column] < parameters$deg.threshold
  deg.table$positive[is.na(deg.table$positive)] <- FALSE ## Replace the NA values by FALSE to avoid problems with selection
  
  ## Select up- and down-regulated genes
  deg.table$up <- deg.table$positive & deg.table$log2FC > 0
  deg.table$down <- deg.table$positive & deg.table$log2FC < 0
  
  ## Summarize the results
  deg.summary <- rbind(
    deg.summary,
    data.frame("tool" = deg.tool,
               "threshold column" = threshold.column,
               "deg.threshold" = parameters$deg.threshold,
               "Up-regulated" = sum(deg.table$up),
               "Down-regulated" = sum(deg.table$down),
               "Positives" = sum(deg.table$positive == TRUE),
               "Negatives" = sum(deg.table$positive == FALSE) - sum(is.na(deg.table[,threshold.column])),
               "NA" = sum(is.na(deg.table[,threshold.column])),
               "nb.genes" = nrow(deg.table)))



  ## Print out the summary table for the Rmd report
  # kable(data.frame(t(deg.summary)), format = "markdown") 
  
  ## Convert gene names into Entrez IDs
  gg <- bitr(row.names(deg.table), fromType="SYMBOL", toType="ENTREZID", annoDb="org.Dm.eg.db")
  ## Add a column with the "entrez.id" column in the deg.table
  deg.table[gg[,1], "entrez.id"] = gg[,2]
  
  ## Run GO enrichment analysis for up-regulated, down-regulated and the union of them.
  groups <- c("up", "down", "positive")  
  group <- "positive" ## Default group for testing
  ## Iterate over groups  
  for (group in groups) {

    ## Iterate over selection criteria (significant genes versus random selections)
    selection.criteria <- c("signif", "random")
    selection.criterion <- "signif"
    for (selection.criterion in selection.criteria) {
      
      ## Selet significant genes
      gene.ids <- deg.table[deg.table[,group], "entrez.id"] ## IDs of the significant genes for the group of interest
      
      if (selection.criterion == "random") {
        ## Select random gene list of the same size as the significant genes
        gene.ids <- sample(x = deg.table$entrez.id, size = length(gene.ids))  
      }
      gene.ids <- na.omit(gene.ids)   ## Filter out the genes for which no EntrezID has been found
      
      ## Iterate over ontology types (biological process, cellular component, molecular function)
      ontologies <- c("BP", "CC", "MF")
      ontology <- "MF" ## Select biological processes
      for (ontology in ontologies) {
        message(paste(sep=" ", "Enrichment analysis for", 
                      ontology, "classfication; ",
                      "DEG set", deg.tool, 
                      length(gene.ids), group, selection.criterion, "genes"))
        
        
        ## Compute the enrichment of the gene list for Biological Processes of the Gene Ontology
        ego <- enrichGO(gene = gene.ids, organism = "fly",
                        pvalueCutoff = parameters$enrich.threshold,
                        qvalueCutoff = parameters$enrich.threshold,
                        ont = ontology, readable = T)
        # nrow(ego@result) ## Count the number of significant classes
        
        
        ## Export the result table
        if (selection.criterion == "random") {
          ego.prefix <- "random"
        } else {
          ego.prefix <- comparison
        }
        ego.prefix <- paste(
          sep="", ego.prefix, "_", group, 
          "_", deg.tool, "_padj", parameters$deg.threshold,
          "_enrichGO_", ego@ontology, "_", ego@pAdjustMethod, ego@pvalueCutoff)
        ego.path <- file.path(dir.enrichment, ego.prefix)
        
        ## Save the result table (one row per significant biological process)
        write.table(ego@result, file=paste(sep="", ego.path, ".tab"), 
                    quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
      
        enrich.summary <- rbind(
          enrich.summary, 
          data.frame(
            prefix=ego.prefix,
            deg.tool=deg.tool,
            threshold.column = threshold.column,
            threshold.value = parameters$deg.threshold,
            select.type = selection.criterion,
            nb.genes = length(ego@gene),
            nb.genes.universe = length(ego@universe),
            ontology = ego@ontology,
            pAdjustMethod = ego@pAdjustMethod,
            pvalueCutoff = ego@pvalueCutoff,
            nb.signif.classes = nrow(ego@result)
          )
        )
        
        ## Barplot with the enrichGO results
        #         png(file=paste(sep="", ego.path, "_barplot.png"), width=10*150, heigh=7*150)
        #         barplot(ego,  showCategory=parameters$barplot_catnb, 
        #                 main=ego.prefix)
        #         tais.toi <- dev.off()

        max.to.plot <- min(nrow(ego@result), parameters$barplot_catnb)
        par.ori <- par()
        pdf(file=paste(sep="", ego.path, "_barplot.pdf"), width=8, height=10)
        par(mar= c(4, 15,4,1))
        barplot(rev(ego@result[1:max.to.plot, "Count"]),
                names.arg = rev(ego@result[1:max.to.plot, "Description"]),
                horiz = TRUE, las=1, cex.names=0.7, width=0.7, cex.main=0.9,
                main=ego.prefix, xlab="Nb genes")
        tais.toi <- dev.off()
        
#         plot(ego@result[,c("p.adjust", "Count")])
#         
#         pdf(file=paste(sep="", ego.path, "_barplot.pdf"), width=10, heigh=7)
#         #      X11(width=10, height=7)
#         barplot(ego, showCategory=parameters$barplot_catnb)
#         tais.toi <- dev.off()
#         
        ## Compute the enrichment of the gene list for Biological Processes
        ## of the Gene Ontology using groupGO, which permits to specify a minimal 
        ## level for the GO annotations. This avoids to have significant enrichment 
        ## for very basic classes (e.g. biological process).
        #     ggo <- groupGO(gene = gene.ids, organism = "fly",
        #                    level=parameters$groupGo.level,
        #                    ont = ontology, readable = T)
        #     
        #     ## Export the result table
        #     ggo.path <- file.path(
        #       dir.enrichment, 
        #       paste(sep="", comparison, "_", group, 
        #             "_", deg.tool, "_padj", parameters$deg.threshold,
        #             "_groupGO_", ggo@ontology, "_level", ggo@level))
        #     write.table(ggo@result, file=paste(sep="", ggo.path, ".tab"), 
        #                 quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
        #     
        #     ## Export a barplot with the groupGO results
        #     pdf(file=paste(sep="", ggo.path, "_barplot.pdf"), width=10, heigh=7)
        #     barplot(ggo,  showCategory=parameters$barplot_catnb)
        #     tais.toi <- dev.off()
        
        
        ## TO DO LATER: Gene Set Enrichment Analysis !!!
        ## ego2 <- gseGO(geneList = ggo,organism = "fly",ont = "BP" )
      }
    }
  } 
}