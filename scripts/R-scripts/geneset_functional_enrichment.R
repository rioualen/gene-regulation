library("knitr")
library('clusterProfiler')
library(GO.db)
#library("org.Dm.eg.db")




## Parameters
dir.main <- "/home/lkhamvongsa/mountrsatlocal/dr-chip-rna-seq"
setwd(dir.main)
## Source : M. Carlson. How To Use GOstats and Category to do Hypergeometric testing with unsupported model organisms. 
org <- "dme"
# org <- "eco"
if (org == "eco") {
  organism <- "Escherichia coli K12 MG1655"
  #annot.schema <- "ECOLI_DB"
  org.db <- "org.EcK12.eg.db"
} else if (org =="dme") {
  organism <- "Drosophila melanogaster"
  annot.schema <- "FLY_DB"
  org.db <- "org.Dm.eg.db"
}

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
  all.entrez.ids <- deg.table$entrez.id[!is.na(deg.table$entrez.id)]
  names(all.entrez.ids) <- row.names(deg.table)[!is.na(deg.table$entrez.id)]
  # length(entrez.ids)
  
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
        gene.ids <- sample(x = deg.table$entrez.id, size = sum(deg.table[,group]))  
      }
      gene.ids <- gene.ids[!is.na(gene.ids)] ## Filter ou unidentified genes
      #      gene.ids <- na.omit(gene.ids)   ## Filter out the genes for which no EntrezID has been found
      
      ## Iterate over ontology types (biological process, cellular component, molecular function)
      ontologies <- c("BP", "CC", "MF")
      ontology <- "MF" ## Select biological processes
      for (ontology in ontologies) {
        message(paste(sep=" ", "Enrichment analysis for", 
                      ontology, "classfication; ",
                      "DEG set", deg.tool, 
                      length(gene.ids), group, selection.criterion, "genes"))
        
        ## Just a wrapper around various functional enrichment packages
        enrich.result <- functional.enrichment(geneset=gene.ids, 
                                               allgenes=all.entrez.ids, 
                                               db=org.db, 
                                               ontology="BP",
                                               thresholds = c("evalue"=1, "qvalue"=0.05),
                                               select.positives=FALSE,
                                               run.GOstats = TRUE,
                                               run.clusterProfiler = TRUE,
                                               clusterProfiler.org="fly",
                                               plot.adjust = TRUE,
                                               verbosity=1)
        
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

        
        ## Prepare a result table with some added columns for convenience of plot generation
        ego.result.table <- ego@result
        
        ## Barplot with the enrichGO results. 
        max.to.plot <- min(nrow(ego@result), parameters$barplot_catnb)
        pdf(file=paste(sep="", ego.path, "_barplot.pdf"), width=8, height=10)
        par.ori <- par()
        par(mar= c(4, 15,4,1))
        barplot(rev(ego@result[1:max.to.plot, "Count"]),
                names.arg = rev(ego@result[1:max.to.plot, "Description"]),
                horiz = TRUE, las=1, cex.names=0.7, width=0.7, cex.main=0.9,
                main=ego.prefix, xlab="Nb genes")
        silence <- dev.off()
        
        ## The customized barplot() for enrichGO return object has a bug: 
        ## when it is inserted in a loop, the exported pdf or png files are empty.
        #         png(file=paste(sep="", ego.path, "_barplot.png"), width=10*150, heigh=7*150)
        #         barplot(ego,  showCategory=parameters$barplot_catnb, 
        #                 main=ego.prefix)
        #         tais.toi <- dev.off()
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


################################################################
## Build custom annotation table

## Source : M. Carlson. How To Use GOstats and Category to do Hypergeometric testing with unsupported model organisms. 
org <- "dme"
# org <- "eco"
if (org == "eco") {
  organism <- "Escherichia coli K12 MG1655"
  #annot.schema <- "ECOLI_DB"
  org.db <- "org.EcK12.eg.db"
} else if (org =="dme") {
  organism <- "Drosophila melanogaster"
  annot.schema <- "FLY_DB"
  org.db <- "org.Dm.eg.db"
}

library(org.db, character.only = TRUE)

annot.prefix <- sub(org.db, pattern = ".db$", replacement = "")
annot.db.go <- paste(sep="", annot.prefix, "GO")
annot.db.kegg <- paste(sep="", annot.prefix, "PATH")


## ??? Do we need the following piece of code ? It seems to be only valid for Homo sapiens
library("AnnotationForge")
# available.dbschemas()
frame = toTable(get(annot.db.go))
frame$Evidence[frame$Evidence=="-"] <- "ND"
goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)
# View(goframeData)
# dim(goframeData)

## Prepare GO gene mapping
goFrame=GOFrame(goframeData,organism=organism)
goAllFrame=GOAllFrame(goFrame)

## Generate a geneSetCollection object
library("GSEABase")
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

## Setting up the parameter object
library("GOstats")
universe = Lkeys(get(annot.db.go))
# length(universe)
genes = universe[1:100] ## Quick test with arbitrary genes
params <- GSEAGOHyperGParams(
  name=paste(goFrame@organism, "GSEA based annot Params"),
  geneSetCollection=gsc,
  geneIds = genes,
  universeGeneIds = universe,
  ontology = "BP",
  pvalueCutoff = 0.05,
  conditional = FALSE,
  testDirection = "over")
Over <- hyperGTest(params)
head(summary(Over))

## Preparing KEGG to gene mapping
kegg.frame = toTable(get(annot.db.kegg))
keggframeData = data.frame(kegg.frame$path_id, kegg.frame$gene_id)
# head(keggframeData)
keggFrame=KEGGFrame(keggframeData,organism=organism)
gsc <- GeneSetCollection(keggFrame, setType = KEGGCollection())
universe = Lkeys(get(annot.db.kegg))
genes = universe[1:100]

kparams <- GSEAKEGGHyperGParams(
  name=paste(goFrame@organism, "GSEA based annot Params"),
  geneSetCollection=gsc,
  geneIds = genes,
  universeGeneIds = universe,
  pvalueCutoff = 0.05,
  testDirection = "over")

kOver <- hyperGTest(kparams)
head(summary(kOver))

toLatex(sessionInfo())

## ROntoTool Vignette
require(graph)
require(ROntoTools)
korg <- "eco" ## Organism abbreviation in KEGG
kpg <- keggPathwayGraphs(korg, verbose = FALSE)
#kpg <- keggPathwayGraphs(korg, updateCache = TRUE, verbose = TRUE)  ## This may be slow
head(names(kpg))

kpg <- setEdgeWeights(kpg, edgeTypeAttr = "subtype",
                      edgeWeightByType = list(activation = 1, inhibition = -1,
                                              expression = 1, repression = -1),
                      defaultWeight = 0)

load(system.file("extdata/E-GEOD-21942.topTable.RData", package = "ROntoTools"))
head(top)
fc <- top$logFC[top$adj.P.Val <= .01]
names(fc) <- top$entrez[top$adj.P.Val <= .01]
pv <- top$P.Value[top$adj.P.Val <= .01]
names(pv) <- top$entrez[top$adj.P.Val <= .01]
head(fc)
head(pv)
fcAll <- top$logFC
names(fcAll) <- top$entrez
pvAll <- top$P.Value
names(pvAll) <- top$entrez
ref <- top$entrez
kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)
head(nodeWeights(kpg[[names(kpg[1])]]))
peRes <- pe(x = fc, graphs = kpg, ref = ref,  nboot = 200, verbose = FALSE)
head(Summary(peRes))


## ReactomePA vignette
require(DOSE)
data(geneList)
de <- names(geneList)[abs(geneList) > 1.5]
head(de)


require(ReactomePA)
x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
head(summary(x))
barplot(x, showCategory=36)
dotplot(x, showCategory=36)
enrichMap(x)
cnetplot(x, categorySize="pvalue", foldChange=geneList)

require(clusterProfiler)
data(gcSample)
res <- compareCluster(gcSample, fun="enrichPathway")
plot(res)


## Gene Set Enrichment Analysis (GSEA)
y <- gsePathway(geneList, nPerm=100, 
                minGSSize=120, pvalueCutoff=0.05, 
                pAdjustMethod="BH", verbose=FALSE)
res <- summary(y)
head(res)
enrichMap(y)

gseaplot(y, geneSetID = "1280215")
viewPathway("E2F mediated regulation of DNA replication", readable=TRUE, foldChange=geneList)

