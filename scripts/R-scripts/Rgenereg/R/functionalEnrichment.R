#' @title Run functional enrichment for a given gene set
#'
#' @author Jacques van Helden (\email{Jacques.van-Helden@@univ-amu.fr})
#'
#' @description Run functional enrichment for a given gene set.
#' This function is a compilation of pieces of code from the vignettes
#' of different functional enrichment analysis package.
#'
#' @details
#' First version: 2015-08
#' Last modification: 2015-08.
#'
#' @param geneset   Entrez IDs of the gene set of interest.
#' @param allgenes  EntrezIDs for all genes of the analyzed genome (the "universe").
#' @param db an annotation database
#' @param thresholds["evalue"]=0.1 Upper threshold on the e-value for GOstats.
#' Note that we perform the selection with a p-value threshold of 1, in order to obtain the
#' stats for all the classes. The e-value is calculated a posteriori and indicated in the table.
#' @param select.positives=TRUE If TRUE, the result table only contains the positive tests,
#' i.e. those passing all the thresholds. Otherwise, all rows are returned, with a column "positive"
#' indicating whether or not they did or not pass all thresholds. Additional columns indicate the
#' result of the individual thresholds.
#' @param organism organism name in various tools. E.g:
#'
#' @examples
#'
#' @export
functionalEnrichment <- function(geneset,
                                 allgenes,
                                 db,
                                 organism.names,
                                 ontology="BP",
                                 thresholds = thresholds,
                                 select.positives=TRUE,
                                 run.GOstats = TRUE,
                                 run.clusterProfiler = FALSE,
                                 plot.adjust = TRUE) {
  #library("ALL")
  verbose(paste(sep="", "Go over-representation analysis. ",
                length(geneset), " input genes (among ", length(allgenes), ")"), 2)

  ## Prepare the result
  result <- list()
  result$ngenes.test <- length(geneset)
  result$ngenes.all <- length(allgenes)
  result$db = db
  result$ontology = ontology
  result$thresholds <- thresholds
  result$select.positives <- select.positives

  ## Create a data frame with results per gene (one row per gene, one column per attribute)
  result.per.gene <- data.frame("entrez.id" = allgenes)
  row.names(result.per.gene) <- allgenes
  result.per.gene$selected <- 0
  result.per.gene[geneset,"selected"] <- 1

  ## Load annotation library
  library("GO.db")
  library("annotate")
  library(db, character.only = TRUE)
  envir <- c(
    db=db,
    prefix = sub(db, pattern = ".db", replacement = ""))
  for (suffix in c("GO", "ENTREZID")) {
    envir[suffix] <- paste(sep="", envir["prefix"],suffix)
  }

  ################################################################
  ## Run over-representation analysis with GOstats
  if (run.GOstats) {
    library("GOstats")

    ## Get go annotations per gene
    go.annot.geneset <- mget(geneset, envir=get(envir["GO"]))
    go.annot.allgenes <- mget(allgenes, envir=get(envir["GO"]))
    result.per.gene$nb.go.annot <- as.vector(unlist(lapply(go.annot.allgenes, length)  ))

    ## Compute number of genes with at least one annotation in GO.
    ## Note that this is for all the GO ontologies. For the hypergeometric test
    ## below we need to take into account the ontology-specific annotations.
    nb.allgenes.with.go <- sum(result.per.gene$nb.go.annot >0) ## Total number of genes with GO annotations
    nb.geneset.with.go <- sum(result.per.gene[geneset, "nb.go.annot"] >0) ## Number of genes with GO annotations in the test set
    # head(result.per.gene)
    # table(result.per.gene$nb.go.annot)

    ## Define parameters for GOstats analysis
    gostats.params <- new("GOHyperGParams",
                          geneIds=geneset,
                          universeGeneIds=allgenes,
                          annotation=db,
                          ontology=ontology,
                          pvalueCutoff=1,
                          conditional=FALSE,
                          testDirection="over")

    ## Run over-representation hypergeometric test with GOstats
    over.result <- hyperGTest(gostats.params)
    # plot(goDag(over.result))

    result$geneIds <-  length(geneIds(over.result)) # Number of mapped genes
    result$geneIdUniverse <- length(geneIdUniverse(over.result))

    ## Collect the statistics for all the tested GO classes.
    ## For this we set the p-value threshold to 1.1, so we select
    ## all classes, even non-significant!
    go.enrich.table <- summary(over.result, pvalue=2)
    go.enrich.table$expectedCounts = expectedCounts(over.result)

    go.enrich.table <- complete.enrich.table(go.enrich.table,
                                             pvalue.column="Pvalue",
                                             thresholds=thresholds,
                                             select.positives=select.positives,
                                             plot.adjust=plot.adjust)

    ## Add  GO table to the result
    result$go.bp.table <- go.enrich.table


    ## For the sake of understanding, we can re-compute the p-value
    ## with the hypergeometric distribution. We obtain the same result as GOstats.
    N <- universeMappedCount(over.result) ## Universe, i.e. nb of genes with at least one annotation in the selected ontology
    m <- go.enrich.table$Size ## genes marked as belonging to the considered GO class
    n <- N -m ## Number of "non-marked" genes, i.e. not belonging to the considered GO class
    k <- length(over.result@geneIds) ## Number of test genes with at least one annotation in the selected ontology
    Pvalue.check <-
      phyper(q = go.enrich.table$Count -1, m=m, n=n, k = k, lower=FALSE)
    #   plot(go.enrich.table[,pvalue.column], Pvalue.check, log="xy")
    ExpCount.check <- m * k/N
    # plot(go.enrich.table$ExpCount, ExpCount.check)
    #   abline(a=0, b=1)

    # View(go.enrich.table)
    # names(go.enrich.table)
    # dim(go.enrich.table)

  }    ## End of GOstats analysis


  ## Run GO enrichment analysis with clusterProfiler
  if (run.clusterProfiler) {
    library("clusterProfiler")
    ## Check the organism
    if (is.na(organism.names["clusterProfiler"])) {
      message("Skipping clusterProfiler analysis because organism.names[\"clusterProfiler\"] is not specified")
    } else {
      ## Compute the enrichment of the gene list for Biological Processes of the Gene Ontology
      ## We intently set all cutoffs to a value > than the possible max, in order to get the full table.
      ## We then perform the filter only if requested.
      ego <- enrichGO(gene = geneset,
                      organism = organism.names["clusterProfiler"],
                      pvalueCutoff = 2,
                      qvalueCutoff = 2,
                      ont = ontology, readable = T)
      ego.table <- data.frame(ego@result)


      ego.table <- complete.enrich.table(
        ego.table,
        pvalue.column="pvalue",
        thresholds=thresholds,
        select.positives=select.positives,
        plot=TRUE)

      ## Add  GO table to the result
      result$clusterProfiler.bp.table <- ego.table
    }
  }

  ## Att the results per gene to the result object
  result$result.per.gene <- result.per.gene

  # plot(hist(go.enrich.table[,pvalue.column], breaks=20))
  #dim(go.enrich.table)
  # View(go.enrich.table)
  # print(over.result.cond)
  return(result)
}
