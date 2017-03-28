#' @title Descriptive statistics from a read count table.
#' @author Jacques van Helden
#' @description Compute descriptive statistics for each sample of a 
#' RNA-seq count table. The table contains one row per feature 
#' (gene, transcript) and one column per sample. 
#' @param config.file a yaml-formatted file describing the parameters 
#' of the workflow. 
#' @param count.file file containing the counts of reads per feature, 
#' tab-separated value text file.
#' @export
RnaSeqExploreCounts <- function(config.file, count.file) {
  
  ## Initialize the result object
  result <- list()
  result$config.file <- config.file
  result$count.file <- count.file
  
  ## Read the configuration file
  library(yaml)
  param <- yaml.load_file(config.file)
  
  
  ################################################################
  ## Read the count table
  verbose("Loading count table", 2)
  counts <- read.delim(count.file, row.names=1, sep="\t")

  ## Filter out the rows corresponding to non-assigned counts, 
  ## e.g. __not_aligned, __ambiguous, __too_low_qAual, __not_aligned
  not.feature <- grep(rownames(counts), pattern = "^__")
  if (length(not.feature) > 0) {
    counts <- counts[-not.feature,]
  }
  
  ## Check that some required parameters are well defined
  ## A trick: to enable log-scaled plots for 0 values, I add an epsilon increment
  epsilon <- as.logical(param$deg$epsilon)
  if (!is.numeric(epsilon)) {
    param$deg$epsilon <- 0.1
  }
  
  
  ## Default method for the selection of the final list of DEG
  selectionCriterion <- param$deg$selectionCriterion
  if (is.null(selectionCriterion)) {
    param$deg$selectionCriterion <- "edgeR"
  }
  
  ## Check that selection criterion is valid
  acceptedSelectionCriteria <- c("edgeR", "DESeq2", "union", "intersection")
  if (!selectionCriterion %in% acceptedSelectionCriteria) {
    stop("Invalid selection criterion. Supported: ", 
            paste(acceptedSelectionCriteria, collapse=", "), ".")
  }
  
  gtf.file <- param$genome$gtf_file
  
  ################################################################
  ## Load gene information from the GTF file
  ## (should be the same as used to count tags per gene)
  if (!is.null(gtf.file)) {
    ## Jacques or Lucie, TO DO : check if the method makeTxDbFromGFF allows 
    ## to get gene names and descriptions.
    ## 
    ## In principle it should be possible since this info is in the GTF file. 
    ## If not, we migh rewite a parser for GTF. 
    
    message("Loading gene annotations from GTF file: ", gtf.file)
    gtf.file.path <- file.path(dir.main, param$dir$genome, gtf.file)
    library(rtracklayer)  # for the import() function
    #   gr <- import(gtf.file)
    txdb <- makeTxDbFromGFF(file=gtf.file.path)
    #organism=organism.name, ## OPTION organism DOES NOT WORK WITH DESULFO because not supported in this library. However the GTF is correctly loaded
    #                          dataSource=gtf.source)
    #seqlevels(txdb)
    #   transcripts <- transcriptsBy(txdb, by="gene")
    #   transcript.table <- as.data.frame(transcripts@unlistData)
    #   cds <- cdsBy(txdb, by="gene")
    #   cds.table <- as.data.frame(cds@unlistData)
    #   exons <- exonsBy(txdb, by="gene")
    #   exon.table <- as.data.frame(exons@unlistData)
    all.genes <- genes(txdb)
    gene.info <- as.data.frame(all.genes)
    gene.info$name <- gene.info$gene_id
    #  gene.info$entrez.id <- NA
    gene.info$description <- "no description"
  } else {
    verbose(paste("No GTF file has been specified"))
    g <- nrow(counts)
    gene.info <- data.frame("seqnames"=rep(NA, times=g),
                            "start"=rep(NA, times=g),
                            "end"=rep(NA, times=g),
                            "width"=rep(NA, times=g),
                            "strand"=rep(NA, times=g),
                            "gene_id"=row.names(counts),
                            "name"=row.names(counts),
                            #                          "entrez.id" = rep(NA, times=g),
                            "description"=rep("no description", times=g))
    row.names(gene.info) <- row.names(counts)
  }
    
  ## Logical variable indicating whether or not the data tables should be exported 
  ## in Excel format.
  export.excel.files <- as.logical(param$deg$export_xlsx)
  if (is.na(export.excel.files)) {
    param$deg$export_xlsx <- FALSE
  }

  ## Output directory for the differential expression analysis
  dir.diffexpr <- param$dir$diffexpr
  if (is.null(dir.diffexpr))  {
    param$dir$diffexpr <- file.path("results", "diffexpr")
  }
  ## Prefix and suffix for output files concerning the whole count table (all samples together)
  ## prefix["general.file"] <- sub(pattern = ".tab", replacement="", count.file)
  suffix <- vector()
  suffix["DEG"] <- 'DEG' ## Suffix for the output files
  suffix["DESeq2"] <- paste(sep="_", suffix["DEG"], "DEG_DESeq2")
  suffix["edgeR"] <- paste(sep="_", suffix["DEG"], "DEG_edgeR")
  prefix <- vector()
  prefix["general.file"] <- file.path(dir.diffexpr, suffix["DEG"]) ## Path prefix for the general files
  
  
  ## Check that thresholds have been defined
  if (is.null(param$deg$thresholds)) {
    param$deg$thresholds <- list()
  } 
  if (is.null(param$deg$thresholds$padj)) {
    param$deg$thresholds$padj <- 0.05 # Threshold on adjusted p-value
  }
  if (is.null(param$deg$thresholds$FC)) {
    param$deg$thresholds$FC <- 1  # If not specified, do not impose any threshold on the fold-change
  }
  if (is.null(param$deg$thresholds$outliers)) {
    param$deg$thresholds$outliers <- 8.5 ## Somewhat arbitrary threshold to discard outliers
  }
  
  ## Read the sample description file, which indicates the 
  ## condition associated to each sample ID.
  sample.description.file <- file.path(dir.main, param$metadata$samples)
  message("Reading sample description file: ", sample.description.file)
  sample.desc <- read.delim(sample.description.file, sep="\t", 
                            comment=";", header = TRUE, row.names=1)
  sample.ids <- row.names(sample.desc)
  
  ## Experimental conditions
  sample.conditions <- as.vector(sample.desc[,1]) ## Condition associated to each sample
  names(sample.conditions) <- sample.ids
  # print(sample.conditions)
  
  ## Build sample labels by concatenating their ID and condition
  sample.desc$label <- paste(sep="_", sample.ids, sample.conditions)
  
  ## Define a specific color for each distinct condition
  conditions <- unique(sample.conditions) ## Set of distinct conditions
  cols.conditions <- brewer.pal(max(3, length(conditions)),"Dark2")[1:length(conditions)]
  names(cols.conditions) <- conditions
  # print(cols.conditions)
  
  ## Define a color per sample according to its condition
  sample.desc$color <- cols.conditions[sample.conditions]
  # names(cols.samples) <- sample.ids
  # print(cols.samples)
  
  
  ## Check that the header of counts match the sample IDs
  ids.not.found <- setdiff(sample.ids, names(counts)) ## Identify sample IDs with no column in the count table
  if ((!is.null(ids.not.found)) && (length(ids.not.found) > 0)) {
    stop(paste(length(ids.not.found), "Missing columns in count table", count.file), paste(sep="; ", ids.not.found))
  }
  
  ################################################################
  ## Restrict the count table to the sample IDs found in the sample description file
  counts <- counts[, sample.ids]
  # names(counts)
  # dim(counts)
  
  ########################################################################
  ## Treatment of 0 values.
  ## Add an epsilon to 0 values only, in order to enable log-transform and display on logarithmic axes.
  verbose(paste("Treating zero-values by adding epsilon =", epsilon), 2)
  counts.epsilon <- counts
  counts.epsilon[counts==0] <- epsilon
  
  ## Log-transformed data for some plots. 
  counts.log10 <- log10(counts.epsilon)
  counts.log2 <- log2(counts.epsilon)
  
  
  ################################################################
  ## Compute sample-wise statistics on mapped counts
  ################################################################
  stats.per.sample <- calc.stats.per.sample(sample.desc, counts)
  
  
  ################################################################
  ## Compute the counts per million reads 
  ################################################################
  #message("Computing CPMs")
  ## Note: the default normalization criterion (scaling by libbrary sum) 
  ## is questionable because it is stronly sensitive to outliers 
  ## (very highly expressed genes).  A more robust normalisation criterion 
  ## is to use the 75th percentile, or the median. We use the median, somewhat arbitrarily, 
  ## beause it gives a nice alignment on the boxplots.
  cpms.libsum <- cpm(counts.epsilon)    ## Counts per million reads, normalised by library sum
  cpms.perc75 <- cpm(counts.epsilon, lib.size = stats.per.sample$perc75)    ## Counts per million reads, normalised by 75th percentile
  cpms.perc95 <- cpm(counts.epsilon, lib.size = stats.per.sample$perc95)    ## Counts per million reads, normalised by 95th percentile
  cpms.median <- cpm(counts.epsilon, lib.size = stats.per.sample$median)    ## Counts per million reads, normalised by sample-wise median count
  cpms <- cpms.median ## Choose one normalization factor for the CPMs used below
  # cpms <- cpms.perc75 ## Choose one normalization factor for the CPMs used below
  cpms.log10 <- log10(cpms) ## Log-10 transformed CPMs, with the epsilon for 0 counts
  cpms.log2 <- log2(cpms) ## Log-10 transformed CPMs, with the epsilon for 0 counts
  
  ## Compute Trimmed Means of M Values (TMM): TO BE DONE
  stats.per.sample$cpm.mean <- apply(cpms, 2, mean)
  stats.per.sample$log2.cpm.mean <- apply(cpms.log2, 2, mean)
  stats.per.sample$log10.cpm.mean <- apply(cpms.log10, 2, mean)
  
  ## Detect outliers, i.e. genes with a very high number of reads (hundreds of thousands), most of which result from problems with ribodepletion.
  outliers <- (apply(cpms.log10, 1, max) > param$deg$thresholds$outliers)
  # rownames(cpms.log10[outliers,])
  # sum(outliers)
  
  
  result$param <- param ## Include the parameters in the result
  result$gene.info <- gene.info
  result$prefix <- prefix
  result$suffix <- suffix
  result$sample.desc <- sample.desc
  result$sample.ids <- sample.ids
  result$conditions <- conditions
  result$counts <- counts
  result$counts.epsilon <- counts.epsilon
  result$counts.log10 <- counts.log10
  result$counts.log2 <- counts.log2
  result$cpms <- cpms
  result$cpms.log10 <- cpms.log10
  result$cpms.log2 <- cpms.log2
  result$stats.per.sample <- stats.per.sample
  result$outliers <- outliers
  return(result)
}