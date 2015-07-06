################################################################
## R code to run differential analysis
##
## Authors: Jeanne Chèneby & Justine Long
################################################################

library("edgeR")
library("limma")

## setwd("../../") ## Only for testing when the test starts from the scripts directory

## Temporary: set up manually the files and folder.
## We will later improve this by reading the config file from arguments (argv) 
r.params.path <- commandArgs(trailingOnly = FALSE)[6]
source(r.params.path)
dir.main <- "."

for (i in 1:length(comparisons)){
	
## Create a specific result directory for this differential analysis
dir.results <- file.path(dir.main, "results", "DEG", paste(sep="", comparisons.cond1[i], "_vs_", comparsions.cond2[i]))
dir.create(path = dir.results, showWarnings = FALSE, recursive = TRUE)
dir.figures <- file.path(dir.results, "figures")
dir.create(path = dir.figures, showWarnings = FALSE, recursive = TRUE)

## Read data from counts files

counts <- readDGE(counts.files.per.comparisons[[i]])$counts

## Remove summary lines from HTseq files
# noint <- rownames(counts) %in% c("no_feature","ambiguous","too_low_aQual", "not_aligned","alignment_not_unique")
info.rows <- rownames(counts) %in% grep(rownames(counts),
                                        pattern = "__", 
                                        value = TRUE)

## Remove all features that have less than 1 reads per millions of reads
cpms <- cpm(counts)    ## Counts per million reads


## Only keep genes detected in at least n.rep samples, where n.rep is the number of replicates
keep <- (rowSums(cpms > 1) >= n.rep[i]) & !info.rows  
counts <- counts[keep,]

# Preparation of tables
colnames(counts) <- names.per.comparisons[[i]]

########
## Convert the count table in a DGEList structure
d <- DGEList(counts=counts, group=condition.per.comparisons[[i]])
d <- calcNormFactors(d)

## Export MDS plot
pdf(file=file.path(dir.figures, paste(sep="", "MDS_plot_", comparisons.cond1[i], "_vs_", comparsions.cond2[i], "_edgeR.pdf")))
plotMDS(d, labels=names.per.comparisons[[i]], col=c("darkgreen","blue")[factor(condition.per.comparisons[[i]])])
dev.off()

d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)

pdf(file= file.path(dir.figures, paste(sep = "", "plotMeanVar_", comparisons.cond1[i], "_VS_", comparsions.cond2[i], "_edgeR.pdf")))
plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
dev.off()

pdf(file= file.path(dir.figures, paste(sep = "", "plotBCV_", comparisons.cond1[i], "_VS_", comparsions.cond2[i], "_edgeR.pdf")))
plotBCV(d)
dev.off()

## Detect differentially expressed genes by applying the exact test
de <- exactTest(d, pair = c(comparisons[[i]]))

## Tabular summary of the DE stats
tt <- topTags(de, n=nrow(d))

nc <- cpm(d, normalized.lib.sizes=TRUE)
rn <- rownames(tt$table)

deg <- rn[tt$table$FDR < FDR.threshold]

pdf(file=file.path(dir.figures, paste(sep = "", "plotSmear_", comparisons.cond1[i], "_VS_", comparsions.cond2[i], "_edgeR.pdf")))
plotSmear(d, de.tags = deg)
dev.off()


## TEMP output <- sub(pattern="data/1258-BRM", replacement = dir.results, output)
# write.csv(tt$table, output[i])

## adding "_edgeR" to the output files
temp_output <- strsplit(output[i], "[.]")
final_output <- paste(sep='', temp_output[[1]][1], "_edgeR.", temp_output[[1]][2])

write.table(tt$table, final_output, sep="\t")
# write.table(tt$table, output[i], sep="\t")
# write.delim(tt$table, output[i], row.names = TRUE, sep = "\t")
}