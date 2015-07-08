################################################################
## R code to run differential analysis
##
## Authors: Jeanne Chèneby & Justine Long
################################################################

library('DESeq2')

# directory<-"results/DEG/"

## Importing R_params file
r.params.path <- commandArgs(trailingOnly = FALSE)[6]
source(r.params.path)
dir.main <- "."

for (i in 1:length(comparisons)){

## Create a specific result directory for this differential analysis
dir.results <- file.path(dir.main, "results", "DEG", paste(sep="", comparisons.cond1[i], "_vs_", comparsions.cond2[i]))
dir.create(path = dir.results, showWarnings = FALSE, recursive = TRUE)
dir.figures <- file.path(dir.results, "figures")
dir.create(path = dir.figures, showWarnings = FALSE, recursive = TRUE)

sampleFiles <- counts.files.per.comparisons[[i]]

# condition.per.comparisons  <- c('N', 'N', 'N','NN','NN','NN') 
 
sampleCondition <- condition.per.comparisons[[i]]

sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)


ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=dir.main, design=~condition)
# ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, design=~condition)


colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=comparisons[[i]])



dds<-DESeq(ddsHTSeq)
res<-results(dds)
res<-res[order(res$padj),]
head(res)


# dev.copy(png, paste(sep="", data.root, "DEG/", comparisons.cond1[i], "_VS_", comparsions.cond2[i], "_deseq2_MAplot.png"))
pdf(file=file.path(dir.figures, paste(sep="", comparisons.cond1[i], "_VS_", comparsions.cond2[i], "_MAplot_DESeq2.pdf")))
plotMA(dds,ylim=c(-2,2),main="DESeq2")
dev.off()


mcols(res,use.names=TRUE)

# Now we want to transform the raw discretely distributed counts so that we can do clustering.
rld <- rlogTransformation(dds, blind=TRUE)

# dev.copy(png, paste(sep="", data.root, "DEG/", comparisons.cond1[i], "_VS_", comparsions.cond2[i], "deseq2_pca_DESeq2.png"))
pdf(file=file.path(dir.figures, paste(sep="", comparisons.cond1[i], "_VS_", comparsions.cond2[i], "_pca_DESeq2.pdf")))
print(plotPCA(rld, intgroup=c("condition")))
dev.off()

col_names <- colnames(res)
col_names[1] <- paste(sep='', "#gene_id\t", col_names[1])
colnames(res) <- col_names
print(colnames(res))
print(as.data.frame(res)[1,])

## adding "_DESeq2" to the output files
temp_output <- strsplit(output[i], "[.]")
final_output <- paste(sep='', temp_output[[1]][1], "_DESeq2.", temp_output[[1]][2])
#write the table to a tab delimited file
write.table(as.data.frame(res), final_output, sep="\t")


}