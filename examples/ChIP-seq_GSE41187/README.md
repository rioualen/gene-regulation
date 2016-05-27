# Study case for the development of a Bacterial ChIP-seq workflow


## Data source

Reference: Myers KS, Yan H, Ong IM, Chung D et al. Genome-scale analysis of escherichia coli FNR reveals complex features of transcription factor binding. PLoS Genet 2013 Jun;9(6):e1003565. PMID: [23818864](http://www.ncbi.nlm.nih.gov/pubmed/23818864)


GEO series: [GSE41187](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41187)

<!--
mkdir -p ~/mydisk/GSM41187_ChIP-seq/source/GSM1010224 ~/mydisk/GSM41187_ChIP-seq/source/GSM1010219 
cd ~/mydisk/GSM41187_ChIP-seq/source/GSM1010224
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX189%2FSRX189778/SRR576938/SRR576938.sra
cd ~/mydisk/GSM41187_ChIP-seq/source/GSM1010219
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX189%2FSRX189773/SRR576933/SRR576933.sra
-->

## File organisation

The short reads were downloaded from SRA ftp site using the links in the Excel file samples.xlsx. In this dataset, there is a single run per sample, and thus a single sra file per GSM ID. 

Short read files (.sra) are stored in sample-specific folders named data/source/[GSM_ID].


