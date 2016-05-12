# Study case for the development of a Bacterial RNA-seq workflow


## Data source

Reference: Myers KS, Yan H, Ong IM, Chung D et al. Genome-scale analysis of escherichia coli FNR reveals complex features of transcription factor binding. PLoS Genet 2013 Jun;9(6):e1003565. PMID: [23818864](http://www.ncbi.nlm.nih.gov/pubmed/23818864)


GEO series: [GSE41190](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41190)

<!--
mkdir -p /data/source/RNA-seq_GSE41190/GSM1010244 /data/source/RNA-seq_GSE41190/GSM1010246 /data/source/RNA-seq_GSE41190/GSM1010245 /data/source/RNA-seq_GSE41190/GSM1010247
cd /data/source/RNA-seq_GSE41190/GSM1010244
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059774/SRR191809/SRR191809.sra
cd /data/source/RNA-seq_GSE41190/GSM1010246
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059767/SRR191812/SRR191812.sra
cd /data/source/RNA-seq_GSE41190/GSM1010245
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059768/SRR191805/SRR191805.sra
cd /data/source/RNA-seq_GSE41190/GSM1010247
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX116%2FSRX116381/SRR400301/SRR400301.sra
-->

## File organisation

The short reads were downloaded from SRA ftp site using the links in the Excel file samples.xlsx. In this dataset, there is a single run per sample, and thus a single sra file per GSM ID. 

Short read files (.sra) are stored in sample-specific folders named data/source/[GSM_ID].


