# Study case for the development of a Bacterial RNA-seq workflow


## Data data

Reference: Myers KS, Yan H, Ong IM, Chung D et al. Genome-scale analysis of escherichia coli FNR reveals complex features of transcription factor binding. PLoS Genet 2013 Jun;9(6):e1003565. PMID: [23818864](http://www.ncbi.nlm.nih.gov/pubmed/23818864)


GEO series: [GSE41190](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41190)

<!--
mkdir -p data/GSM1010244 data/GSM1010246 data/GSM1010245 data/GSM1010247
cd data/GSM1010244
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059774/SRR191809/SRR191809.sra
cd data/GSM1010246
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059767/SRR191812/SRR191812.sra
cd data/GSM1010245
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059768/SRR191805/SRR191805.sra
cd data/GSM1010247
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX116%2FSRX116381/SRR400301/SRR400301.sra
-->

## File organisation

The short reads were downloaded from SRA ftp site using the links in the Excel file samples.xlsx. In this dataset, there is a single run per sample, and thus a single sra file per GSM ID. 

Short read files (.sra) are stored in sample-specific folders named data/[GSM_ID].

# Running the example workflow

## Downloading the data

```
##================================================================
## Create the base directory for the workflow 
## (can be adapted to your local configuration)
export TESTBASE=~/mydisk/GSM41190_RNA-seq/
mkdir -p ${TESTBASE}


##================================================================
## Download samples
##
## BEWARE: the samples are heavy (>1Gb per sample). 

cd ${TESTBASE}; mkdir -p data/GSM1010244
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059774/SRR191809/SRR191809.sra -P data/GSM1010244

cd ${TESTBASE}; mkdir -p data/GSM1010246
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059767/SRR191812/SRR191812.sra -P data/GSM1010246

cd ${TESTBASE}; mkdir -p data/GSM1010245
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059768/SRR191805/SRR191805.sra -P data/GSM1010245

## The last sample has formating problems related to an ambiguity between 
## single and paired-end, we temporarily ignore it
# cd ${TESTBASE}; mkdir -p data/GSM1010247
# wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX116%2FSRX116381/SRR400301/SRR400301.sra -P data/GSM1010247

```

## Running the workflow on the command line

```
## Go to the base directory
cd ${TESTBASE}; ln -fs ~/gene-regulation .

## Check that the worklow is consistent by printing all commands without running them (option -n)
snakemake -p -s gene-regulation/scripts/snakefiles/workflows/rna-seq_workflow_pe.py  --configfile gene-regulation/examples/RNA-seq_GSE41190/RNA-seq_GSE41190.yml -n

## Run the workflow
snakemake -p -s gene-regulation/scripts/snakefiles/workflows/rna-seq_workflow_pe.py  --configfile gene-regulation/examples/RNA-seq_GSE41190/RNA-seq_GSE41190.yml

```

## Running the workflow with a job scheduler

Snaemake automatically manage the distribution of taks on a multiprocessor computer or/and a PC cluster via the *qsub* job scheduler. 


```
## Go to the base directory
cd ${TESTBASE}; ln -fs ~/gene-regulation .

## Check that the worklow is consistent by printing all commands without running them (option -n)
snakemake -p -c "qsub {params.qsub}" -j 12 \
  -s gene-regulation/scripts/snakefiles/workflows/rna-seq_workflow_pe.py \
  --configfile gene-regulation/examples/RNA-seq_GSE41190/RNA-seq_GSE41190.yml

```


