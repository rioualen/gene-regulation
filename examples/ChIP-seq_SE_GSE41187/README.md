# Study case for the development of a Bacterial ChIP-seq workflow


## Data source

Reference: Myers KS, Yan H, Ong IM, Chung D et al. Genome-scale analysis of escherichia coli FNR reveals complex features of transcription factor binding. PLoS Genet 2013 Jun;9(6):e1003565. PMID: [23818864](http://www.ncbi.nlm.nih.gov/pubmed/23818864)


GEO series: [GSE41187](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41187)



## File organisation

Short read files (.sra) are stored in sample-specific folders named data/source/[GSM_ID].
In this dataset, there is a single run per sample, and thus a single sra file per GSM ID. 


## Setup analysis environment

export ANALYSIS_DIR=/your/preferred/dir/ChIP-seq_SE_GSE41187

mkdir -p ${ANALYSIS_DIR}
mkdir -p ${ANALYSIS_DIR}/data 
mkdir -p ${ANALYSIS_DIR}/genome

cd ${ANALYSIS_DIR}
git clone https://github.com/rioualen/gene-regulation.git
<!--TODO: replace with tar.gz download of gene-regulation-v3.0-->

## Download data

mkdir -p ${ANALYSIS_DIR}/data/GSM1010224 ${ANALYSIS_DIR}/data/GSM1010219 ${ANALYSIS_DIR}/data/GSM1010220
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX189%2FSRX189778/SRR576938/SRR576938.sra -P ${ANALYSIS_DIR}/data/GSM1010224
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX189%2FSRX189773/SRR576933/SRR576933.sra -P ${ANALYSIS_DIR}/data/GSM1010219
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX189/SRX189774/SRR576934/SRR576934.sra -P ${ANALYSIS_DIR}/data/GSM1010220


## Download reference genome & annotations

wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/fasta/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna.genome.fa.gz -P ${ANALYSIS_DIR}/genome
wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gff3/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gff3.gz -P ${ANALYSIS_DIR}/genome
wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gtf/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gtf.gz -P ${ANALYSIS_DIR}/genome
gunzip ${ANALYSIS_DIR}/genome/*.gz


## Execute workflow

snakemake -s ${ANALYSIS_DIR}/gene-regulation/scripts/snakefiles/workflows/ChIP-seq_workflow_SE.py -p -j 5 --configfile ${ANALYSIS_DIR}/gene-regulation/examples/ChIP-seq_SE_GSE41187/config.yml
