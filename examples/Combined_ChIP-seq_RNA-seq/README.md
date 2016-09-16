# Study case combined ChIP-seq & RNA-seq workflow on E. coli data

## Data source

Reference: Myers KS, Yan H, Ong IM, Chung D et al. Genome-scale analysis of escherichia coli FNR reveals complex features of transcription factor binding. PLoS Genet 2013 Jun;9(6):e1003565. PMID: [23818864](http://www.ncbi.nlm.nih.gov/pubmed/23818864)


GEO series: 
    [GSE41187](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41187)
    [GSE41190](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41190)

GEO Superseries:
    [GSE41195](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41195)


## Setup ChIP-seq analysis environment

export CHIP_DIR=/data/analyses/ChIP-seq_SE_GSE41187

mkdir -p ${CHIP_DIR}
mkdir -p ${CHIP_DIR}/data 
mkdir -p ${CHIP_DIR}/genome

cd ${CHIP_DIR}
<!--git clone https://github.com/rioualen/gene-regulation.git-->
<!--TODO: replace with tar.gz download of gene-regulation-v3.0-->
ln -s ~/Desktop/workspace/gene-regulation

### Download data

mkdir -p ${CHIP_DIR}/data/GSM1010224 ${CHIP_DIR}/data/GSM1010219 ${CHIP_DIR}/data/GSM1010220
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX189%2FSRX189778/SRR576938/SRR576938.sra -P ${CHIP_DIR}/data/GSM1010224
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX189%2FSRX189773/SRR576933/SRR576933.sra -P ${CHIP_DIR}/data/GSM1010219
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX189/SRX189774/SRR576934/SRR576934.sra -P ${CHIP_DIR}/data/GSM1010220


### Download reference genome & annotations

wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/fasta/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna.genome.fa.gz -P ${CHIP_DIR}/genome
wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gff3/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gff3.gz -P ${CHIP_DIR}/genome
wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gtf/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gtf.gz -P ${CHIP_DIR}/genome
gunzip ${CHIP_DIR}/genome/*.gz


## Setup RNA-seq analysis environment

export RNA_DIR=/data/analyses/RNA-seq_PE_GSE41190

mkdir -p ${RNA_DIR}
mkdir -p ${RNA_DIR}/data 
mkdir -p ${RNA_DIR}/genome

cd ${RNA_DIR}
<!--git clone https://github.com/rioualen/gene-regulation.git -->
<!--TODO: replace with tar.gz download of gene-regulation-v3.0-->
ln -s ~/Desktop/workspace/gene-regulation


### Download data

mkdir -p ${RNA_DIR}/data/GSM1010244 ${RNA_DIR}/data/GSM1010245 ${RNA_DIR}/data/GSM1010246 ${RNA_DIR}/data/GSM1010247
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059774/SRR191809/SRR191809.sra -P ${RNA_DIR}/data/GSM1010244
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059768/SRR191805/SRR191805.sra -P ${RNA_DIR}/data/GSM1010245
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059767/SRR191812/SRR191812.sra -P ${RNA_DIR}/data/GSM1010246
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059767/SRR191812/SRR191812.sra -P ${RNA_DIR}/data/GSM1010247

<!--Note: sample GSM1010247 is oddly formatted, so for we use a trick to run the workflow, by duplicating GSM1010245-->
<!--wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX116%2FSRX116381/SRR400301/SRR400301.sra -P ${RNA_DIR}/data/GSM1010247-->


### Download reference genome & annotations

wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/fasta/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna.genome.fa.gz -P ${RNA_DIR}/genome
wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gff3/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gff3.gz -P ${RNA_DIR}/genome
wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gtf/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gtf.gz -P ${RNA_DIR}/genome
gunzip ${RNA_DIR}/genome/*.gz


## Run the combined workflow

export COMBINED_DIR=/data/analyses/Combined_ChIP-seq_RNA-seq
mkdir -p ${COMBINED_DIR}
ln -s ~/Desktop/workspace/gene-regulation ${COMBINED_DIR}/gene-regulation
snakemake -s ${COMBINED_DIR}/gene-regulation/scripts/snakefiles/workflows/combined_ChIP-seq_RNA-seq.py -j 10 -p --nolock --configfile ${COMBINED_DIR}/gene-regulation/examples/Combined_ChIP-seq_RNA-seq/config.yml

