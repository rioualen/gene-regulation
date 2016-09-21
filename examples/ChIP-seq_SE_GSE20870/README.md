## Study case

The telomere-binding protein Tbf1 demarcates snoRNA gene promoters in Saccharomyces cerevisiae

Reference
    Preti M, Ribeyre C, Pascali C, Bosio MC et al. 
    The telomere-binding protein Tbf1 demarcates snoRNA gene promoters in Saccharomyces cerevisiae. 
    Mol Cell 2010 May 28;38(4):614-20. PMID: 20513435

Access link
    http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20870

## Setup analysis environment

export ANALYSIS_DIR=/data/analyses/ChIP-seq_SE_GSE20870

<!--mkdir -p ${ANALYSIS_DIR}-->
<!--mkdir -p ${ANALYSIS_DIR}/data -->
<!--mkdir -p ${ANALYSIS_DIR}/genome-->

cd ${ANALYSIS_DIR}
<!--git clone https://github.com/rioualen/gene-regulation.git-->
ln -s ~/Desktop/workspace/gene-regulation gene-regulation

## Download data

<!--mkdir -p ${ANALYSIS_DIR}/data/GSM521934 ${ANALYSIS_DIR}/data/GSM521935-->
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021358/SRR051929/SRR051929.sra -P ${ANALYSIS_DIR}/data/GSM521934
wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021359/SRR051930/SRR051930.sra -P ${ANALYSIS_DIR}/data/GSM521935


## Download reference genome & annotations

wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.30.dna.genome.fa.gz -P ${ANALYSIS_DIR}/genome
wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.30.gff3.gz -P ${ANALYSIS_DIR}/genome
wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.30.gtf.gz -P ${ANALYSIS_DIR}/genome
gunzip ${ANALYSIS_DIR}/genome/*.gz


## Execute workflow

snakemake -s ${ANALYSIS_DIR}/gene-regulation/scripts/snakefiles/workflows/ChIP-seq_workflow_SE.py -p -j 5 --configfile ${ANALYSIS_DIR}/gene-regulation/examples/ChIP-seq_SE_GSE20870/config.yml
