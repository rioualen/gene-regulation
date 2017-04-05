Study case combined ChIP-seq & RNA-seq workflow on E. coli data
===============================================================

Data source
-----------

Reference: Myers KS, Yan H, Ong IM, Chung D et al. Genome-scale analysis
of escherichia coli FNR reveals complex features of transcription factor
binding. PLoS Genet 2013 Jun;9(6):e1003565. PMID:
`23818864 <http://www.ncbi.nlm.nih.gov/pubmed/23818864>`__

GEO series:
`GSE41187 <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41187>`__
`GSE41190 <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41190>`__

GEO Superseries:
`GSE41195 <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41195>`__

Setup ChIP-seq analysis environment
-----------------------------------

export CHIP\_DIR=/data/analyses/ChIP-seq\_SE\_GSE41187

mkdir -p ${CHIP\_DIR} mkdir -p ${CHIP\_DIR}/data mkdir -p
${CHIP\_DIR}/genome

cd ${CHIP\_DIR} ln -s ~/Desktop/workspace/gene-regulation

Download data
~~~~~~~~~~~~~

mkdir -p ${CHIP\_DIR}/data/GSM1010224 ${CHIP\_DIR}/data/GSM1010219
${CHIP\_DIR}/data/GSM1010220 wget --no-clobber
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX189%2FSRX189778/SRR576938/SRR576938.sra
-P ${CHIP\_DIR}/data/GSM1010224 wget --no-clobber
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX189%2FSRX189773/SRR576933/SRR576933.sra
-P ${CHIP\_DIR}/data/GSM1010219 wget --no-clobber
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX189/SRX189774/SRR576934/SRR576934.sra
-P ${CHIP\_DIR}/data/GSM1010220

Download reference genome & annotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

wget -nc
ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/fasta/bacteria\_22\_collection/escherichia\_coli\_str\_k\_12\_substr\_mg1655/dna/Escherichia\_coli\_str\_k\_12\_substr\_mg1655.GCA\_000005845.1.21.dna.genome.fa.gz
-P ${CHIP\_DIR}/genome wget -nc
ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gff3/bacteria\_22\_collection/escherichia\_coli\_str\_k\_12\_substr\_mg1655/Escherichia\_coli\_str\_k\_12\_substr\_mg1655.GCA\_000005845.1.21.gff3.gz
-P ${CHIP\_DIR}/genome wget -nc
ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gtf/bacteria\_22\_collection/escherichia\_coli\_str\_k\_12\_substr\_mg1655/Escherichia\_coli\_str\_k\_12\_substr\_mg1655.GCA\_000005845.1.21.gtf.gz
-P ${CHIP\_DIR}/genome gunzip ${CHIP\_DIR}/genome/\*.gz

Setup RNA-seq analysis environment
----------------------------------

export RNA\_DIR=/data/analyses/RNA-seq\_PE\_GSE41190

mkdir -p ${RNA\_DIR} mkdir -p ${RNA\_DIR}/data mkdir -p
${RNA\_DIR}/genome

cd ${RNA\_DIR} ln -s ~/Desktop/workspace/gene-regulation

Download data
~~~~~~~~~~~~~

mkdir -p ${RNA\_DIR}/data/GSM1010244 ${RNA\_DIR}/data/GSM1010245
${RNA\_DIR}/data/GSM1010246 ${RNA\_DIR}/data/GSM1010247 wget
--no-clobber
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059774/SRR191809/SRR191809.sra
-P ${RNA\_DIR}/data/GSM1010244 wget --no-clobber
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059768/SRR191805/SRR191805.sra
-P ${RNA\_DIR}/data/GSM1010245 wget --no-clobber
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059767/SRR191812/SRR191812.sra
-P ${RNA\_DIR}/data/GSM1010246 wget --no-clobber
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059767/SRR191812/SRR191812.sra
-P ${RNA\_DIR}/data/GSM1010247

.. raw:: html

   <!--Note: sample GSM1010247 is oddly formatted, so for we use a trick to run the workflow, by duplicating GSM1010245-->

<!--wget --no-clobber
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX116%2FSRX116381/SRR400301/SRR400301.sra
-P ${RNA\_DIR}/data/GSM1010247-->

Download reference genome & annotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

wget -nc
ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/fasta/bacteria\_22\_collection/escherichia\_coli\_str\_k\_12\_substr\_mg1655/dna/Escherichia\_coli\_str\_k\_12\_substr\_mg1655.GCA\_000005845.1.21.dna.genome.fa.gz
-P ${RNA\_DIR}/genome wget -nc
ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gff3/bacteria\_22\_collection/escherichia\_coli\_str\_k\_12\_substr\_mg1655/Escherichia\_coli\_str\_k\_12\_substr\_mg1655.GCA\_000005845.1.21.gff3.gz
-P ${RNA\_DIR}/genome wget -nc
ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gtf/bacteria\_22\_collection/escherichia\_coli\_str\_k\_12\_substr\_mg1655/Escherichia\_coli\_str\_k\_12\_substr\_mg1655.GCA\_000005845.1.21.gtf.gz
-P ${RNA\_DIR}/genome gunzip ${RNA\_DIR}/genome/\*.gz

Run the combined workflow
-------------------------

export COMBINED\_DIR=/data/analyses/Combined\_ChIP-seq\_RNA-seq mkdir -p
${COMBINED\_DIR} ln -s ~/Desktop/workspace/gene-regulation
${COMBINED\_DIR}/gene-regulation snakemake -s
${COMBINED\_DIR}/gene-regulation/scripts/snakefiles/workflows/combined\_ChIP-seq\_RNA-seq.py
-j 10 -p --nolock --configfile
${COMBINED\_DIR}/gene-regulation/examples/Combined\_ChIP-seq\_RNA-seq/config.yml
