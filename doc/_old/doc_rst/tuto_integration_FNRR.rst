

Genome-scale analysis of *Escherichia coli* FNR
----------------------------------------------------------------

Presentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Reference**

Myers KS, Yan H, Ong IM, Chung D et al. Genome-scale analysis of
Escherichia coli FNR reveals complex features of transcription factor
binding. PLoS Genet 2013 Jun;9(6):e1003565. PMID:
`23818864 <http://www.ncbi.nlm.nih.gov/pubmed/23818864>`__

**GEO series**

- ChIP-seq: `GSE41187 <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41187>`__
- RNA-seq: `GSE41190 <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41190>`__

Download reference genome & annotations
****************************************************************

::

    mkdir ${GENOME_DIR}/Ecoli-K12
    wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/fasta/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna.genome.fa.gz -P ${GENOME_DIR}/Ecoli-K12
    wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gff3/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gff3.gz -P ${GENOME_DIR}/Ecoli-K12
    wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gtf/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gtf.gz -P ${GENOME_DIR}/Ecoli-K12
    gunzip ${GENOME_DIR}/Ecoli-K12/*.gz

Setup analysis environment
****************************************************************

::

    ANALYSIS_DIR=${HOME}/Integrated_analysis


Workflow 'ChIP-seq'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Setup analysis environment
****************************************************************

::

    ANALYSIS_DIR_CHIP=${ANALYSIS_DIR}/ChIP-seq_GSE41187
    mkdir -p ${ANALYSIS_DIR_CHIP} 
    ln -s ${GENE_REG_PATH} ${ANALYSIS_DIR_CHIP}/gene-regulation             # ${GENE_REG_PATH} should have been defined beforehand
    ln -s ${GENOME_DIR} ${ANALYSIS_DIR_CHIP}/genome                         # ${GENOME_DIR} should have been defined beforehand
    CONFIG_CHIP=${ANALYSIS_DIR_CHIP}/gene-regulation/examples/ChIP-seq_SE_GSE41187/config.yml

Download ChIP-seq data 
****************************************************************

::

    mkdir -p ${ANALYSIS_DIR_CHIP}/data/GSM1010224 ${ANALYSIS_DIR_CHIP}/data/GSM1010219 ${ANALYSIS_DIR_CHIP}/data/GSM1010220
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX189%2FSRX189778/SRR576938/SRR576938.sra -P ${ANALYSIS_DIR_CHIP}/data/GSM1010224
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX189%2FSRX189773/SRR576933/SRR576933.sra -P ${ANALYSIS_DIR_CHIP}/data/GSM1010219
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX189/SRX189774/SRR576934/SRR576934.sra -P ${ANALYSIS_DIR_CHIP}/data/GSM1010220

Workflow execution
****************************************************************

Your directory should now look like this:


.. figure:: ../img/arbo_tuto_FNR_ChIP.png
   :alt: 

And you should be able to execute it like this:

::

    snakemake -s ${ANALYSIS_DIR_CHIP}/gene-regulation/scripts/snakefiles/workflows/import_from_sra.wf -p --configfile ${CONFIG_CHIP}
    snakemake -s ${ANALYSIS_DIR_CHIP}/gene-regulation/scripts/snakefiles/workflows/quality_control.wf -p --configfile ${CONFIG_CHIP}
    snakemake -s ${ANALYSIS_DIR_CHIP}/gene-regulation/scripts/snakefiles/workflows/ChIP-seq.wf -p --configfile ${CONFIG_CHIP}

Workflow 'RNA-seq' DEG
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Setup analysis environment
****************************************************************

::

    ANALYSIS_DIR_RNA=${ANALYSIS_DIR}/RNA-seq_GSE41190
    mkdir ${ANALYSIS_DIR_RNA}
    ln -s ${GENE_REG_PATH} ${ANALYSIS_DIR_RNA}/gene-regulation
    ln -s ${GENOME_DIR} ${ANALYSIS_DIR_RNA}/genome
    CONFIG_RNA=${ANALYSIS_DIR_RNA}/gene-regulation/examples/RNA-seq_PE_GSE41190/config.yml

Download RNA-seq data
****************************************************************

::

    mkdir -p ${ANALYSIS_DIR_RNA}/data/GSM1010244 ${ANALYSIS_DIR_RNA}/data/GSM1010245 ${ANALYSIS_DIR_RNA}/data/GSM1010246 ${ANALYSIS_DIR_RNA}/data/GSM1010247
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX264/SRX2641374/SRR5344681/SRR5344681.sra -P ${ANALYSIS_DIR_RNA}/data/GSM1010244
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX264/SRX2641375/SRR5344682/SRR5344682.sra -P ${ANALYSIS_DIR_RNA}/data/GSM1010245
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX264/SRX2641376/SRR5344683/SRR5344683.sra -P ${ANALYSIS_DIR_RNA}/data/GSM1010246
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX264/SRX2641377/SRR5344684/SRR5344684.sra -P ${ANALYSIS_DIR_RNA}/data/GSM1010247

Workflow execution
****************************************************************

::

    snakemake -s ${ANALYSIS_DIR_RNA}/gene-regulation/scripts/snakefiles/workflows/import_from_sra.wf -p --configfile ${CONFIG_RNA}
    snakemake -s ${ANALYSIS_DIR_RNA}/gene-regulation/scripts/snakefiles/workflows/quality_control.wf -p --configfile ${CONFIG_RNA}
    snakemake -s ${ANALYSIS_DIR_RNA}/gene-regulation/scripts/snakefiles/workflows/RNA-seq_workflow_PE.py -p --configfile ${CONFIG_RNA}


.. figure:: ../examples/RNA-seq_PE_GSE41190/rulegraph.png
   :alt: 


Workflow 'integrated_ChIP_RNA'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*todo*


