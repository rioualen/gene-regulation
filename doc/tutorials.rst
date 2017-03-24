
Tutorials
================================================================

Run with test data
----------------------------------------------------------------

Worklfow quality
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



Workflow ChIP-seq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


ChIP-seq study case in *S. cerevisiae*
****************************************************************




**Reference**

Preti M, Ribeyre C, Pascali C, Bosio MC et al. The telomere-binding
protein Tbf1 demarcates snoRNA gene promoters in Saccharomyces
cerevisiae. Mol Cell 2010 May 28;38(4):614-20. PMID: 20513435

**Access link**

http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20870

**ID**

GSE20870


Setup analysis environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    ANALYSIS_DIR=$HOME/ChIP-seq_GSE20870
    mkdir -p ${ANALYSIS_DIR}
    cd ${ANALYSIS_DIR}
    ln -s ${HOME}/gene-regulation-4.0 gene-regulation


Download data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    mkdir -p ${ANALYSIS_DIR}/data/GSM521934 ${ANALYSIS_DIR}/data/GSM521935
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021358/SRR051929/SRR051929.sra -P ${ANALYSIS_DIR}/data/GSM521934
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX021%2FSRX021359/SRR051930/SRR051930.sra -P ${ANALYSIS_DIR}/data/GSM521935

Download reference genome & annotations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    mkdir -p ${ANALYSIS_DIR}/genome
    wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.30.dna.genome.fa.gz -P ${ANALYSIS_DIR}/genome
    wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.30.gff3.gz -P ${ANALYSIS_DIR}/genome
    wget -nc ftp://ftp.ensemblgenomes.org/pub/fungi/release-30/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.30.gtf.gz -P ${ANALYSIS_DIR}/genome
    gunzip ${ANALYSIS_DIR}/genome/*.gz

Execute workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    snakemake -s ${ANALYSIS_DIR}/gene-regulation/scripts/snakefiles/workflows/ChIP-seq.wf -p -j 5 --configfile ${ANALYSIS_DIR}/gene-regulation/examples/ChIP-seq_SE_GSE20870/config.yml


.. figure:: rulegraph.png
   :alt: 


Genome-scale analysis of escherichia coli FNR
****************************************************************


**Reference**

Myers KS, Yan H, Ong IM, Chung D et al. Genome-scale analysis of
escherichia coli FNR reveals complex features of transcription factor
binding. PLoS Genet 2013 Jun;9(6):e1003565. PMID:
`23818864 <http://www.ncbi.nlm.nih.gov/pubmed/23818864>`__

**GEO series**

`GSE41187 <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41187>`__

File organisation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Short read files (.sra) are stored in sample-specific folders named
data/source/[GSM\_ID]. In this dataset, there is a single run per
sample, and thus a single sra file per GSM ID.

Setup analysis environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    ln -s ${HOME}/gene-regulation-4.0 gene-regulation
    ANALYSIS_DIR=$HOME/ChIP-seq_GSE41187
    mkdir -p ${ANALYSIS_DIR}
    cd ${ANALYSIS_DIR}
    CONFIG=gene-regulation/examples/ChIP-seq_SE_GSE41187/config.yml

Download data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    mkdir -p ${ANALYSIS_DIR}/data/GSM1010224 ${ANALYSIS_DIR}/data/GSM1010219 ${ANALYSIS_DIR}/data/GSM1010220
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX189%2FSRX189778/SRR576938/SRR576938.sra -P ${ANALYSIS_DIR}/data/GSM1010224
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX189%2FSRX189773/SRR576933/SRR576933.sra -P ${ANALYSIS_DIR}/data/GSM1010219
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX189/SRX189774/SRR576934/SRR576934.sra -P ${ANALYSIS_DIR}/data/GSM1010220

Download reference genome & annotations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/fasta/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna.genome.fa.gz -P ${ANALYSIS_DIR}/genome
    wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gff3/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gff3.gz -P ${ANALYSIS_DIR}/genome
    wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gtf/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gtf.gz -P ${ANALYSIS_DIR}/genome
    gunzip ${ANALYSIS_DIR}/genome/*.gz

Execute workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    snakemake -s ${ANALYSIS_DIR}/gene-regulation/scripts/snakefiles/workflows/import_to_fastq.wf -p --configfile ${ANALYSIS_DIR}/${CONFIG}
    snakemake -s ${ANALYSIS_DIR}/gene-regulation/scripts/snakefiles/workflows/quality_control.wf -p --configfile ${ANALYSIS_DIR}/${CONFIG}
    snakemake -s ${ANALYSIS_DIR}/gene-regulation/scripts/snakefiles/workflows/ChIP-seq.wf -p --configfile ${ANALYSIS_DIR}/${CONFIG}

.. figure:: ../img/rulegraph.png
   :alt: 

Workflow RNA-seq DEG
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

41190
****************************************************************

Study case for the development of a Bacterial RNA-seq workflow
==============================================================

Data source
-----------

**Reference**

Myers KS, Yan H, Ong IM, Chung D et al. Genome-scale analysis of
escherichia coli FNR reveals complex features of transcription factor
binding. PLoS Genet 2013 Jun;9(6):e1003565. PMID:
`23818864 <http://www.ncbi.nlm.nih.gov/pubmed/23818864>`__

**GEO series**

`GSE41190 <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41190>`__

File organisation
-----------------

Short read files (.sra) are stored in sample-specific folders named
data/source/[GSM\_ID]. In this dataset, there is a single run per
sample, and thus a single sra file per GSM ID.

Setup analysis environment
--------------------------

::

    ANALYSIS_DIR=/your/preferred/dir/RNA-seq_PE_GSE41190

    mkdir -p ${ANALYSIS_DIR}

    cd ${ANALYSIS_DIR}
    git clone https://github.com/rioualen/gene-regulation.git 
    <!--TODO: replace with tar.gz download of gene-regulation-v3.0-->

Download data
-------------

<!-- Old data (still downloadable)

::

    mkdir -p ${ANALYSIS_DIR}/data/GSM1010244 ${ANALYSIS_DIR}/data/GSM1010245 ${ANALYSIS_DIR}/data/GSM1010246 ${ANALYSIS_DIR}/data/GSM1010247
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059774/SRR191809/SRR191809.sra -P ${ANALYSIS_DIR}/data/GSM1010244
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059768/SRR191805/SRR191805.sra -P ${ANALYSIS_DIR}/data/GSM1010245
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059767/SRR191812/SRR191812.sra -P ${ANALYSIS_DIR}/data/GSM1010246
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX059%2FSRX059767/SRR191812/SRR191812.sra -P ${ANALYSIS_DIR}/data/GSM1010247

-->

.. raw:: html

   <!--Note: sample GSM1010247 is oddly formatted, so for we use a trick to run the workflow, by duplicating GSM1010245-->

<!--wget --no-clobber
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX116%2FSRX116381/SRR400301/SRR400301.sra
-P ${ANALYSIS\_DIR}/data/GSM1010247-->

::

    mkdir -p ${ANALYSIS_DIR}/data/GSM1010244 ${ANALYSIS_DIR}/data/GSM1010245 ${ANALYSIS_DIR}/data/GSM1010246 ${ANALYSIS_DIR}/data/GSM1010247
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX264/SRX2641374/SRR5344681/SRR5344681.sra -P ${ANALYSIS_DIR}/data/GSM1010244
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX264/SRX2641375/SRR5344682/SRR5344682.sra -P ${ANALYSIS_DIR}/data/GSM1010245
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX264/SRX2641376/SRR5344683/SRR5344683.sra -P ${ANALYSIS_DIR}/data/GSM1010246
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX264/SRX2641377/SRR5344684/SRR5344684.sra -P ${ANALYSIS_DIR}/data/GSM1010247

Download reference genome & annotations
---------------------------------------

::

    wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/fasta/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna.genome.fa.gz -P ${ANALYSIS_DIR}/genome
    wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gff3/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gff3.gz -P ${ANALYSIS_DIR}/genome
    wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gtf/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gtf.gz -P ${ANALYSIS_DIR}/genome
    gunzip ${ANALYSIS_DIR}/genome/*.gz

Execute workflow
----------------

::

    snakemake -s ${ANALYSIS_DIR}/gene-regulation/scripts/snakefiles/workflows/RNA-seq_workflow_PE.py -p -j 5 --configfile ${ANALYSIS_DIR}/gene-regulation/examples/RNA-seq_PE_GSE41190/config.yml

Rulegraph
---------

.. figure:: rulegraph.png
   :alt: 

71562
*****

Study case for the development of a Bacterial RNA-seq workflow
==============================================================

We show here how to run a template workflow for the detection of
differentially expressed genes in the bacteria *Escherichia coli K12*
submittted to a change from anaerobic to aerobic conditions.

Data
----

**Reference**

GEO series:
`GSE71562 <http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71562>`__

File organisation
-----------------

Short read files (.sra) are stored in sample-specific folders named
data/source/[GSM\_ID]. In this dataset, there is a single run per
sample, and thus a single sra file per GSM ID.

Setup analysis environment
--------------------------

We will create a specific folder for the analysis of this series, and
include in the same folder a clone of the *gene-regulation* repository,
which contains the *snakemake* librairies, workflows and configuration
files required to run the whole analytic workflow. The inclusion of hte
*gene-regulation* package within the analysis fodler ensures consistency
between the code and the results.

::

    ANALYSIS_DIR=${HOME}/analyses/RNA-seq_SE_GSE71562

    mkdir -p ${ANALYSIS_DIR}
    cd ${ANALYSIS_DIR}
    git clone https://github.com/rioualen/gene-regulation.git
    <!--TODO: replace with tar.gz download of gene-regulation-v3.0-->

Download data
-------------

We download the short read files from the SRA dtabase. Beware, the sra
files require 3Gb.

::

    mkdir -p ${ANALYSIS_DIR}/data/GSM1838496 ${ANALYSIS_DIR}/data/GSM1838502 ${ANALYSIS_DIR}/data/GSM1838508 
    mkdir -p ${ANALYSIS_DIR}/data/GSM1838499 ${ANALYSIS_DIR}/data/GSM1838505 ${ANALYSIS_DIR}/data/GSM1838511
    mkdir -p ${ANALYSIS_DIR}/data/GSM1838501 ${ANALYSIS_DIR}/data/GSM1838507 ${ANALYSIS_DIR}/data/GSM1838513

    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX112/SRX1125282/SRR2135663/SRR2135663.sra -P ${ANALYSIS_DIR}/data/GSM1838496
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX112/SRX1125288/SRR2135669/SRR2135669.sra -P ${ANALYSIS_DIR}/data/GSM1838502
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX112/SRX1125294/SRR2135675/SRR2135675.sra -P ${ANALYSIS_DIR}/data/GSM1838508

    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX112/SRX1125285/SRR2135666/SRR2135666.sra -P ${ANALYSIS_DIR}/data/GSM1838499
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX112/SRX1125291/SRR2135672/SRR2135672.sra -P ${ANALYSIS_DIR}/data/GSM1838505
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX112/SRX1125297/SRR2135678/SRR2135678.sra -P ${ANALYSIS_DIR}/data/GSM1838511

    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX112/SRX1125287/SRR2135668/SRR2135668.sra -P ${ANALYSIS_DIR}/data/GSM1838501
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX112/SRX1125293/SRR2135674/SRR2135674.sra -P ${ANALYSIS_DIR}/data/GSM1838507
    wget --no-clobber ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX112/SRX1125299/SRR2135680/SRR2135680.sra -P ${ANALYSIS_DIR}/data/GSM1838513

Download reference genome & annotations
---------------------------------------

We download the reference genome sequence and annotations.

.. raw:: html

   <!-- TO DO: integrate this in the rule, and the full URLs should be in the config file-->

::

    wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/fasta/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.dna.genome.fa.gz -P ${ANALYSIS_DIR}/genome
    wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gff3/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gff3.gz -P ${ANALYSIS_DIR}/genome
    wget -nc ftp://ftp.ensemblgenomes.org/pub/release-21/bacteria/gtf/bacteria_22_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.1.21.gtf.gz -P ${ANALYSIS_DIR}/genome
    gunzip ${ANALYSIS_DIR}/genome/*.gz

Execute workflow
----------------

::

    cd ${ANALYSIS_DIR}
    snakemake -s gene-regulation/scripts/snakefiles/workflows/RNA-seq_workflow_SE.wf -p -j 5 --configfile ${ANALYSIS_DIR}/gene-regulation/examples/RNA-seq_SE_GSE71562/config.yml

Rulegraph
---------

.. figure:: rulegraph.png
   :alt: 
...


Workflow RNA-seq infer transcripts from bam
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Workflow orthologs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With no ref genome, several ref genomes...

Integrated workflow ChIP-RNA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

todo



Run with your own data
----------------------------------------------------------------

File organization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ref genome(s)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Gene-reg library
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

workflow.wf or custom
****************************************************************

samples.tab
****************************************************************

design.tab
****************************************************************


config.yml
****************************************************************


