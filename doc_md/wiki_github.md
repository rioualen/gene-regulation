# Versioning, code sharing

## GitHub 

[Link](http://github.com)

A few more links:

**TODO**

## SourceForge 

[Link](http://sourceforge.net)

## BitBucket
 
[Link](http://bitbucket.org/)

## SourcesSup Renater 

[Link](http://sourcesup.renater.fr)


## How to...

### Find answers to probably-already-asked questions?

Check out these Q & A sites:

* [SeqAnswers](http://seqanswers.com/)
* [Biostars](https://www.biostars.org/)
* [Biostars Galaxy](https://biostar.usegalaxy.org/)

For questions related to computing problems:

* [Stack Overflow](http://stackoverflow.com/)
* [Ask Ubuntu](http://askubuntu.com/)
* [Super User](http://superuser.com/)

### Get the best out of quality control?

* [QC Fail Sequencing](https://sequencing.qcfail.com/)

* [FastQC results interpretation](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)

### Download and install the SPP peak-caller?

I've struggled with this one so I'll leave a few notes here: 

* Download page can be found [here](http://compbio.med.harvard.edu/Supplements/ChIP-seq/), better chose version `1.11`.
* SPP requires the Bioconductor library [Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html) to be installed beforehand.
* Unix packages `gcc` and `libboost` (or equivalents) must be installed. 
* You can find a few more notes [here](http://seqanswers.com/forums/archive/index.php/t-22653.html).
* Good luck!

Here's the procedure I use on Ubuntu 14.04, in this very order:

In unix shell:
```
# unix libraries
apt-get update
apt-get -y install r-base
apt-get -y install libboost-dev zlibc zlib1g-dev
```

In R shell:
```
# Rsamtools
source("http://bioconductor.org/biocLite.R")
biocLite("Rsamtools")
```

In unix shell:
```
# spp
wget http://compbio.med.harvard.edu/Supplements/ChIP-seq/spp_1.11.tar.gz
sudo R CMD INSTALL spp_1.11.tar.gz
```

## Find the exact size of genomes of model organisms?

* [Link](http://users.rcn.com/jkimball.ma.ultranet/BiologyPages/G/GenomeSizes.html)


# File formats

*NB* Should be merged with glossary section?

A list of formats maintained by the [UCSC](http://genome.ucsc.edu/FAQ/FAQformat.html).

 * **fastq:** raw sequences + quality ([more](http://maq.sourceforge.net/fastq.shtml)).
 * **sam:** aligned reads ([more](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full.pdf)).
 * **bam:** compressed sam ([more](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full.pdf)).
 * **cigar:** alignment ([more]()).
 * **bed:** genomic coordinates ([more](http://genome.ucsc.edu/FAQ/FAQformat.html#format1)).
 * **wig:** coverage / density of some signal along genome ([more](http://genome.ucsc.edu/goldenPath/help/wiggle.html)).
 * **bedgraph:** positions + density values ([more](http://genome.ucsc.edu/goldenPath/help/bedgraph.html)).
 * **gff:** genome feature file - annotations ([more](http://www.sanger.ac.uk/resources/software/gff/spec.html)).
 * **gtf:** variant of GFF, with two fields for annotation ([more](http://www.ensembl.org/info/website/upload/gff.html)).
 * **gft2:** Gene annotation ([more](http://mblab.wustl.edu/GTF22.html)).
 * **fastq:** raw reads + quality scores ([more](http://maq.sourceforge.net/fastq.shtml)).
 * **pileup:** base-pair information at each chromosomal position ([more](http://samtools.sourceforge.net/pileup.shtml)).
 * **vcf:** variant call format ([more](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41)).


# French Bioinformatics Institute (IFB)

The [IFB cloud](http://cloud.france-bioinformatique.fr) proposes many services such as: 

  * Galaxy
  * RSAT
  * R

[Documentation](http://www.france-bioinformatique.fr/?q=fr/core/cellule-infrastructure/documentation-cloud)

# Glossary


We define hereafter a series of abbreviations, terms and concepts which appear recurrently in the litterature about NGS analysis. This document aims at providing a support for the interpretation of analysis reports. 

Other resources:

* [link](https://github.com/fidelram/deepTools/wiki/Glossary)

## A

## B

  * **bam (file format):** compressed sam ([more](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full.pdf)).
  * **bed (file format):** genomic coordinates ([more](http://genome.ucsc.edu/FAQ/FAQformat.html#format1)).
  * **bedgraph (file format):** positions + density values ([more](http://genome.ucsc.edu/goldenPath/help/bedgraph.html)).
  * **bin:** 
  * **Bonferroni correction:** used in **multiple testing**. Consists in adapting the alpha threshold rather than correcting the **p-value**.
  
## C

  * **ChIP-exo:** 
  * **ChIP-seq:** 
  * **cigar (file format):** alignment.
  * **Cloud:** 
  * **Copy number variation:** 
  * **Core:** 
  
## D

  * **DEG/Differentially Expressed Gene:** 
  
## E

  * **e-value (E):** indicates the number of false positives expected by chance, for a given threshold of **p-value**. It is a number that can exceed 1, it is thus not a probability, and thus, not a p-value.

E = \<FP\> = P \. m

Where **m** is the number of tests (e.g. genes), FP the number of false positives, the notation \< \> denotes the random expectation, and P is the nominal p-value of the considered gene.

Note that the e-value is a positive number ranging from 0 to m (number of tests). It is thus not a p-value, since probabilities are by definition comprized between 0 and 1. 
  
## F

  * **Family-wise error rate (FWER):**  indicates the probability to observe at least one false positive among the multiple tests. 

FWER = P(FP >= 1)

  * **fastq (file format):** raw sequences + quality ([more](http://maq.sourceforge.net/fastq.shtml)).
  * **False discovery rate (FDR):** indicates the expected proportion of false positives *among the cases declared positive*. For example, if a differential analysis reports 200 differentially expressed genes with an FDR threshold of 0.05, we should expect to have 0.05 x 200=10 false positive among them. 

  
## G

  * **genome (file format):**
  * **genomic input:**
  * **gff (file format):** genome feature file - annotations ([more](http://www.sanger.ac.uk/resources/software/gff/spec.html)). See also `gtf`.
  * **gtf (file format):** variant of GFF, with two fields for annotation ([more](http://www.ensembl.org/info/website/upload/gff.html)).
  * **gft2 (file format):** Gene annotation ([more](http://mblab.wustl.edu/GTF22.html)).
  * **GSM:** Gene Expression Omnibus Sample identifier.
  * **GSE:** Gene Expression Omnibus Series identifier (a collection of samples related to the same publication or thematics).

## H

## I

  * **input:** Pour le peak-calling, le mot "input" est utilisé dans un sens tout à fait particulier, pour désigner un jeu de séquences servant à estimer les densités de reads attendues au hasard en fonction de la position génomique. Les méthodes typiques sont l'input génomique (actuellement le plus généralement utilisé) et le mock.
  
## J

## K

## L

  * **Library:** Terme utilisé de façon parfois ambiguë selon le contxte. 
Les biologistes se réfèrent à une librairie d'ADN pour désigner ... (à définir). Les bioinformaticiens parlent de librairie de séquences pour désigner l'nsemble des fragments de lectures provenant du séquençage d'un même échantillon. Les informaticiens appellent ""library"" (bibliothèque, librairies ?) des modules de code regroupant une série de fonctions et procédures.

## M

  * **m:** number of tests in a multiple-testing schema (e.g. number of genes in differential expression analysis).
  * **Mapped read:** 
  * **Mapping:** Identifying genomic positions for the raw reads of a sequence library.
  * **mock:** type of control for the peak-calling in ChIP-seq. It is an input obtained by using a non-specific antibody (eg. anti-GFP) for the immunoprecipitation. *afin d'estimer le taux de séquençage aspécifique pour chaque région génomique. L'intérêt du mock est qu'il constitue un contrôle réalisé dans les mêmes conditions que le ChIP-seq spécifique. La faiblesse est que les tailles de librairries sont parfois tellement faibles que l'estimation du backgroun est très peu robuste.
  * **motif:** 
  * **Multiple testing:** the multiple testing problem arises from the application of a given statistical test to a large number of cases. For example, in differential expression analysis, each gene/transcript is submitted to a test of equality between two conditions. A single analysis thus typically involves  several tens of thousands tests. 
The general problem of multiple testing is that the risk of false positive indicated by the nominal **p-value** will be challenged for each element. Various types of corrections for multiple testing have been defined (**Bonferroni**, **e-value**, **FWER**, **FDR**).


## N

  * **Negative control:** 
  * **Next Generation Sequencing:**

## O

## P

  * **p-value (P):** the **nominal p-value** is the p-value attached to one particular element in a series of multiple tests. For example, in differential analysis, one nominal p-value is computed for each gene. This p-value indicates the risk to obtain an effect at least as important as our observation *under the null hypothesis*, i.e. in the absence of regulation. 
  * **padj (abbr.):** adjusted p-value. Statistics derived from the nominal **p-value** in order to correct for the effects of **multiple testing** (see **Bonferroni correction**, **e-value**). 

The most usual correction is the FDR, which can be estimated in various ways. 

  * **Paired end:** 
  * **Peak:** 
  * **Peak-calling:** 
  * **pileup (file format):** base-pair information at each chromosomal position ([more](http://samtools.sourceforge.net/pileup.shtml)).

## Q

  * **q-value:** 

## R

  * **RAM:** 
  * **Raw read:** non-aligned read.
  * **Read:** short sequence (typically 25-75bp) obtained by high-throughput sequencing.
  * **Region-calling:** 
  * **Replicate:** ... distinguer réplicat technique et réplicat biologique
  * **RNA-seq:** 

## S

  * **sam (file format):** aligned reads ([more](http://bioinformatics.oxfordjournals.org/content/25/16/2078.full.pdf)).
  * **Single end:**
  * **Single nucleotide polymorphism:** 
  * **SRA:** Sequence Read Archive (SRA). Database maintained by the [NCBI](www.ncbi.nlm.nih.gov/sra).
  * **SRX:** Short Read Experiment. See [documentation](www.ncbi.nlm.nih.gov/books/NBK56913/#search.the_entrez_sra_search_response_pa).
  * **SRR:** Short Read Run. See [documentation](www.ncbi.nlm.nih.gov/books/NBK56913/#search.each_srx_entry_in_the_entrez_sra).

## T

## U

## V
 
  * **vcf (file format):** variant call format ([more](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41)).
  * **Virtual machine:** 
  
## W

  * **wig (file format):** coverage / density of some signal along genome ([more](http://genome.ucsc.edu/goldenPath/help/wiggle.html)).

## X

## Y

## Z




## Abbreviations

| Abbreviation | Meaning |
|--------------|---------|
| NGS | Next Generation Sequencing |
| DEG | Differentially expressed gene(s) |
| padj | Adjusted p-value |
| FDR | False Discovery Rate |
| | |

## Symbols

| Symbol | Meaning |
|--------------|---------|
| $m$ | Number of tests in a multiple-testing schema (e.g. number of genes in differential expression analysis) |
| $P$ | p-value |
| $E$ | e-value |


* * * *

## Multiple testing corrections

### The problem with multiple tests

The **multiple testing** problem arises from the application of a given statistical test to a large number of cases.  For example, in differential expression analysis, each gene/transcript is submitted to a test of equality between two conditions. A single analysis thus typically involves  several tens of thousands tests. 



The general problem of **multiple testing** is that the risk of false positive indicated by the nominal p-value will be challenged for each element. 


### P-value and derived multiple testing corrections



### P-value (nominal p-value)

The **nominal p-value** is the p-value attached to one particular element in a series of multiple tests. For example, in differential analysis, one nominal p-value is computed for each gene. This p-value indicates the risk to obtain an effect at least as important as our observation *under the null hypothesis*, i.e. in the absence of regulation. 

### Bonferroni correction

### E-value

The **e-value** indicates the number of false positives expected by chance, for a given threshold of p-value. 

$E = <FP> = P \cdot m$

Where $m$ is the number of tests (e.g. genes), $FP$ the number of false positives, the notation $< >$ denotes the random expectation, and $P$ is the nominal p-value of the considered gene.

Note that the e-value is a positive number ranging from $0$ to $m$ (number of tests). It is thus not a p-value, since probabilities are by definition comprized between 0 and 1. 

### Family-wise error rate (FWER)

The Family-Wise Error Rate (**FWER**) indicates the probability to observe at least one false positive among the multiple tests. 

$FWER = P(FP >= 1)$

### False Discovery Rate (FDR)

The **False Discovery Rate** (**FDR**) indicates the expected proportion of false positives *among the cases declared positive*. For example, if a differential analysis reports 200 differentially expressed genes with an FDR threshold of 0.05, we should expect to have $0.05 \cdot 200=10$ false positive among them. 


### What is an adjusted p-value?

An **adjusted p-value** is a statistics derived from the nominal p-value in order to correct for the effects of multiple testing. 

Various types of corrections for multiple testing have been defined (Bonferoni, e-value, FWER, FDR). Note that some of these corrections are not actual "adjusted p-values". 

- the original Bonferoni correction consists in adapting the $\alpha$ threshold rather than correcting the p-value.
- the e-value is a number that can exceed 1, it is thus not a probability, and thus, not a p-value.

The most usual correction is the FDR, which can be estimated in various ways. 


* * * 

# NGS analysis tools

This list is far from exhaustive. You can find other lists:

* [Sequencing (OmicTools)](http://omictools.com/sequencing-c152-p1.html)
* [ChIP-seq (OmicTools)](http://omictools.com/chip-seq-c1215-p1.html)
* [ChIP-seq](https://github.com/crazyhottommy/ChIP-seq-analysis)

* Elixir's [Tools and Data Services Registry](https://bio.tools/)

## Quality assessment

### FastQC

[Link](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

The main functions of FastQC are:

* Import of data from BAM, SAM or FastQ files (any variant)
* Providing a quick overview to tell you in which areas there may be problems
* Summary graphs and tables to quickly assess your data
* Export of results to an HTML based permanent report
* Offline operation to allow automated generation of reports without running the interactive application


## Trimming

The quality of the reads generated by high-throughput sequencing technologies tend to decrease at their ends. Trimming consists in cutting out theses ends, and thus better the quality of reads before the mapping.

### Sickle

[Link](https://github.com/najoshi/sickle)

Sickle uses sliding windows computing sequencing quality along the reads. When the quality falls below a chose q-value threshold, the reads is cut. If the size of the remaining read is too short, it is completely removed. 
Sickle takes into account three different types of read quality: Illumina, Solexa, Sanger. 

### Cutadapt

## Alignment/mapping

*note* find that presentation explaining the difference between the 2 of em

*to translate*

Le but de l'alignement est de replacer les //reads// issus du séquençage à leur emplacement sur un génome de référence. Lorsque le //read// est suffisamment long, il peut généralement être //mappé// sur le génome avec une bonne certitude, en tolérant une certain quantité de //mismatches//, c'est-à-dire de nucléotides mal appariés. Néanmoins certaines séquences répétées du génome peuvent s'avérer plus difficiles à aligner. 
On désigne par l'expression "profondeur de séquençage" (ou //sequencing depth//) le nombre moyen de //reads// alignés par position sur le génome. Plus cette profondeur est importante, meilleure est la qualité de l'alignement, et plus les analyses ultérieures seront de qualité. 

### BWA

[Link](http://bio-bwa.sourceforge.net/)

[Manual](http://bio-bwa.sourceforge.net/bwa.shtml)

[Publication](http://www.ncbi.nlm.nih.gov/pubmed/19451168): Li H. and Durbin R. (2009). Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60.

### Bowtie

### Others 

A list on [Wikipedia](https://en.wikipedia.org/wiki/List_of_sequence_alignment_software)

## Peak-calling

**TODO** section à étoffer, voir protocole d'install snakemake.

### MACS14/MACS2
### SWEMBL
### HOMER (findPeaks)
### SPP
### bPeaks
### SICER

## Motif analysis

### Regulatory Sequence Analysis Tools (RSAT)

[Link](http://rsat.eu/)

*to translate*

Suite logicielle spécialisée pour l'analyse de motifs cis-régulateurs, développée par les équipes de Morgane Thomas-Chollier (ENS, Paris) et Jacques van Helden (TAGC, Marseille). Inclut des outils spécifiques pour l'analyse de données de ChIP-seq.

More: see the tutorials section in `resources.md`.

### MEME

[Link](http://meme.ebi.edu.au/meme/doc/meme-chip.html)

*to translate*

Suite logicielle spécialisée pour l'analyse de motifs cis-régulateurs, développée par l'équipe de Tim Bailey. Inclut des outils spécifiques pour l'analyse de données de ChIP-seq. 



## File conversion

### SamTools

[Link](http://samtools.sourceforge.net/)

### BamTools

### BedTools

[Link](http://bedtools.readthedocs.org/en/latest/)

### SRA Toolkit

[Documentation](http://www.ncbi.nlm.nih.gov/Traces/sra/?view=toolkit_doc)

Set of tools for the conversion of `*.sra` files (sequence read archive) into several formats. `fastq-dump` converts to `*.fastq` files.

* More info on the [SRA format](http://www.ncbi.nlm.nih.gov/Traces/sra/)
* fastq-dump [manual](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump)
* [Download](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software)
* [Install guide](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std)

## Miscellaneous

* [MICSA](http://bioinfo-out.curie.fr/software.html): peak-calling & motifs discovery ([publication](http://bioinformatics.oxfordjournals.org/content/26/20/2622.long)).
* [ChIPMunk](http://line.imb.ac.ru/ChIPMunk): deep and wide digging for binding motifs in ChIP-Seq data ([publication](http://bioinformatics.oxfordjournals.org/content/26/20/2622.long)).
* [HMCan](http://www.cbrc.kaust.edu.sa/hmcan/): a method for detecting chromatin modifications in cancer samples using ChIP-seq data ([publication](http://bioinformatics.oxfordjournals.org/content/29/23/2979.long)).
* seqMINER
* [Crunch project](http://crunch.unibas.ch/fcgi/crunch.fcgi)
* CSDeconv
* ...



# Resources

## Research articles

## Protocols

## ChIP-seq guidelines

* [Bailey et al., 2013](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003326). Practical Guidelines for the Comprehensive Analysis of ChIP-seq Data.

* [ENCODE & modENCODE consortia, 2012](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/). ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia.

## Tutorials

### French

* [Thomas-Chollier et al. 2012](http://www.nature.com/nprot/journal/v7/n8/full/nprot.2012.088.html). A complete workflow for the analysis of full-size ChIP-seq (and similar) data sets using peak-motifs.
* **TODO** add JvH & MTC tutos
* **TODO** Roscoff bioinformatics school: [link](http://ecole-bioinfo-aviesan.sb-roscoff.fr/archives-2014)
* [RNA-seq tutorial](http://bioinfo-fr.net/lanalyse-de-donnees-rna-seq-mode-demploi)

### English

* [Galaxy tutorial](https://wiki.galaxyproject.org/Learn)

## Databases

* [Wikipedia list](https://en.wikipedia.org/wiki/List_of_biological_databases)




# Sequencing technologies

cf slides 

# Virtualization

...

# Workflow development tools
