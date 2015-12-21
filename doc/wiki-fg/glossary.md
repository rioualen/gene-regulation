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


