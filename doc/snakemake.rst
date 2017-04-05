Snakemake
================================================================

This tutorial was developed assuming a unix-like architecture (Ubuntu
14.04).

Introduction
----------------------------------------------------------------

Snakemake concepts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Inspired by GNU Make: system of **rules** & **targets**
-  A rule is the *recipe* for a target
-  Rules are combined by *matching their inputs and outputs*

Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    sudo apt-get -y install python3-pip
    sudo pip3 install snakemake




Downloads for practical exercises
----------------------------------------------------------------

Ubuntu libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    sudo apt-get -y install zlib1g-dev                          # samtools (1-6)
    sudo apt-get -y install libncurses5-dev libncursesw5-dev    # samtools (1-6)

    sudo apt-get -y install r-base-core                         # Rsamtools (4-6)
    sudo pip3 install "rpy2<2.5.6"                              # Rsamtools (4-6)

    sudo pip3 install pyyaml                                    # Config management (5-6)

Tuto material
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    wget https://github.com/rioualen/gene-regulation/archive/1.0.tar.gz
    tar xvzf 1.0.tar.gz
    cd gene-regulation-1.0/doc/snakemake_tutorial

Samtools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    wget -nc http://sourceforge.net/projects/samtools/files/samtools/1.3/samtools-1.3.tar.bz2
    bunzip2 -f samtools-1.3.tar.bz2
    tar xvf samtools-1.3.tar
    cd samtools-1.3
    make
    sudo make install
    cd gene-regulation-1.0/doc/snakemake_tutorial

Rsamtools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    R

``{r eval=FALSE} source("http://bioconductor.org/biocLite.R") biocLite("Rsamtools") quit()``

Demo workflows
----------------------------------------------------------------

Workflow 1: Rules and targets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Only the first **rule** is executed by default
-  Rule ``all`` defines the **target**
-  Rule ``sam_to_bam`` automatically produces the target

\`\`\`{python, eval=FALSE} # file: workflow1.py rule all: input:
"GSM521934.bam"

rule sam\_to\_bam: input: "GSM521934.sam" output: "GSM521934.bam" shell:
"samtools view {input} > {output}" \`\`\`

In the terminal:

::

    snakemake -s workflow1/workflow1.py

Workflow 2: Introducing wildcards
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  **Wildcards** can replace variables
-  Workflow applies to list of files or samples
-  Use of the **expand** function

\`\`\`{python, eval=FALSE} # file: workflow2.py SAMPLES = ["GSM521934",
"GSM521935"]

rule all: input: expand("{sample}.bam", sample = SAMPLES)

rule sam\_to\_bam: input: "{file}.sam" output: "{file}.bam" shell:
"samtools view {input} > {output}" \`\`\`

In the terminal:

::

    snakemake -s workflow2/workflow2.py

Workflow 3: Keywords
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Rules can use a variety of **keywords**
-  An exhaustive list can be found
   `here <https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-grammar>`__

\`\`\`{python, eval=FALSE} # file: workflow3.py SAMPLES = ["GSM521934",
"GSM521935"]

rule all: input: expand("{sample}.bam", sample = SAMPLES)

rule sam\_to\_bam: input: "{file}.sam" output: "{file}.bam" params:
threads = 2 log: "{file}.log" benchmark: "{file}.json" shell: "(samtools
view -bS --threads {params.threads} {input} > {output}) > {log}" \`\`\`

In the terminal:

::

    snakemake -s workflow3/workflow3.py

Workflow 4: Combining rules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Dependencies are handled implicitly, by matching filenames
-  Commands can be executed by keywords ``run`` or ``shell``
-  Several languages: ``R``, ``bash``, ``python``

\`\`\`{python, eval=FALSE} # file: workflow4.py from snakemake.utils
import R

SAMPLES = ["GSM521934", "GSM521935"]

rule all: input: expand("{sample}\_sorted.bam", sample = SAMPLES)

rule sam\_to\_bam: input: "{file}.sam" output: "{file}.bam" params:
threads = 2 log: "{file}.log" benchmark: "{file}.json" shell: "(samtools
view -bS --threads {params.threads} {input} > {output}) > {log}"

rule bam\_sorted: input: "{file}.bam" output: "{file}\_sorted.bam" run:
R(""" library(Rsamtools) sortBam("{input}", "{output}") """) \`\`\`

In the terminal:

::

    snakemake -s workflow4/workflow4.py

Workflow 5: Configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Can be in ``json`` or in ``yml`` format
-  Acessible through the global variable **config**

\`\`\`{python, eval=FALSE} # file: workflow5.py from snakemake.utils
import R

configfile: "config.yml"

SAMPLES = config["samples"].split() OUTDIR = config["outdir"]

rule all: input: expand(OUTDIR + "{sample}\_sorted.bam", sample =
SAMPLES)

rule sam\_to\_bam: input: "{file}.sam" output: "{file}.bam" params:
threads = config["samtools"]["threads"] log: "{file}.log" benchmark:
"{file}.json" shell: "(samtools view -bS --threads {params.threads}
{input} > {output}) > {log}"

rule bam\_sorted: input: "{file}.bam" output: "{file}\_sorted.bam" run:
R(""" library(Rsamtools) sortBam("{input}", "{output}") """) \`\`\`

``{yaml, eval=FALSE} # file: config.yml   samples: "GSM521934 GSM521935"    outdir: "gene-regulation-1.0/doc/snakemake_tutorial/results/"   samtools:     threads: "2"``

In the terminal:

::

    snakemake -s workflow5/workflow5.py

Workflow 6: Separated files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  The keyword ``include`` is used to import rules

\`\`\`{python, eval=FALSE} # file: workflow6.py from snakemake.utils
import R

configfile: "config.yml"

SAMPLES = config["samples"].split() OUTDIR = config["outdir"]

include: "sam\_to\_bam.rules" include: "bam\_sorted.rules"

rule all: input: expand(OUTDIR + "{sample}\_sorted.bam", sample =
SAMPLES) \`\`\`

``{python, eval=FALSE} # file: sam_to_bam.rules rule sam_to_bam:     input: "{file}.sam"     output: "{file}.bam"     params: threads = config["samtools"]["threads"]     log: "{file}.log"     benchmark: "{file}.json"     shell: "(samtools view -bS --threads {params.threads} {input} > {output}) > {log}"``

``{python, eval=FALSE} # file: bam_sorted.rules rule bam_sorted:     input: "{file}.bam"     output: "{file}_sorted.bam"     run:         R("""         library(Rsamtools)         sortBam("{input}", "{output}")         """)``

In the terminal:

::

    snakemake -s workflow6/workflow6.py

Workflow 7: The keyword Ruleorder todo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

<!-- - ...

\`\`\`{python, eval=FALSE} # file: workflow4.py from snakemake.utils
import R

SAMPLES = ["GSM521934", "GSM521935"]

rule all: input: expand("{sample}\_sorted.bam", sample = SAMPLES)

rule sam\_to\_bam: input: "{file}.sam" output: "{file}.bam" params:
threads = 2 log: "{file}.log" benchmark: "{file}.json" shell: "(samtools
view -bS --threads {params.threads} {input} > {output}) > {log}"

rule bam\_sorted: input: "{file}.bam" output: "{file}\_sorted.bam" run:
R(""" library(Rsamtools) sortBam("{input}", "{output}") """) \`\`\` -->

Workflow 8: Combining wildcards with zip
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Workflow 9: Combining wildcards selectively
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Workflow 10: Using regular expression in wildcards
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Other
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  temp()
-  touch()
-  target/all

Bonus: generating flowcharts
----------------------------------------------------------------

::

    snakemake -s workflow6/workflow6.py --dag | dot -Tpng -o d.png
    snakemake -s workflow6/workflow6.py --rulegraph | dot -Tpng -o r.png

|Direct Acyclic Graph (DAG)| |Rulegraph|

More on snakemake...
----------------------------------------------------------------

Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `Manual <https://bitbucket.org/snakemake/snakemake/wiki/Documentation>`__
-  `FAQ <https://bitbucket.org/snakemake/snakemake/wiki/FAQ>`__
-  `Forum <https://groups.google.com/forum/#!forum/snakemake>`__

.. raw:: html

   <!--[![Snakemake](https://img.shields.io/badge/snakemake-≥3.5.2-brightgreen.svg?style=flat-square)](http://snakemake.bitbucket.org)-->

Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    apt-get install python3-pip
    pip3 install snakemake

Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Köster, Johannes and Rahmann, Sven. "Snakemake - A scalable
bioinformatics workflow engine". Bioinformatics 2012.

.. |Direct Acyclic Graph (DAG)| image:: img/dag.png
.. |Rulegraph| image:: img/rulegraph.png
