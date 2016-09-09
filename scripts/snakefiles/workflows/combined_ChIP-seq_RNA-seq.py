## This is just a draft

# To be put in separate config file
config["metadata"]["ChIP_snakefile"] = "gene-regulation/scripts/snakefiles/workflows/ChIP-seq_workflow_SE.py"
config["metadata"]["RNA_snakefile"] = "gene-regulation/scripts/snakefiles/workflows/RNA-seq_workflow_PE.py"

config["dir"]["ChIP_analysis"] = "/data/analyses/ChIP-seq_SE_GSE41187"
config["dir"]["RNA_analysis"] = "/data/analyses/RNA-seq_PE_GSE41190"

# Subworkflows
subworkflow ChIP-seq_workflow_SE:
    workdir: config["dir"]["ChIP_analysis"]
    snakefile: config["metadata"]["ChIP_snakefile"]

subworkflow RNA-seq_workflow_PE:
    workdir: config["dir"]["RNA_analysis"]
    snakefile: config["metadata"]["RNA_snakefile"]

# Combined workflow
rule all:
    input:
        ChIP-seq_workflow_SE("genes.txt")
        RNA-seq_workflow_SE("genes.txt")
    output: # some statistics, venn diagram...
    shell:  #...

