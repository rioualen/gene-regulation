## TEMPORARY FOR DEBUGGING: 
dir.main <- "~/BeatriceRoche/"
setwd(dir.main)
r.params.path <- "results/DEG/sickle_pe_q20_bowtie2_pe_sorted_name_params.R"  
org.db <- "org.EcK12.eg.db" ## Should be added to parameters
gene.info.file <- "genome/Escherichia_coli_str_k_12_substr_mg1655_GCA_000005845.2_gene_info.tab"
organism.names <- c("name" = "Escherichia coli",
                    "clusterProfiler" = NA,
                    "kegg"="eco")
gtf.file <- "genome/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.2.28.gtf"
gtf.source <- "ftp://ftp.ensemblgenomes.org/pub/bacteria/release-28/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/"
                                        #   pet.gene <- "b2531"
genes.of.interest <- c("b2531")
go.map.file <- 'genome/Escherichia_coli_str_k_12_substr_mg1655_GCA_000005845.2_gene_GO.tab'
go.description.file <- "genome/GO_description.tab"
