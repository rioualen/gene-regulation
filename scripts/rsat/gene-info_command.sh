ORG=Escherichia_coli_str_k_12_substr_mg1655_GCA_000005845.2
ORG_NCBI=Escherichia_coli_K_12_substr__MG1655_uid57779
#ORG=Escherichia_coli_K_12_substr__MG1655_uid57779
cut -f 1 results/DEG/sickle_pe_q20_bowtie2_pe_sorted_name_allcounts.tab \
    | perl -pe  's|^gene_id|#gene_id|' \
    | add-gene-info -org ${ORG_NCBI} -info id -feattype gene \
    | add-gene-info -v 1 -org ${ORG} -feattype gene \
    -info id,name,ctg,left,right,strand,descr,names \
    -o genome/${ORG}_gene_info.tab
echo genome/${ORG}_gene_info.tab

awk '$4=="GO"' ${RSAT}/public_html/data/genomes/${ORG}/genome/cds_names.tab \
    | add-gene-info -org Escherichia_coli_str_k_12_substr_mg1655_GCA_000005845.2 -feattype mrna -info id,name  \
    | add-gene-info -org Escherichia_coli_str_k_12_substr_mg1655_GCA_000005845.2 -col 6 -feattype gene -info id \
    | awk '{print $7"\t"$2"\t"$6"\t"$5}'  \
    | add-gene-info -org Escherichia_coli_K_12_substr__MG1655_uid57779  -info id  \
    | add-gene-info -org Escherichia_coli_K_12_substr__MG1655_uid57779 -info id -feattype gene \
    | sort > genome/${ORG}_gene_GO.tab
echo genome/${ORG}_gene_GO.tab
