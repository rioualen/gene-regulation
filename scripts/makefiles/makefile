################################################################################
## Reminder for the main commands (which should actually be executed by snakemake)

## Print the list of targets
help:
	@echo "Supported targets"
	@echo "	list_param		List parameters"
	@echo "	snakemake_list		List snakemake targets (call snakemake with the option -l)"
	@echo "	merge_lanes		Merge lanes to obtain one fastq per sample (more precisely one for each read size, R1, R2)"
	@echo "	bowtie2_index		Build index for bowtie2"
	@echo "	subread_index		Build index for subread"
	@echo "	one_target		Run one target specified as TARGET"
	@echo "	all			Run all tasks"
	@echo "	random_genome		Generate a randomized genome based on a Markov model trained on a reference genome."
	@echo "	merge_counts		Temporary way to merge counts for all samples"
	@echo "	clean			Remove ghost files (e.g. temporary files from emacs)"
	@echo "	fix_modes		Fix modes, which are systematically set to 777 on the ExtFat USB drive"
	@echo "	to_server_one_dir	Synchronize one directory from local machine to server"
	@echo "	from_server_one_dir	Synchronize one directory from server to local machine"

## Organism-specific parameters.
## This project is particular because each sample has to be mapped against multiple genomes
ORG=Drosophila_melanogaster
CONFIG_FILE=metadata/rna-seq/rna-seq_workflow_config.yaml

## Define snakemake command
SNAKEFILE=gene-regulation/scripts/snakefiles/workflows/rna-seq_workflow_se.py
SNAKEMAKE=snakemake -p -s ${SNAKEFILE} --configfile ${CONFIG_FILE} ${OPT}

## qsub parameters (cluster)
QSUB_JOBS=30
QSUB_QUEUE=medium
QSUB_OPT=-c "qsub -q ${QSUB_QUEUE} {params.qsub}" -j ${QSUB_JOBS}

## List parameters
list_param:
	@echo "Parameters"
	@echo "	ORG		${ORG}"
	@echo "	CONFIG_FILE	${CONFIG_FILE}"
	@echo "	SNAKEFILE	${SNAKEFILE}"
	@echo "	SNAKEMAKE	${SNAKEMAKE}"
	@echo "	QSUB_OPT	${QSUB_OPT}"

## Draw flowcharts (DAG and rulegraph)
flowcharts:
	${SNAKEMAKE} --config verbosity=0 --dag all \
		| perl -pe 's/graph\[/graph\[rankdir="LR", /' > reports/dag_${ORG}.dot; \
		dot -Tpng -o reports/dag_${ORG}.png reports/dag_${ORG}.dot; \
		dot -Tpdf -o reports/dag_${ORG}.pdf reports/dag_${ORG}.dot
	${SNAKEMAKE}  --config verbosity=0  --rulegraph all \
		>  reports/rulegraph.dot; \
		dot -Tpng -o reports/rulegraph.png  reports/rulegraph.dot; \
		dot -Tpdf -o reports/rulegraph.pdf  reports/rulegraph.dot
#	${SNAKEMAKE} flowcharts

## List snakemake targets
snakemake_list:
	${SNAKEMAKE} -l

## Build index for bowtie2
bowtie2_index:
	${SNAKEMAKE} bowtie2_build

## Build index for subread
subread_index:
	${SNAKEMAKE} subread_buildindex

## Merge lanes to obtain one fastq per sample (more precisely one for each read size, R1, R2)
merge_lanes:
	${SNAKEMAKE} merge_lanes

TARGET=all
one_target:
	${SNAKEMAKE} ${TARGET}

## Run all tasks
all:
	${SNAKEMAKE} ${QSUB_OPT} all


################################################################
## Generate a randomized genome based on a Markov model trained on a
## reference genome.
RAND_MKV=6
REF_FASTA=`awk '$$1=="fasta:" {print $$2}' ${CONFIG_FILE} | perl -pe 's|\"||g'`
REF_GTF=`awk '$$1=="features_gtf:" {print $$2}' ${CONFIG_FILE} | perl -pe 's|\"||g'`
RAND_ORG=randomized_mkv${RAND_MKV}_${ORG}
RAND_DIR=genome/${RAND_ORG}
RAND_FASTA=${RAND_DIR}/${RAND_ORG}.dna.genome.fa
RAND_GTF=${RAND_DIR}/${RAND_ORG}.gtf
RAND_CONFIG=metadata/random-genome_rna-seq_${RAND_ORG}.yml
random_genome:
	@echo
	@echo "Generating random genome"
	@echo "	ORG		${ORG}"
	@echo "	CONFIG_FIE	${CONFIG_FILE}"
	@echo "	REF_FASTA	${REF_FASTA}"
	@echo "	REF_GTF		${REF_GTF}"
	@echo "	RAND_ORG	${RAND_ORG}"
	@echo "	RAND_DIR	${RAND_DIR}"
	@echo "	RAND_FASTA	${RAND_FASTA}"
	@echo "	RAND_GTF	${RAND_GTF}"
	@echo "	RAND_CONFIG	${RAND_CONFIG}"
	@mkdir -p ${RAND_DIR}
	@cp ${REF_GTF} ${RAND_GTF}
	@echo "		${RAND_GTF}"
	random-seq -v ${V} -i ${REF_FASTA} \
		-template_format fasta -markov ${RAND_MKV} \
		-o ${RAND_FASTA}
	@echo "		${RAND_FASTA}"
	 perl -pe 's|${REF_FASTA}|${RAND_FASTA}|' ${CONFIG_FILE} \
	 	| perl -pe 's|${REF_GTF}|${RAND_GTF}|' \
	 	> ${RAND_CONFIG}
	@echo "		${RAND_CONFIG}"

## TEMPORARY: quick and dirty way to merge the count files
RESULT_DIR=results/rna-seq_results
DEG_DIR=${RESULT_DIR}/DEG
ALL_COUNTS=${DEG_DIR}/${ORG}_all_counts.tab
merge_counts:
	@mkdir -p ${DEG_DIR}
	find ${RESULT_DIR} -name '*featurecounts.tab' \
		| xargs compare-scores -v 0 -sc 7 -ic 1 -files \
		| perl -pe 's|_subread-align_pe.bam||g' \
		| perl -pe 's|\t\S+\/|\t|g' \
		> ${ALL_COUNTS}
	@echo ${ALL_COUNTS}


################################################################
## Synchronisation between local machine and server

## Remove ghost files (e.g. temporary files from emacs)
clean:
	find . -name '._*' -exec rm {} \;
	find . -name '*~' -exec rm {} \;

## Fix modes, which are systematically set to 777 on the ExtFat USB drive
fix_modes:
	find . -type d -exec chmod 755 {} \;
	find . -type f -exec chmod 644 {} \;

## Synchronize one directory from local machine to server
DIR=metadata
EXCLUDED=--exclude '*~'  \
	--exclude '._*'  \
	--exclude subread  \
	--exclude bowtie2  \
	--exclude '*.fastq'  \
	--exclude '*.fastq.gz'   \
	--exclude '*.sam' \
	--exclude '*sorted_pos.0*.bam' #\
#	--exclude '*.bam'  \
#	--exclude '*.bam.bai' #\
#	--exclude '*.bedgraph' \
#	--exclude '*.tdf'to
#SERVER_BASE=/workspace/sftp_projects/geiger/201504011_Anne_Geiger/
SERVER_BASE=/workspace/tgml/GLOSSINE/
to_server_one_dir:
	rsync -ruptvl -R ${EXCLUDED} ${OPT} ${DIR} tgml@rsat-tagc.univ-mrs.fr:${SERVER_BASE}

## Synchronize one directory from server to local machine
from_server_one_dir:
	rsync -ruptvl -R ${EXCLUDED} ${OPT} tgml@rsat-tagc.univ-mrs.fr:${SERVER_BASE}/${DIR} ..


