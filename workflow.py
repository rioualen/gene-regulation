
# Sample IDs
SAMPLES = ["sample1", "sample2"]
CONTROL = ["sample1"]
TREATMENT = ["sample2"]

RESULTS = "."

rule all:
    input: expand("{treatment}_vs_{control}.bed", treatment=TREATMENT, control=CONTROL)

rule peak_calling:
    input: control="{control}.sam", treatment="{treatment}.sam"
    output: "{treatment}_vs_{control}.bed"
    shell: "touch {output}"

rule mapping:
    input: "{samples}.fastq"
    output: "{samples}.sam"
    shell: "cp {input} {output}"

rule import_fastq:
    input: "{samples}.sra"
    output: "{samples}.fastq"
    shell:"cp {input} {output}"
