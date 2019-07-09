configfile: "config.json"

rule all:
    input:
        config['path']['annotation'],
        expand("data/qc/{srr_id}_1_fastqc.html", srr_id=config['datasets'].values()),
        expand("results/kallisto/{srr_id}", srr_id=config['datasets'].values()),
        "results/DESeq2"

rule all_downloads:
    input:
        config['path']['annotation'],
        config['path']['genome'],
        config['path']['transcriptome'],
        expand("data/reads/{srr_id}_1.sra", srr_id=list(config['datasets'].values())[:2])

# Adding rules downloading datasets and references
include: "rules/downloads.smk"

# Adding rules for the RNA-seq pipeline
include: "rules/rnaseq.smk"

# Adding rules to run DESeq2
include: "rules/diff_expression.smk"
