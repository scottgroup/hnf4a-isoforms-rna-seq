configfile: "config.json"

rule all:
    input:
        config['path']['annotation'],


rule all_downloads:
    input:
        config['path']['annotation'],
        config['path']['genome'],
        config['path']['transcriptome']
        # expand("data/reads/{srr_id}_1.sra", srr_id=config['datasets'].values())

# Adding rules downloading datasets and references
include: "rules/downloads.smk"

# Adding rules for the RNA-seq pipeline
include: "rules/rnaseq.smk"
