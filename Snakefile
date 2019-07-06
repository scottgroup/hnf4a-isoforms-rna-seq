configfile: "config.json"

rule all:
    input:
        config['path']['annotation'],
        config['path']['genome']

# Adding rules downlaoding datasets and references
include: "rules/downloads.smk"
