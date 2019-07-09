

rule DESeq2:
    input:
        counts = rules.combine_gene_quantification.output.est_counts,
        samples = "data/references/design.tsv"
    output:
        results = directory("results/DESeq2")
    log:
        "logs/DESeq2.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/DESeq2.R"
