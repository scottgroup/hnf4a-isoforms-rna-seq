
rule trimming:
    input:
        fq1 = rules.split_sra_datasets.output.fq1,
        fq2 = rules.split_sra_datasets.output.fq2
    output:
        fq1 = "data/trimmed/{srr_id}_1.fastq",
        fq2 = "data/trimmed/{srr_id}_2.fastq",
        unpaired_fq1 = "data/trimmed/{srr_id}_1.unpaired.fastq",
        unpaired_fq2 = "data/trimmed/{srr_id}_2.unpaired.fastq"
    params:
        options = [
            "ILLUMINACLIP:data/adapters.fa:2:30:10", "LEADING:25",
            "TRAILING:25", "MINLEN:45"
        ]
    log:
        "logs/trimmomatic/{srr_id}.log"
    threads:
        32
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        "trimmomatic PE "
        "-threads {threads} "
        "-phred33 "
        "{input.fq1} {input.fq2} "
        "{output.fq1} {output.fq2} "
        "{output.unpaired_fq1} {output.unpaired_fq2} "
        "{params.options} "
        "&> {log}"

rule qc:
    input:
        fq1 = rules.trimming.output.fq1,
        fq2 = rules.trimming.output.fq2,
        unpaired_fq1 = rules.trimming.output.unpaired_fq1,
        unpaired_fq2 = rules.trimming.output.unpaired_fq2,
    output:
        fq1_out = "data/qc/{srr_id}_1_fastqc.html"
    params:
        out_dir = "data/qc"
    log:
        "logs/fastqc/{srr_id}.log"
    threads:
        32
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc "
        "--outdir {params.out_dir} "
        "--format fastq "
        "--threads {threads} "
        "{input.fq1} {input.fq2} "
        "{input.unpaired_fq1} {input.unpaired_fq2} "
        "&> {log}"
