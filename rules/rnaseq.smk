
rule trimming:
    input:
        fq1 = rules.split_sra_datasets.output.fq1,
        fq2 = rules.split_sra_datasets.output.fq2
    output:
        fq1 = "data/trimmed/{srr_id}_1.fastq.gz",
        fq2 = "data/trimmed/{srr_id}_2.fastq.gz",
        unpaired_fq1 = "data/trimmed/{srr_id}_1.unpaired.fastq.gz",
        unpaired_fq2 = "data/trimmed/{srr_id}_2.unpaired.fastq.gz"
    params:
        options = [
            "ILLUMINACLIP:data/adapters.fa:2:30:10", "LEADING:5",
            "TRAILING:5", "MINLEN:45"
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
        "{output.fq1} {output.unpaired_fq1}  "
        "{output.fq2} {output.unpaired_fq2} "
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


rule kallisto_index:
    input:
        transcriptome = rules.create_transcriptome.output.seqs
    output:
        idx = "data/references/kallisto.idx"
    params:
        kmer = "31"
    log:
        "logs/kallisto/index.log"
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index "
        "--index={output.idx} "
        "--kmer-size={params.kmer} "
        "{input.transcriptome} "
        "&> {log}"


rule kallisto_quant:
    input:
        idx = rules.kallisto_index.output.idx,
        fq1 = rules.trimming.output.fq1,
        fq2 = rules.trimming.output.fq2
    output:
        quant = "results/kallisto/{srr_id}/abundance.tsv"
    params:
        bootstrap = "50",
	outdir = "results/kallisto/{srr_id}"
    log:
        "logs/kallisto/{srr_id}.log"
    threads:
        32
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant "
        "--bias "
        "--index={input.idx} "
        "--output-dir={params.outdir} "
        "--bootstrap-samples={params.bootstrap} "
        "--threads={threads} "
        "{input.fq1} {input.fq2} "
        "&> {log}"


rule combine_gene_quantification:
    input:
        datasets = expand(
            "results/kallisto/{srr_id}/abundance.tsv",
            srr_id=config['datasets'].values()
        ),
        map = rules.generate_transcriptID_geneName.output.map
    output:
        tpm = "results/kallisto/tpm.tsv",
        est_counts = "results/kallisto/est_counts.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/combine_gene_quantification.py"
