
def sra_download_link(wildcards):
    """ Returns all FTP url for SRA download """
    url = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR"
    return '/'.join([
        url, wildcards.srr_id[:6],
        wildcards.srr_id, wildcards.srr_id + '.sra'
    ])


rule download_annotation:
    """ Downloads RefSeq GFF3 annotation """
    output:
        gff3 = temp("data/references/temp_annotation.gff3")
    params:
        link = config['download']['annotation']
    shell:
        "wget --quiet -O {output.gff3}.gz {params.link} && "
        "gzip -d {output.gff3}.gz"


rule download_datasets:
    """ Downloads the datasets from the SRA FTP server """
    output:
        sra_file = "data/reads/{srr_id}.sra"
    params:
        link = sra_download_link
    conda:
        "../envs/sra_tools.yaml"
    shell:
        "wget --quiet -O {output.sra_file} {params.link}"


rule download_genome:
    """ Downloads the genome from RefSeq FTP servers """
    output:
        genome = config['path']['genome']
    params:
        link = config['download']['genome']
    shell:
        "wget --quiet -O {output.genome}.gz {params.link} && "
        "gzip -d {output.genome}.gz "


rule modify_annotation:
    """ Transforms annotation to GTF and remove HNF4alpha definitions """
    input:
        gff3 = rules.download_annotation.output.gff3
    output:
        gtf = config['path']['annotation']
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread -TF {input.gff3} -o- | "
        "awk '/^NC/' | "
        "awk '!/HGNC:5024,/' > {output.gtf}"


rule create_transcriptome:
    """ Uses gffread to generate a transcriptome """
    input:
        genome = rules.download_genome.output.genome,
        gtf = rules.modify_annotation.output.gtf
    output:
        seqs = config['path']['transcriptome']
    params:
        hnf4a_sequences = "data/hnf4a.fa"
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread {input.gtf} -g {input.genome} -w {output.seqs}.temp && "
        "cat {output.seqs}.temp {params.hnf4a_sequences} >> {output.seqs} && "
        "rm {output.seqs}.temp"


rule split_sra_datasets:
    """ Splits the raw SRA files into the two corresponding FASTQ files """
    input:
        sra_file = rules.download_datasets.output.sra_file
    output:
        fq1 = "data/reads/{srr_id}_1.fastq.gz",
        fq2 = "data/reads/{srr_id}_2.fastq.gz"
    params:
        path = "data/reads"
    conda:
        "../envs/sra_tools.yaml"
    log:
        "../log/fastq-dump/{srr_id}.log"
    shell:
        "fastq-dump "
        "--split-files "
        "--gzip "
        "-O {params.path} {input.sra_file} "
        "&> {log}"


rule generate_transcriptID_geneName:
    """
    Generating a two-column text file containing the gene -> transcript
    relationship
    """
    input:
        gtf = rules.modify_annotation.output.gtf
    output:
        map = "data/references/gene_name.txt"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/generate_transcriptID_geneName.py"
