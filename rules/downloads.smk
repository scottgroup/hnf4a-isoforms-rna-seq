
def sra_download_link(wildcards):
    """ Dude """
    url = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR"
    return '/'.join([
        url, wildcards.srr_id[:6],
        wildcards.srr_id, wildcards.srr_id + '.sra'
    ])


rule download_genome:
    """ Download the genome from RefSeq FTP servers """
    output:
        genome = config['path']['genome']
    params:
        link = config['download']['genome']
    shell:
        "wget --quiet -O {output.genome}.gz {params.link} && "
        "gzip -d {output.genome}.gz "


rule download_annotation:
    """ Download RefSeq GFF3 annotation """
    output:
        gff3 = temp("data/temp/annotation.gff3")
    params:
        link = config['download']['annotation']
    shell:
        "wget --quiet -O {output.gff3}.gz {params.link} && "
        "gzip -d {output.gff3}.gz"


rule modify_annotation:
    """ Transform annotation to GTF and remove HNF4alpha definitions """
    input:
        gff3 = rules.download_annotation.output.gff3
    output:
        gtf = config['path']['annotation']
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread -TF {input.gff3} | "
        "awk '!/HGNC:5024,/' > {output.gtf}"


rule create_transcriptome:
    """ """
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


rule download_datasets:
    """ """
    output:
        sra_file = "data/reads/{srr_id}.sra"
    params:
        link = sra_download_link
    conda:
        "../envs/sra_tools.yaml"
    shell:
        "wget -O {output.sra_file} {params.link}"


rule split_sra_datasets:
    """ """
    input:
        sra_file = rules.download_datasets.output.sra_file
    output:
        fq1 = "data/reads/{srr_id}_1.fastq",
        fq2 = "data/reads/{srr_id}_2.fastq"
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
