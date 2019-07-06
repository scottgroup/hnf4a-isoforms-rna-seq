
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

rule download_datasets:
