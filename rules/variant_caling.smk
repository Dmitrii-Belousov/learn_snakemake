rule freebayes:
    input:
        bam = "results/{sample}.sorted.bam",
        bai = "results/{sample}.sorted.bam.bai",
        ref = config["genome"],
        ref_fai = "resources/genome.fasta.fai"
    threads: 1
    output:
        "results/{sample}.vcf"
    log:
        "logs/freebayes/{sample}.log"
    benchmark:
        "benchmarks/freebayes/{sample}.tsv"
    shell:
        get_docker_cmd(
            image="quay.io/biocontainers/freebayes:1.3.6--hb089aa1_0",
            cmd="freebayes -f {input.ref} {input.bam} > {output}"
        )


rule snpeff_annotate_vcf:
    input:
        vcf = "results/{sample}.vcf"
    output:
        vcf = "results/{sample}.annotated.vcf",
        stats = "results/{sample}.snpeff.csv"
    threads: 1
    resources:
        mem_mb = 4096
    shell:
        get_docker_cmd(
            image="quay.io/biocontainers/snpeff:5.1--hdfd78af_4",
            cmd="snpEff -Xmx4g -csvStats {output.stats} {config[snpeff_db]} {input.vcf} > {output.vcf}"
        )


rule unpack_annotation:
    input:
        vcf = "results/{sample}.annotated.vcf"
    output:
        tsv = "results/{sample}.variants.tsv"
    shell:
        get_docker_cmd(
            image="quay.io/biocontainers/snpsift:5.4.0a--hdfd78af_0",
            cmd="""SnpSift extractFields {input.vcf} \
            CHROM POS REF ALT \
            "ANN[*].GENE" "ANN[*].EFFECT" "ANN[*].IMPACT" \
            "ANN[*].HGVS_P" "ANN[*].HGVS_C" \
            > {output.tsv}"""
        )