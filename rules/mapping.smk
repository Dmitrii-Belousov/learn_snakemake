rule bwa_index:
    input:
        config["genome"]
    output:
        multiext(config["genome"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    shell:
        get_docker_cmd(
            image="quay.io/biocontainers/bwa:0.7.17--hed695b0_7",
            cmd="bwa index {input}"
        )


rule bwa_map:
    input:
        r1 = "results/{sample}_R1.filtered.fastq.gz",
        r2 = "results/{sample}_R2.filtered.fastq.gz",
        ref = config["genome"],
        idx = multiext(config["genome"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        "results/{sample}.sam"
    threads: 4
    benchmark:
        "benchmarks/bwa/{sample}.tsv"
    log:
        "logs/bwa/{sample}.log"
    shell:
        get_docker_cmd(
            image="quay.io/biocontainers/bwa:0.7.17--hed695b0_7",
            cmd="bwa mem {input.ref} {input.r1} {input.r2} > {output} -t {threads} 2> {log}"
        )

rule samtools_faidx:
    input:
        config["genome"]
    output:
        "resources/genome.fasta.fai"
    shell:
        get_docker_cmd(
            image="quay.io/biocontainers/samtools:1.15--h3843a85_0",
            cmd="samtools faidx {input}"
        )


rule samtools_dict:
    input:
        config["genome"]
    output:
        "resources/genome.dict"
    shell:
        get_docker_cmd(
            image="quay.io/biocontainers/samtools:1.15--h3843a85_0",
            cmd="samtools dict {input} -o {output}"
        )


rule samtools_sort:
    input:
        "results/{sample}.sam"
    output:
        "results/{sample}.sorted.bam"
    threads: 4
    shell:
        get_docker_cmd(
            image="quay.io/biocontainers/samtools:1.15--h3843a85_0",
            cmd="samtools sort -o {output} {input} -@ {threads}"
        )

rule samtools_idex:
    input:
        "results/{sample}.sorted.bam"
    output:
        "results/{sample}.sorted.bam.bai"
    shell:
        get_docker_cmd(
            image="quay.io/biocontainers/samtools:1.15--h3843a85_0",
            cmd="samtools index {input}"
        )


rule samtools_stats:
    input:
        bam = "results/{sample}.sorted.bam",
        bai = "results/{sample}.sorted.bam.bai"
    output:
        "results/{sample}.stats"
    shell:
        get_docker_cmd(
            image="quay.io/biocontainers/samtools:1.15--h3843a85_0",
            cmd="samtools stats {input.bam} > {output}"
        )
