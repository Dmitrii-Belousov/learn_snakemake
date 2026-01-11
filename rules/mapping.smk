rule bwa_index:
    input:
        config["genome"]
    output:
        multiext(config["genome"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    container: 
        "docker://quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
    shell:
        "bwa index {input}"


rule bwa_map:
    input:
        r1 = "results/{sample}/{sample}_R1.filtered.fastq.gz",
        r2 = "results/{sample}/{sample}_R2.filtered.fastq.gz",
        ref = config["genome"],
        idx = multiext(config["genome"], ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        "results/{sample}/{sample}.sorted.bam"
    threads: 4
    params:
        mem_gb = "1G"
    benchmark:
        "benchmarks/bwa/{sample}.tsv"
    log:
        "logs/bwa/{sample}.log"
    container:
        "docker://quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:f45ad9036aa41bb10f875a330fa877d8869018a1-0"
        #"docker://quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
    shell:
        """
	bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | \
	samtools sort -@ {threads} -m {params.mem_gb} -o {output} -
	"""

rule samtools_faidx:
    input:
        config["genome"]
    output:
        "resources/genome.fasta.fai"
    container:
        "docker://quay.io/biocontainers/samtools:1.15--h3843a85_0"
    shell:
        "samtools faidx {input}"


rule samtools_dict:
    input:
        config["genome"]
    output:
        "resources/genome.dict"
    container: 
        "docker://quay.io/biocontainers/samtools:1.15--h3843a85_0"
    shell:
        "samtools dict {input} -o {output}"


#rule samtools_sort:
#    input:
#        "results/{sample}/{sample}.sam"
#    output:
#        "results/{sample}/{sample}.sorted.bam"
#    threads: 4
#    container:
#        "docker://quay.io/biocontainers/samtools:1.15--h3843a85_0"
#    shell:
#        "samtools sort -o {output} {input} -@ {threads}"

rule samtools_idex:
    input:
        "results/{sample}/{sample}.sorted.bam"
    output:
        "results/{sample}/{sample}.sorted.bam.bai"
    container:
        "docker://quay.io/biocontainers/samtools:1.15--h3843a85_0"
    shell:
        "samtools index {input}"


rule samtools_stats:
    input:
        bam = "results/{sample}/{sample}.sorted.bam",
        bai = "results/{sample}/{sample}.sorted.bam.bai"
    output:
        "results/{sample}/{sample}.stats"
    container:
        "docker://quay.io/biocontainers/samtools:1.15--h3843a85_0"
    shell:
        "samtools stats {input.bam} > {output}"
