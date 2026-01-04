rule download_sra:
    output:
        fq1 = "results/{sample}/{sample}_R1.fastq.gz",
        fq2 = "results/{sample}/{sample}_R2.fastq.gz"
    params:
        sra_id = lambda wildcards: samples_table.loc[wildcards.sample, "SRA_ID"]
    threads: 4
    container:
        "docker://quay.io/biocontainers/sra-tools:3.2.1--h4304569_1"
    shell:
        """
        fasterq-dump {params.sra_id} \
            --outdir results/{wildcards.sample}/ \
            --split-3 \
            --threads {threads} \
            --temp /tmp \
            --progress

        mv results/{wildcards.sample}/{params.sra_id}_1.fastq \
        results/{wildcards.sample}/{wildcards.sample}_R1.fastq
        mv results/{wildcards.sample}/{params.sra_id}_2.fastq \
        results/{wildcards.sample}/{wildcards.sample}_R2.fastq

        gzip results/{wildcards.sample}/{wildcards.sample}_R1.fastq
        gzip results/{wildcards.sample}/{wildcards.sample}_R2.fastq
        """