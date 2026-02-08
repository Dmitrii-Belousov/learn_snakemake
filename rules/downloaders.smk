rule genome_exists:
    """
    Dummy rule to declare that the genome file exists in Azure storage.
    The file was pre-uploaded and doesn't need to be generated.
    """
    output:
        config["genome"]
    localrule: True
    shell:
        "true"  # No-op command


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

rule get_common_snp:
    output:
        vcf = "resources/common_vcf/common_all_20180418.vcf.gz",
        tbi = "resources/common_vcf/common_all_20180418.vcf.gz.tbi"
    params:
        vcf_path = "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz",
        tbi_path = "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz.tbi"
    container:
        "docker://alpine:latest"
    shell:
        """
        wget -O {output.vcf} {params.vcf_path}
        wget -O {output.tbi} {params.tbi_path}
        """