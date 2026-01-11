rule freebayes:
    input:
        bam = "results/{sample}/{sample}.sorted.bam",
        bai = "results/{sample}/{sample}.sorted.bam.bai",
        ref = config["genome"],
        ref_fai = "resources/genome.fasta.fai"
    threads: 1
    output:
        "results/{sample}/{sample}.vcf"
    log:
        "logs/freebayes/{sample}.log"
    benchmark:
        "benchmarks/freebayes/{sample}.tsv"
    container:
        "docker://quay.io/biocontainers/freebayes:1.3.6--hb089aa1_0"
    shell:
        "freebayes -f {input.ref} {input.bam} > {output}"


rule filter_vcf:
    input:
        "results/{sample}/{sample}.vcf"
    output:
        "results/{sample}/{sample}.filtered.vcf"
    container:
        "docker://quay.io/biocontainers/bcftools:1.16--haef29d1_2"
    shell:
        "bcftools filter -i 'QUAL>20 && DP>10' {input} > {output}"
    

rule get_snpeff_db:
    output:
        directory("resources/snpeff_data/GRCh38.86")
    container:
        "docker://quay.io/biocontainers/snpeff:5.1--hdfd78af_4"
    shell:
        """
        snpEff -Xmx4g download -dataDir \"$PWD/resources/snpeff_data\" -v GRCh38.86
        """



rule snpeff_annotate_vcf:
    input:
        vcf = "results/{sample}/{sample}.filtered.vcf",
        db = "resources/snpeff_data/GRCh38.86"
    output:
        vcf = "results/{sample}/{sample}.annotated.vcf",
        stats = "results/{sample}/{sample}.snpeff.csv"
    threads: 1
    resources:
        mem_mb = 4096
    container:
        "docker://quay.io/biocontainers/snpeff:5.1--hdfd78af_4"
    shell:
        """
        snpEff -dataDir \"$PWD/resources/snpeff_data\" \
        -Xmx4g -csvStats {output.stats} {config[snpeff_db]} \
        {input.vcf} > {output.vcf}
        """


rule unpack_annotation:
    input:
        vcf = "results/{sample}/{sample}.annotated.vcf"
    output:
        tsv = "results/{sample}/{sample}.variants.tsv"
    container:
        "docker://quay.io/biocontainers/snpsift:5.4.0a--hdfd78af_0"
    shell:
        """
        SnpSift extractFields {input.vcf} \
        CHROM POS REF ALT \
        "ANN[*].GENE" "ANN[*].EFFECT" "ANN[*].IMPACT" \
        "ANN[*].HGVS_P" "ANN[*].HGVS_C" \
        > {output.tsv}
        """
