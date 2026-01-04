rule fast_qc:
    input:
        unpack(get_fastq)
    output:
        o1 = "results/{sample}_R1.filtered.fastq.gz",
        o2 = "results/{sample}_R2.filtered.fastq.gz",
        json = "results/{sample}_fastp.json",
        html = "results/{sample}_fastp.html"
    container:
        "docker://quay.io/biocontainers/fastp:0.20.1--h8b12597_0"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.o1} -O {output.o2} -h {output.html} -j {output.json}"


rule multiqc:
    input:
        expand("results/{sample}_fastp.json", sample=samples_list),
        expand("results/{sample}.stats", sample=samples_list),
        expand("results/{sample}.snpeff.csv", sample=samples_list)
    output:
        "results/multiqc_report.html"
    container:
        "docker://staphb/multiqc:1.8"
    shell:
        "multiqc results/ -n multiqc_report.html -o results/ -f"