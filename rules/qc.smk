rule fast_qc:
    input:
        fq1 = "results/{sample}/{sample}_R1.fastq.gz",
        fq2 = "results/{sample}/{sample}_R2.fastq.gz"
    output:
        o1 = "results/{sample}/{sample}_R1.filtered.fastq.gz",
        o2 = "results/{sample}/{sample}_R2.filtered.fastq.gz",
        json = "results/{sample}/{sample}_fastp.json",
        html = "results/{sample}/{sample}_fastp.html"
    container:
        "docker://quay.io/biocontainers/fastp:0.20.1--h8b12597_0"
    shell:
        "fastp -i {input.fq1} -I {input.fq2} -o {output.o1} -O {output.o2} -h {output.html} -j {output.json}"


rule multiqc:
    input:
        expand("results/{sample}/{sample}_fastp.json", sample=samples_list),
        expand("results/{sample}/{sample}.stats", sample=samples_list),
        expand("results/{sample}/{sample}.snpeff.csv", sample=samples_list)
    output:
        html="results/multiqc_report.html",
        archive="results/multiqc_data.html"
    container:
        "docker://staphb/multiqc:1.8"
    shell:
        """
        multiqc {input} -n multiqc_report.html -o . -f
        tar -czf {output.archive} multiqc_report_data/
        """