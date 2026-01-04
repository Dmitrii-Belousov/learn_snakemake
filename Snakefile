import pandas as pd

configfile: "config.yml"

samples_table = pd.read_csv(
    config["samples"], sep=r"\s+", index_col=0
)

samples_list = samples_table.index.to_list()

# TO LEGACY

# def get_fastq(wildcards):
#     return {
#         "r1": samples_table.loc[wildcards.sample, "fq1"],
#         "r2": samples_table.loc[wildcards.sample, "fq2"]
#     }

# DOCKER_WRAPPER = "docker run --platform linux/amd64 --rm -v $(pwd):/io -w /io {image} {cmd}"

# def get_docker_cmd(image, cmd):
#     return DOCKER_WRAPPER.format(image=image, cmd=cmd)


rule all:
    input:
        expand("results/{sample}/{sample}.vcf", sample=samples_list),
        "results/multiqc_report.html",
        expand("results/{sample}/{sample}.variants.tsv", sample=samples_list)


include: "rules/downloaders.smk"
include: 'rules/qc.smk'
include: 'rules/mapping.smk'
include: 'rules/variant_caling.smk'


rule clean:
    shell:
        "rm -rf results/ logs/ benchmarks/ dag.png report.html"