import pandas as pd
from snakemake.utils import min_version

storage:
    provider="azure",
    account_name="snakemakestore",
    default_azure_credential=True
            

configfile: "config.yml"

samples_table = pd.read_csv(
    config["samples"], sep=r"\s+", index_col=0
)

samples_list = samples_table.index.to_list()

rule all:
    input:
        expand("results/{sample}/{sample}.somatic.vcf", sample=samples_list),
        "results/multiqc_report.html",
        expand("results/{sample}/{sample}.variants.tsv", sample=samples_list)


include: "rules/downloaders.smk"
include: "rules/qc.smk"
include: "rules/mapping.smk"
include: "rules/variant_calling.smk"


rule clean:
    shell:
        "rm -rf results/ logs/ benchmarks/ dag.png report.html"