# Snakemake Variant Calling Workflow

A reproducible Snakemake workflow for variant calling from Next-Generation Sequencing (NGS) data, with support for Azure Batch execution.

## Overview

This workflow performs end-to-end variant calling analysis from raw FASTQ files:

1. **Download**: Fetch reference genome and sequencing data from SRA
2. **Quality Control**: Adapter trimming and quality filtering with FastP
3. **Mapping**: Align reads to reference genome using BWA-MEM
4. **Variant Calling**: Call variants with BCFtools
5. **Annotation**: Annotate variants with SnpEff
6. **Quality Metrics**: Generate QC reports with MultiQC

## Features

- **Containerized**: All tools run in Docker containers for reproducibility
- **Cloud-Ready**: Azure Batch integration for scalable execution
- **Modern Stack**: Uses uv for fast Python dependency management
- **Benchmarking**: Performance tracking for all major steps

## Prerequisites

- Python 3.12+
- [uv](https://github.com/astral-sh/uv) package manager
- Docker (for local execution)
- Azure account (for cloud execution, optional)

## Installation

### 1. Clone the repository

```bash
git clone <repository-url>
cd learn_snakemake
```

### 2. Install dependencies

Using uv (recommended):

```bash
uv sync
```

Or using pip:

```bash
pip install -e .
```

## Configuration

### 1. Configure workflow settings

Edit `config.yml` to specify:

```yaml
samples: "samples.tsv"
genome: "resources/genome.fasta"
```

### 2. Prepare sample sheet

Edit `samples.tsv` with your sample information:

```tsv
sample	SRA_ID
sample1	SRR123456
sample2	SRR123457
```

Columns:
- `sample`: Sample name (used for output file naming)
- `SRA_ID`: NCBI SRA accession number

## Usage

### Local Execution

Run the complete workflow:

```bash
snakemake --cores 8 --use-singularity
```

Generate a workflow DAG:

```bash
snakemake --dag | dot -Tpng > dag.png
```

Clean up output files:

```bash
snakemake clean
```

### Azure Batch Execution

For large-scale analysis, the workflow supports Azure Batch:

```bash
snakemake --cores all \
  --executor azure-batch \
  --default-resources azure_batch_account=<account> \
                        azure_batch_pool_id=<pool>
```

Configure Azure credentials using environment variables or Azure CLI authentication.

## Workflow Structure

```
learn_snakemake/
├── Snakefile              # Main workflow definition
├── config.yml             # Configuration file
├── samples.tsv            # Sample metadata
├── rules/                 # Modular rule definitions
│   ├── downloaders.smk   # Data download rules
│   ├── qc.smk            # Quality control rules
│   ├── mapping.smk       # Read mapping rules
│   └── variant_calling.smk  # Variant calling rules
├── resources/             # Reference data (generated)
├── results/               # Analysis outputs (generated)
├── logs/                  # Execution logs (generated)
└── benchmarks/            # Performance metrics (generated)
```

## Output Files

For each sample, the workflow generates:

- `results/{sample}/{sample}.vcf` - Called variants
- `results/{sample}/{sample}.variants.tsv` - Annotated variant table
- `results/{sample}/{sample}.stats` - Alignment statistics
- `results/multiqc_report.html` - Comprehensive QC report

## Pipeline Steps

### 1. Download Reference Genome

Downloads the reference genome specified in `config.yml`.

### 2. Download Sequencing Data

Fetches FASTQ files from NCBI SRA using `fasterq-dump`.

### 3. Quality Control

- Adapter trimming with FastP
- Quality filtering (Q20 threshold)
- Generates QC metrics

### 4. Read Mapping

- BWA index creation
- BWA-MEM alignment
- SAM to sorted BAM conversion
- BAM indexing
- Alignment statistics

### 5. Variant Calling

- Reference genome indexing (faidx, dict)
- Variant calling with BCFtools mpileup + call
- VCF filtering (quality ≥20, depth ≥10)

### 6. Variant Annotation

- SnpEff annotation
- Variant effect prediction
- Summary statistics

### 7. Quality Reporting

- MultiQC aggregates all QC metrics
- Generates interactive HTML report

## Container Images

All tools use versioned Bioconda containers:

- BWA: `quay.io/biocontainers/bwa:0.7.17`
- SAMtools: `quay.io/biocontainers/samtools:1.15`
- BCFtools: `quay.io/biocontainers/bcftools:1.15.1`
- FastP: `staphb/fastp:0.23.2`
- SnpEff: `quay.io/biocontainers/snpeff:5.1`
- MultiQC: `staphb/multiqc:1.19`

## Troubleshooting

### Memory Issues

If you encounter memory errors with SnpEff, the workflow includes memory limits:

```python
resources:
    mem_mb=8000  # 8GB memory limit
```

### Missing Dependencies

Ensure all containers are pulled:

```bash
snakemake --use-singularity --singularity-pull
```

### Azure Authentication

Configure Azure credentials:

```bash
az login
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history.

## License

See [LICENSE](LICENSE) for licensing information.

## Citation

If you use this workflow in your research, please cite:

- Snakemake: Mölder et al., F1000Research 2021
- BWA: Li & Durbin, Bioinformatics 2009
- SAMtools: Li et al., Bioinformatics 2009
- BCFtools: Danecek et al., GigaScience 2021
- SnpEff: Cingolani et al., Fly 2012

## Support

For issues and questions, please open an issue on GitHub.
