# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Comprehensive README.md with installation, usage, and workflow documentation
- LICENSE file (MIT License)
- CONTRIBUTING.md with contribution guidelines
- This CHANGELOG.md file

### Changed
- Updated project description in pyproject.toml
- Renamed `variant_caling.smk` to `variant_calling.smk` (fixed typo)
- Renamed `samtools_idex` rule to `samtools_index` (fixed typo)
- Standardized quote style in Snakefile includes (now using double quotes)

### Removed
- Removed legacy commented code from Snakefile
- Removed artifact files from git tracking (fastp.json, snpEff_summary.html)

## [0.1.0] - 2026-02-08

### Added
- Initial Snakemake workflow for variant calling
- Support for SRA data download
- Quality control with FastP
- Read mapping with BWA-MEM
- Variant calling with BCFtools
- Variant annotation with SnpEff
- Quality reporting with MultiQC
- Docker containerization for all tools
- Azure Batch execution support
- Python package configuration with pyproject.toml
- Environment management with uv

### Features
- Modular rule organization in separate files
- Automated benchmarking and logging
- Resource management (threads, memory)
- Clean rule for removing outputs

[Unreleased]: https://github.com/your-username/learn_snakemake/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/your-username/learn_snakemake/releases/tag/v0.1.0
