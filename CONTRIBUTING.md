# Contributing to Snakemake Variant Calling Workflow

Thank you for your interest in contributing to this project! This document provides guidelines for contributing.

## Getting Started

1. Fork the repository
2. Clone your fork: `git clone https://github.com/your-username/learn_snakemake.git`
3. Create a new branch: `git checkout -b feature/your-feature-name`
4. Set up the development environment: `uv sync`

## Development Guidelines

### Code Style

- Follow [PEP 8](https://pep8.org/) for Python code
- Use double quotes for strings in Python and Snakemake files
- Keep line length to 88 characters (Black formatter standard)
- Use descriptive variable and rule names

### Snakemake Best Practices

- **Rule names**: Use lowercase with underscores (e.g., `bwa_map`, `samtools_index`)
- **File paths**: Use forward slashes, even on Windows
- **Containers**: Always specify exact versions for reproducibility
- **Resources**: Specify threads and memory requirements
- **Logs**: Direct stderr/stdout to log files
- **Benchmarks**: Add benchmark directives for performance tracking

### Adding New Rules

When adding a new rule:

1. Place it in the appropriate file under `rules/`
2. Use containerized tools when possible
3. Add logging and benchmarking
4. Document input/output expectations
5. Test with sample data

Example:

```python
rule example_tool:
    input:
        "results/{sample}/input.txt"
    output:
        "results/{sample}/output.txt"
    threads: 2
    log:
        "logs/example_tool/{sample}.log"
    benchmark:
        "benchmarks/example_tool/{sample}.tsv"
    container:
        "docker://organization/tool:version"
    shell:
        "tool {input} > {output} 2> {log}"
```

### Testing

Before submitting a pull request:

1. Test the workflow with sample data
2. Ensure all rules execute successfully
3. Check that output files are generated correctly
4. Verify log files for errors
5. Run the clean rule to ensure it works

```bash
# Dry run to check workflow
snakemake -n

# Run with test data
snakemake --cores 4 --use-singularity

# Clean up
snakemake clean
```

### Documentation

- Update README.md if adding new features
- Add inline comments for complex logic
- Update CHANGELOG.md with your changes
- Include docstrings for Python functions

## Submitting Changes

### Pull Request Process

1. Ensure your code follows the style guidelines
2. Update documentation as needed
3. Add a clear description of your changes
4. Reference any related issues
5. Ensure all tests pass

### Pull Request Template

```markdown
## Description
Brief description of changes

## Motivation and Context
Why is this change needed?

## Changes Made
- Change 1
- Change 2

## Testing
How did you test these changes?

## Checklist
- [ ] Code follows style guidelines
- [ ] Documentation updated
- [ ] Tested with sample data
- [ ] CHANGELOG.md updated
```

### Commit Messages

- Use clear, descriptive commit messages
- Start with a verb in present tense (e.g., "Add", "Fix", "Update")
- Reference issues when applicable (e.g., "Fix #123")

Examples:
```
Add support for paired-end reads
Fix memory allocation in variant calling
Update BWA container to version 0.7.17
```

## Reporting Issues

### Bug Reports

Include:
- Description of the bug
- Steps to reproduce
- Expected behavior
- Actual behavior
- Environment (OS, Python version, Snakemake version)
- Error messages and log files

### Feature Requests

Include:
- Description of the feature
- Use case and motivation
- Proposed implementation (if applicable)

## Code of Conduct

- Be respectful and inclusive
- Provide constructive feedback
- Focus on what is best for the community
- Show empathy towards other community members

## Questions?

If you have questions about contributing, please open an issue with the "question" label.

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
