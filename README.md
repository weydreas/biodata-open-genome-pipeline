# BioData Open Genome Pipeline

**Open-source genomic variant analysis and reporting pipeline.**

BioData Open Genome Pipeline is a local-first Python/PowerShell workflow for processing human genomic variation data. It supports Dante-style whole-genome variant bundles, single-sample genotyping-array VCFs (including compatible Genotek-style exports), and supported genotype-table inputs. Available stages cover QC, variant normalization, VEP/ClinVar annotation, deterministic prioritization, reference-knowledge enrichment, and report generation. Some WGS-specific analyses require FASTQ, BAM, CNV, or SV inputs.

**Scope:** Research and information support only. The software is not a medical device, does not independently diagnose disease, and does not substitute for clinical assessment or confirmatory testing.

## Quick start (Windows)

Requirements: Git, Python 3, Windows PowerShell, and Docker Desktop with its engine running for Docker-based analysis steps. Large reference resources are downloaded separately and may need considerable storage.

```powershell
git clone https://github.com/weydreas/biodata-open-genome-pipeline.git
cd biodata-open-genome-pipeline
.\start.cmd
```

The launcher prepares a patient workspace and opens the interactive menu. Original VCF/FASTQ files belong under `patients/<patient name>/input/`. For example, place a compatible single-sample genotyping-array VCF under `patients/Demo 1/input/`. Do not copy actual patient files into a Git commit or issue.

`start.cmd` checks for new commits on a clean `main` checkout. Development work should use feature branches and pull requests. Data files and local settings are excluded via `.gitignore`. Some workflows need GRCh37 references, VEP caches, ClinVar, and optional knowledge resources; use the relevant menu entries to install them.

## Supported analysis scope

| Input | Supported | Not justified by the source alone |
| --- | --- | --- |
| Genotyping-array VCF / supported genotype table | Measured-locus export, carried-variant annotation, prioritization, knowledge enrichment and reports | Whole-genome coverage, unmeasured negatives, read-level calls, or telomere estimates |
| Dante-style WGS VCF bundle | Small-variant annotation and, where supplied, CNV/SV annotation and interpretation | FASTQ-specific alignment, coverage, and telomere analyses without FASTQ |
| FASTQ (where supported alongside the required inputs) | Alignment and reference-genotype reconstruction steps | Interpretation without adequate quality controls |

**Important:** A genotype absent from an array export is *unknown*, not automatically reference/negative. VUS classifications are uncertain and must not be treated as confirmed clinical findings.

## Data privacy and third-party licensing

- Patient genomes, patient metadata, `.env` secrets, local databases, runtime outputs, BAM/FASTQ/VCF, and licensed reference datasets must not be committed.
- The repository contains code, synthetic test fixtures, documentation, and deployment manifests only. Review changes carefully for accidentally embedded identifiers, absolute paths, or credentials.
- Code in this repository is licensed under the **Apache License, Version 2.0**; see [LICENSE](LICENSE). This does **not** relicense external reference datasets, clinical databases, research-only material, or third-party software.
- Some research reference data are expressly non-commercial. Do not publish or use those datasets commercially without reviewing their original licenses. See [DATA_AND_LIMITATIONS.md](DATA_AND_LIMITATIONS.md).
- The pipeline may download public reference resources over the network. Optional functionality may use external services; review configuration and privacy implications before processing personal genomes.

## Development and documentation

- [Contribution guide](CONTRIBUTING.md)
- [Security reporting](SECURITY.md)
- [Reference data and interpretation limits](DATA_AND_LIMITATIONS.md)
- [Technical documentation](docs/)
- [Current deployment manifest](release/pipeline_manifest.json)

Contributions should use a feature branch and pull request to `main`. The GitHub Actions regression gate checks important Python, PowerShell, and isolated runtime paths. This public repository is the canonical upstream for future development.
