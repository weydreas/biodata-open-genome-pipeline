# BioData Open Genome Pipeline

Public source snapshot of the BioData genomic analysis pipeline. Patient data, private credentials and internal Git history are not included.

## Tracked here

- launch/update scripts;
- Python/PowerShell pipeline source needed for active development;
- Stage 24 knowledge-enrichment rules;
- open/research reference installers;
- regression tests and technical documentation;
- `release/pipeline_manifest.json`, the allowlist of files deployed into the local pipeline root.

## Intentionally not tracked

Patient genomes, BAM/VCF/FASTQ, `patients/`, generated findings, local SQLite databases,
large reference datasets, API secrets, and licensed/non-commercial datasets are never committed.
Reference managers download or migrate those artifacts locally according to their manifests.

## Patient workspace privacy

All patient identities in the tracked test fixtures are synthetic (`Demo 1`, `Demo 2`, etc.).
Do not commit real patient names, genomic files, patient manifests, reports or local configuration.
The default workspace for a new installation is `patients/Demo 1/`.
Existing local patient folders, metadata, the active-patient selector and completed migration
markers are preserved by upgrades; the code never renames existing patient directories.
## Local bootstrap

Place `bootstrap_git.cmd` in the existing pipeline root and run it once. It creates a hidden
`.pipeline_repo` code cache and deploys only allowlisted code files. Thereafter `start.cmd`
checks `origin/main` and reference-bundle readiness on every launch.

Environment escape hatches:

- `DANTE_SKIP_GIT_UPDATE=1`
- `DANTE_SKIP_REFERENCE_CHECK=1`
- `DANTE_SKIP_REFERENCE_UPDATE=1`

## Important notice
This software is intended for genomic research and analysis and does not independently establish a medical diagnosis. Third-party reference datasets may have separate licensing restrictions.
