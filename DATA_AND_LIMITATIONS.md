# Data, licensing, and interpretation boundaries

## What Apache-2.0 covers

The [Apache License, Version 2.0](LICENSE) applies to code and original documentation in this repository to the extent contributors have rights to license them. It does **not** override licenses, copyright, trademarks, terms of service or use restrictions for third-party tools, reference datasets and other material downloaded separately.

## Reference knowledge

Reference managers may install or access resources such as ClinVar, dbSNP, VEP caches, annotations and research-only knowledge. These resources can have different attribution, redistribution and commercial-use terms, which may change upstream. Verify the actual terms for each resource and data version before distribution or commercial use.

The pipeline distinguishes open/commercial-safe reporting from **research-only / non-commercial** analysis. Do not move research-only findings into commercial reports or assume that publicly accessible data are necessarily commercially reusable.

The project does not distribute real patient genomic files or licensed reference-dataset bundles.

## Genomic interpretation

- Genotyping arrays assess only the loci actually interrogated; all others are unknown.
- A reference genotype at a measured locus is not evidence that the rest of the genome was sequenced.
- Missing calls, low-quality calls, phasing, CNV/SV, telomere length and WGS coverage require appropriate source data and quality controls.
- VUS findings do not establish pathogenicity or diagnosis. A risk association is not a definitive prediction.
- Reference assembly must be known; do not silently mix GRCh37 and GRCh38 coordinates.
- Medical decisions require qualified review and, where appropriate, confirmatory testing.

## Third-party components

See [docs/BROWSER_THIRD_PARTY_NOTICE.md](docs/BROWSER_THIRD_PARTY_NOTICE.md) for local genome-browser components. Other runtime and reference dependencies retain their own licenses and notices.
