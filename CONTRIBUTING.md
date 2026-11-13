# Contributing to BioData Open Genome Pipeline

Thank you for contributing. **This public repository is the sole active upstream** for new work.

## Workflow

1. Fork the project or create a feature branch from `main`.
2. Keep each change focused. Document scientific assumptions, reference assembly, and expected inputs/outputs.
3. Run the relevant local self-tests and check GitHub Actions on the pull request.
4. Submit a pull request to `main`. Do not push changes from a patient workspace.

## Patient and secrets policy

- Use synthetic or public, properly licensed test fixtures only. Synthetic patient names include `Demo 1` and `Demo 2`.
- **Never** commit patient genotypes, real names, patient directory listings, report excerpts, VCF/FASTQ/BAM files, access tokens, `.env` files, or local migration records.
- Review both added and modified files, commit messages and screenshots before opening an issue or PR.
- Do not attach private genomes to public GitHub issues or discussions.
- Report vulnerabilities privately as described in [SECURITY.md](SECURITY.md).

## Scientific validity

- Preserve exact build/assembly identifiers (GRCh37 versus GRCh38).
- Distinguish measured array positions from unmeasured positions; absence is not a negative call.
- Keep VUS, weak associations and preliminary findings separate from established pathogenic evidence.
- Maintain the technical separation between open/commercial-safe references and research-only licensed sources.
- Add deterministic tests for every scientific rule or pipeline behavior change.

## License

Unless clearly stated otherwise, contributions to the repository are submitted under its [Apache-2.0 license](LICENSE). Do not contribute third-party code, data or research content unless redistribution is permitted and attribution obligations are met.
