# Security and patient privacy

This software handles highly sensitive human genomic data. Treat every patient file, genotype, report and identifying metadata as confidential.

## Report vulnerabilities privately

Do not post exploitable details, API keys or patient data in a public GitHub issue or pull request. If GitHub's **Report a vulnerability** / private security advisory option is available for this repository, use it. Otherwise contact the repository owner (`weydreas`) privately before disclosing details.

For a report, provide a minimal **synthetic** reproduction, affected version or commit, impact and suggested mitigation. Never attach real genomic data.

## Safe operation

- Keep original patient input and generated results outside Git; follow `.gitignore` and review `git status` before every commit.
- Keep `.env`, tokens and local credentials off GitHub. Rotate and revoke any credential that was exposed.
- Avoid attaching raw stdout, stack traces or screenshots that reveal patient names or local identifying paths.
- The local genome browser uses a localhost endpoint; do not publicly expose it without separate security review.
- Some stages download external reference data or use optional online services. Review egress and privacy requirements before using personal genomic data.
