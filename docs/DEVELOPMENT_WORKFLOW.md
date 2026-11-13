# Development and regression workflow

`main` is the deployable branch consumed by `start.cmd`.

All mutable pipeline changes are developed on `develop` (or a short-lived branch) and reach `main` only after the GitHub Actions **Regression gate** passes.

The gate runs on Windows and performs:

1. Python compile check for `src/` and `tests/`.
2. Full `unittest` discovery under `tests/`.
3. Stage 24 knowledge-enrichment self-test.
4. Open-reference manager self-test.
5. Research-reference manager self-test.
6. PowerShell syntax/contract checks via `tests/Test-PowerShellSyntax.ps1`.

The local updater continues to read only `main`, so unverified development commits are never deployed into the active genome workspace.

For analytical changes, add/adjust regression tests for the changed contract before merging. Performance changes still require a controlled end-to-end benchmark in addition to unit/regression tests.
