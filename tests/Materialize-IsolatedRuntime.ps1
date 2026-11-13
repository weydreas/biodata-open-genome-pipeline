param(
    [Parameter(Mandatory=$true)][string]$OutputRoot,
    [string]$RepositoryRoot = ""
)

Set-StrictMode -Version 2.0
$ErrorActionPreference = "Stop"

if (-not $RepositoryRoot) {
    $RepositoryRoot = Split-Path -Parent $PSScriptRoot
}
$RepositoryRoot = [System.IO.Path]::GetFullPath($RepositoryRoot)
$OutputRoot = [System.IO.Path]::GetFullPath($OutputRoot)
$sourceSrc = Join-Path $RepositoryRoot "src"
$stagedSrc = Join-Path $OutputRoot "src"
$canonicalWrapper = Join-Path $sourceSrc "dante_wgs_pipeline.ps1"
$canonicalBefore = (Get-FileHash -LiteralPath $canonicalWrapper -Algorithm SHA256).Hash

if (Test-Path -LiteralPath $OutputRoot) {
    Remove-Item -LiteralPath $OutputRoot -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $OutputRoot | Out-Null
Copy-Item -LiteralPath $sourceSrc -Destination $stagedSrc -Recurse -Force

$stagedWrapper = Join-Path $stagedSrc "dante_wgs_pipeline.ps1"
$materialized = Join-Path $stagedSrc ".ci_materialized_runtime.ps1"
& powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File $stagedWrapper `
    -MaterializeOnly -MaterializedRuntimePath $materialized
if ($LASTEXITCODE -ne 0) {
    throw "Production wrapper materialization failed with exit code $LASTEXITCODE"
}
if (-not (Test-Path -LiteralPath $materialized -PathType Leaf)) {
    throw "Materialized runtime was not created: $materialized"
}

$canonicalAfter = (Get-FileHash -LiteralPath $canonicalWrapper -Algorithm SHA256).Hash
if ($canonicalBefore -ne $canonicalAfter) {
    throw "Canonical source wrapper changed during isolated materialization."
}

Write-Host "ISOLATED MATERIALIZATION PASS" -ForegroundColor Green
Write-Host "Canonical source SHA256: $canonicalAfter"
Write-Host "Materialized runtime: $materialized"
