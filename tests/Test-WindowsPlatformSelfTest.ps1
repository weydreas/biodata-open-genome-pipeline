param(
    [string]$RepositoryRoot = ""
)

Set-StrictMode -Version 2.0
$ErrorActionPreference = "Stop"

if (-not $RepositoryRoot) {
    $RepositoryRoot = Split-Path -Parent $PSScriptRoot
}
$RepositoryRoot = [System.IO.Path]::GetFullPath($RepositoryRoot)

function Invoke-CheckedPowerShell([string]$ScriptPath, [string[]]$Arguments, [int[]]$AllowedExitCodes = @(0)) {
    & powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File $ScriptPath @Arguments | Out-Host
    $code = [int]$LASTEXITCODE
    if ($AllowedExitCodes -notcontains $code) {
        throw "$ScriptPath exited with code $code"
    }
    return $code
}

Write-Host "[SELFTEST] bootstrap_update.ps1 safe execution path"
$oldSkipGit = $env:DANTE_SKIP_GIT_UPDATE
$oldSkipRef = $env:DANTE_SKIP_REFERENCE_CHECK
try {
    $env:DANTE_SKIP_GIT_UPDATE = "1"
    $env:DANTE_SKIP_REFERENCE_CHECK = "1"
    $bootstrap = Join-Path $RepositoryRoot "src\bootstrap_update.ps1"
    Invoke-CheckedPowerShell $bootstrap @("-Root", $RepositoryRoot, "-PythonCommand", "python") | Out-Null
} finally {
    $env:DANTE_SKIP_GIT_UPDATE = $oldSkipGit
    $env:DANTE_SKIP_REFERENCE_CHECK = $oldSkipRef
}

Write-Host "[SELFTEST] deployed wrapper executable exit path"
$patientName = "CI Platform SelfTest"
$patientRoot = Join-Path (Join-Path $RepositoryRoot "patients") $patientName
$createdPatientsRoot = $false
$patientsRoot = Join-Path $RepositoryRoot "patients"
if (-not (Test-Path -LiteralPath $patientsRoot -PathType Container)) {
    New-Item -ItemType Directory -Force -Path $patientsRoot | Out-Null
    $createdPatientsRoot = $true
}
if (Test-Path -LiteralPath $patientRoot) {
    throw "Synthetic selftest patient path unexpectedly already exists: $patientRoot"
}
New-Item -ItemType Directory -Force -Path (Join-Path $patientRoot "input") | Out-Null
try {
    $wrapper = Join-Path $RepositoryRoot "src\dante_wgs_pipeline.ps1"
    $code = Invoke-CheckedPowerShell $wrapper @(
        "-Mode", "exit",
        "-Patient", $patientName,
        "-PerformanceProfile", "balanced",
        "-NoPause"
    ) @(10)
    if ($code -ne 10) { throw "Wrapper exit sentinel was not 10." }
} finally {
    Remove-Item -LiteralPath $patientRoot -Recurse -Force -ErrorAction SilentlyContinue
    if ($createdPatientsRoot) {
        Remove-Item -LiteralPath $patientsRoot -Force -ErrorAction SilentlyContinue
    }
}

Write-Host "[SELFTEST] Windows PowerShell publication READY contract with localized presentation label"
$fixtureRoot = Join-Path ([System.IO.Path]::GetTempPath()) ("dante-report-ready-" + [Guid]::NewGuid().ToString("N"))
$Root = $fixtureRoot
$PatientRoot = Join-Path $fixtureRoot "patients\Demo 7"
$Results = Join-Path $PatientRoot "Dante_results"
$CodeRoot = Join-Path $RepositoryRoot "src"
$KnowledgeEnrichmentDir = Join-Path $Results "knowledge_enrichment"
$KnowledgeEnrichmentState = Join-Path $KnowledgeEnrichmentDir "Dante_knowledge_enrichment.state.json"
$KnowledgeOpenFindings = Join-Path $KnowledgeEnrichmentDir "Dante_knowledge_findings_open.tsv"
$reportsDir = Join-Path $Results "reports"
$riskDir = Join-Path $Results "clinical_risk_evidence"
$script:FixtureManifestPath = Join-Path $reportsDir "Dante_genome_report_COMMERCIAL.manifest.json"

function Run-GenomeReport([string]$Kind) { return }
function Test-KnowledgeEnrichmentReady { return $true }
function Test-NonEmptyFile([string]$Path) {
    return ((Test-Path -LiteralPath $Path -PathType Leaf) -and ((Get-Item -LiteralPath $Path).Length -gt 0))
}
function Get-LatestGenomeReportManifest([string]$Kind) { return $script:FixtureManifestPath }
function Get-Sha256([string]$Path) { return ([string](Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash).ToLowerInvariant() }
function Get-Python { return "python" }
function Invoke-NativeCapture([string]$Program, [object[]]$Arguments) {
    return [pscustomobject]@{ ExitCode = 0; Output = "READY" }
}

try {
    New-Item -ItemType Directory -Force -Path $KnowledgeEnrichmentDir, $reportsDir, $riskDir | Out-Null
    [System.IO.File]::WriteAllText($KnowledgeEnrichmentState, "state", (New-Object System.Text.UTF8Encoding($false)))
    [System.IO.File]::WriteAllText($KnowledgeOpenFindings, "open", (New-Object System.Text.UTF8Encoding($false)))
    $usagePath = Join-Path $KnowledgeEnrichmentDir "Dante_knowledge_usage_manifest.json"
    [System.IO.File]::WriteAllText($usagePath, "usage", (New-Object System.Text.UTF8Encoding($false)))
    $riskJson = Join-Path $riskDir "Dante_clinical_risk_evidence.json"
    $riskTsv = Join-Path $riskDir "Dante_clinical_risk_evidence.tsv"
    $riskState = Join-Path $riskDir "Dante_clinical_risk_evidence.state.json"
    [System.IO.File]::WriteAllText($riskJson, "risk-json", (New-Object System.Text.UTF8Encoding($false)))
    [System.IO.File]::WriteAllText($riskTsv, "risk-tsv", (New-Object System.Text.UTF8Encoding($false)))
    [System.IO.File]::WriteAllText($riskState, "risk-state", (New-Object System.Text.UTF8Encoding($false)))
    $mainPdf = Join-Path $reportsDir "Dante_genome_report_COMMERCIAL.pdf"
    $referencePdf = Join-Path $reportsDir "Dante_genome_report_COMMERCIAL_COMMON_MAJOR_REFERENCE.pdf"
    [System.IO.File]::WriteAllText($mainPdf, "synthetic-pdf", (New-Object System.Text.UTF8Encoding($false)))
    [System.IO.File]::WriteAllText($referencePdf, "synthetic-reference-pdf", (New-Object System.Text.UTF8Encoding($false)))

    # Build the Russian visual label from code points so this selftest source itself
    # does not depend on PowerShell 5.1 decoding a non-ASCII UTF-8 source literal.
    $localizedLabel = -join @(
        [char]0x0432, [char]0x0435, [char]0x0440, [char]0x043D, [char]0x0443,
        [char]0x0442, [char]0x044C, [char]0x0441, [char]0x044F
    )

    $manifest = [ordered]@{
        schema = "dante-genome-report-v2"
        mode = "COMMERCIAL"
        generation_mode = "DETERMINISTIC_STAGE24_PLUS_LOCAL_RISK_EVIDENCE"
        ai_used = $false
        all_stage24_findings_included = $true
        rules_schema = "dante-knowledge-rules-v8"
        research_layer_opened = $false
        publication_contract = "dante-report-clickable-toc-backlinks-v17"
        report_version = 1
        versioning_policy = [ordered]@{ overwrite_existing_outputs = $false }
        clinical_risk_evidence_schema = "dante-clinical-risk-evidence-v1"
        table_of_contents = [ordered]@{
            schema = "dante-report-clickable-toc-v2"
            clickable_internal_links = $true
            pdf_outline_bookmarks = $true
            heading_backlinks_to_contents = $true
            backlink_destination = "dante_toc_home"
            backlink_label = $localizedLabel
        }
        source_files = [ordered]@{
            state = [ordered]@{ sha256 = (Get-Sha256 $KnowledgeEnrichmentState) }
            open_findings = [ordered]@{ sha256 = (Get-Sha256 $KnowledgeOpenFindings) }
            usage = [ordered]@{ sha256 = (Get-Sha256 $usagePath) }
            clinical_risk_evidence = [ordered]@{ sha256 = (Get-Sha256 $riskJson) }
            clinical_risk_state = [ordered]@{ sha256 = (Get-Sha256 $riskState) }
        }
        output = [ordered]@{ name = (Split-Path -Leaf $mainPdf); sha256 = (Get-Sha256 $mainPdf) }
        reference_appendix = [ordered]@{ name = (Split-Path -Leaf $referencePdf); sha256 = (Get-Sha256 $referencePdf) }
    }
    $manifestJson = $manifest | ConvertTo-Json -Depth 12
    [System.IO.File]::WriteAllText($script:FixtureManifestPath, $manifestJson, (New-Object System.Text.UTF8Encoding($false)))

    . (Join-Path $RepositoryRoot "src\clinical_risk_runtime.ps1")
    if (-not (Test-GenomeReportReady "commercial")) {
        throw "Synthetic v17 report manifest failed Test-GenomeReportReady on Windows PowerShell."
    }
} finally {
    Remove-Item -LiteralPath $fixtureRoot -Recurse -Force -ErrorAction SilentlyContinue
}

Write-Host "WINDOWS PLATFORM SELFTEST PASS | bootstrap + wrapper + synthetic publication READY contract" -ForegroundColor Green
exit 0
