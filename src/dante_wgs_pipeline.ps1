param(
    [ValidateSet(
        "menu", "setup-vep-cache", "rsid", "setup-annotsv", "annotsv", "stage1",
        "setup-fasta", "setup-clinvar", "prepare-small", "vep-full", "catalog",
        "shortlist", "stage2", "reference-targets", "reference-source-merge", "reference-align",
        "reference-call", "reference-export", "stage3", "telomere", "knowledge-enrichment", "ai",
        "commercial-report", "research-report", "local-all", "open-reference", "research-reference",
        "genome-browser", "exit"
    )]
    [string]$Mode = "menu",
    [string]$Patient = "",
    [int]$Threads = 0,
    [ValidateSet("auto", "max", "balanced")]
    [string]$PerformanceProfile = "auto",
    [switch]$ForceDownload,
    [switch]$KeepIntermediate,
    [switch]$NoPause,
    [switch]$MaterializeOnly,
    [string]$MaterializedRuntimePath = ""
)

# Deployment shim. The seeded orchestrator remains immutable; publication and
# patient-input profiles are late dot-sourced overlays in the materialized runtime.
# Static contract marker: $state.rules_schema -ne "dante-knowledge-rules-v8"

Set-StrictMode -Version 2.0
$ErrorActionPreference = "Stop"
$Base = Join-Path $PSScriptRoot "dante_wgs_pipeline_base_v7.ps1"
if (-not (Test-Path -LiteralPath $Base -PathType Leaf)) {
    throw "Missing seeded orchestrator base: $Base"
}
$Text = Get-Content -LiteralPath $Base -Raw -Encoding UTF8

$Needle = 'dante-knowledge-rules-v7'
if (-not $Text.Contains($Needle)) { throw "Cannot locate the v7 knowledge READY contract in $Base" }
$Text = $Text.Replace($Needle, 'dante-knowledge-rules-v8')
# Keep the seeded base immutable, but make every historical WGS/menu title resolve
# the release identity at runtime from release/pipeline_manifest.json.
$Text = $Text.Replace('Dante WGS local pipeline v9.0.0', 'Dante WGS local pipeline v$(Get-DanteDisplayVersion)')
$Text = $Text.Replace('Dante genomic pipeline v9.0.0', 'Dante genomic pipeline v$(Get-DanteDisplayVersion)')

$RuntimeModeNeedle = '"knowledge-enrichment", "ai", "local-all", "open-reference", "research-reference", "genome-browser", "exit"'
$RuntimeModeReplacement = '"knowledge-enrichment", "ai", "commercial-report", "research-report", "local-all", "open-reference", "research-reference", "genome-browser", "exit"'
if (-not $Text.Contains($RuntimeModeNeedle)) { throw "Cannot locate the v7 runtime Mode ValidateSet" }
$Text = $Text.Replace($RuntimeModeNeedle, $RuntimeModeReplacement)

# Keep menu statuses aligned.
$MenuPattern = '(?ms)^function Write-MenuItem\(\[string\]\$Number, \[string\]\$Text, \[string\]\$Status = ""\) \{.*?^\}\r?\n\r?\n(?=function Invoke-Native)'
$MenuReplacement = @'
function Write-MenuItem([string]$Number, [string]$Text, [string]$Status = "") {
    $prefix = ("{0,2}. {1}" -f $Number, $Text)
    Write-Host $prefix -NoNewline
    if ($Status) {
        $color = if ($Status -match 'MISMATCH|CORRUPT|ERROR|INVALID|BLOCKED|PAIRING') {
            "Red"
        } elseif ($Status -match 'PARTIAL|DOWNLOADED|PRESENT|INCOMPLETE|INTERRUPTED|PENDING|STALE') {
            "Yellow"
        } elseif ($Status -match 'READY|CURRENT|DONE') {
            "Green"
        } else {
            "DarkGray"
        }
        $statusColumn = 78
        $padding = $statusColumn - $prefix.Length
        $leader = if ($padding -ge 3) { " " + ("." * ($padding - 2)) + " " } else { "  " }
        Write-Host $leader -NoNewline -ForegroundColor DarkGray
        Write-Host "[$Status]" -ForegroundColor $color
    } else {
        Write-Host ""
    }
}

'@
$UpdatedText = [regex]::Replace($Text, $MenuPattern, $MenuReplacement, 1)
if ($UpdatedText -eq $Text) { throw "Cannot locate the v7 Write-MenuItem function in $Base" }
$Text = $UpdatedText

# start.cmd reopens a fresh menu after a completed action; exit 10 is the sentinel.
$ExitNeedle = '        "exit" { }'
$ExitReplacement = '        "exit" { $exitCode = 10 }'
if (-not $Text.Contains($ExitNeedle)) { throw "Cannot locate the menu exit action in $Base" }
$Text = $Text.Replace($ExitNeedle, $ExitReplacement)

# Runtime switch support for deterministic reports.
$SwitchNeedle = '        "knowledge-enrichment" { Run-KnowledgeEnrichment }'
$SwitchReplacement = @'
        "knowledge-enrichment" { Run-KnowledgeEnrichment }
        "commercial-report" { Run-GenomeReport "commercial" }
        "research-report" { Run-GenomeReport "research" }
'@
if (-not $Text.Contains($SwitchNeedle)) { throw "Cannot locate runtime switch for report modes" }
$Text = $Text.Replace($SwitchNeedle, $SwitchReplacement)

# Layer order is deliberate: manifest-backed version identity -> WGS report runtime
# -> clinical-risk publication input -> array VCF -> array QC -> independent
# genotype-table TXT -> generic single-VCF export. Every later layer preserves the
# prior fallback and none changes Stage-24 scientific rules.
$HookNeedle = 'Import-DotEnv $EnvFile'
$HookReplacement = @'
. (Join-Path $CodeRoot "version_runtime.ps1")
. (Join-Path $CodeRoot "report_runtime.ps1")
. (Join-Path $CodeRoot "clinical_risk_runtime.ps1")
. (Join-Path $CodeRoot "array_profile_runtime.ps1")
. (Join-Path $CodeRoot "array_normalization_qc_runtime.ps1")
. (Join-Path $CodeRoot "genotype_txt_profile_runtime.ps1")
. (Join-Path $CodeRoot "vcf_export_runtime.ps1")

Import-DotEnv $EnvFile
'@
if (-not $Text.Contains($HookNeedle)) { throw "Cannot locate runtime overlay insertion point" }
$Text = $Text.Replace($HookNeedle, $HookReplacement)

$Text = $Text.Replace(
    'placeholder or missing; only stage 6 / item 25 is affected',
    'placeholder or missing; legacy explicit AI mode only is affected'
)
$Text = $Text.Replace(
    'if ($Mode -notin @("ai", "research-reference", "knowledge-enrichment", "genome-browser"))',
    'if ($Mode -notin @("ai", "research-reference", "knowledge-enrichment", "commercial-report", "research-report", "genome-browser"))'
)
$DoneNeedle = '        if ($Mode -eq "knowledge-enrichment") { Write-Host "No OpenAI calls were made. Stage 24 may use target-sliced public/reference HTTP access for remote-indexed dbSNP/CADD data." }'
$DoneReplacement = @'
        if ($Mode -eq "knowledge-enrichment") { Write-Host "No OpenAI calls were made. Stage 24 may use target-sliced public/reference HTTP access for remote-indexed dbSNP/CADD data." }
        if ($Mode -eq "commercial-report") { Write-Host "Commercial PDF included the Stage-24 open finding layer plus local clinical-risk evidence only; no research TSV and no OpenAI/API were used." }
        if ($Mode -eq "research-report") { Write-Host "RESEARCH ONLY PDF included the Stage-24 research overlay plus local clinical-risk evidence; no OpenAI/API were used." }
'@
if (-not $Text.Contains($DoneNeedle)) { throw "Cannot locate Done-message insertion point for reports" }
$Text = $Text.Replace($DoneNeedle, $DoneReplacement)

# CI/package validation can request the exact production materialization without
# executing patient initialization or any genomic workload. The caller must supply
# an explicit output path; the normal production path below remains unchanged.
if ($MaterializeOnly) {
    if (-not $MaterializedRuntimePath) {
        throw "-MaterializeOnly requires -MaterializedRuntimePath."
    }
    $materialized = [System.IO.Path]::GetFullPath($MaterializedRuntimePath)
    $materializedParent = Split-Path -Parent $materialized
    if ($materializedParent) { New-Item -ItemType Directory -Force -Path $materializedParent | Out-Null }
    Set-Content -LiteralPath $materialized -Value $Text -Encoding UTF8
    Write-Host "Materialized production runtime: $materialized"
    exit 0
}

$Runtime = Join-Path $PSScriptRoot ".dante_wgs_pipeline_runtime_v8.ps1"
try {
    Set-Content -LiteralPath $Runtime -Value $Text -Encoding UTF8
    $Forward = @(
        "-Mode", $Mode,
        "-Threads", [string]$Threads,
        "-PerformanceProfile", $PerformanceProfile
    )
    if ($Patient) { $Forward += @("-Patient", $Patient) }
    if ($ForceDownload) { $Forward += "-ForceDownload" }
    if ($KeepIntermediate) { $Forward += "-KeepIntermediate" }
    if ($NoPause) { $Forward += "-NoPause" }
    & powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File $Runtime @Forward
    $Code = $LASTEXITCODE
} finally {
    Remove-Item -LiteralPath $Runtime -Force -ErrorAction SilentlyContinue
}
exit $Code