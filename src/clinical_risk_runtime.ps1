# Clinical-risk evidence publication layer for Dante v9.
# Dot-sourced after report_runtime.ps1 and before patient-input profile overlays.
# It does not change Stage-24 scientific findings; it makes reports depend on a
# deterministic local risk-evidence artifact built from Stage 24 + local ClinGen.

$script:DanteOriginalRunGenomeReport = (Get-Command Run-GenomeReport -CommandType Function).ScriptBlock

$ClinicalRiskEvidenceDir = Join-Path $Results "clinical_risk_evidence"
$ClinicalRiskEvidenceJson = Join-Path $ClinicalRiskEvidenceDir "Dante_clinical_risk_evidence.json"
$ClinicalRiskEvidenceTsv = Join-Path $ClinicalRiskEvidenceDir "Dante_clinical_risk_evidence.tsv"
$ClinicalRiskEvidenceState = Join-Path $ClinicalRiskEvidenceDir "Dante_clinical_risk_evidence.state.json"
$ClinicalRiskEvidenceSchema = "dante-clinical-risk-evidence-v1"
$ClinicalRiskEvidenceStateSchema = "dante-clinical-risk-evidence-state-v1"
$ClinicalRiskPublicationContract = "dante-report-clickable-toc-backlinks-v17"
$ClinicalRiskTocSchema = "dante-report-clickable-toc-v2"
$ClinicalRiskGenerationMode = "DETERMINISTIC_STAGE24_PLUS_LOCAL_RISK_EVIDENCE"

function Get-ClinicalRiskEvidenceHelper {
    $path = Join-Path $CodeRoot "clinical_risk_evidence_hotfix1.py"
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "Clinical-risk evidence helper is missing: $path"
    }
    return $path
}

function Test-ClinicalRiskEvidenceReady {
    if (-not (Test-KnowledgeEnrichmentReady)) { return $false }
    foreach ($required in @($ClinicalRiskEvidenceJson, $ClinicalRiskEvidenceTsv, $ClinicalRiskEvidenceState)) {
        if (-not (Test-NonEmptyFile $required)) { return $false }
    }
    try {
        $helper = Get-ClinicalRiskEvidenceHelper
        $python = Get-Python
        $probe = Invoke-NativeCapture $python @(
            $helper, "status",
            "--root", $Root,
            "--patient-root", $PatientRoot
        )
        return (($probe.ExitCode -eq 0) -and ($probe.Output -match "READY"))
    } catch {
        return $false
    }
}

function Ensure-ClinicalRiskEvidence {
    if (Test-ClinicalRiskEvidenceReady) {
        Write-Host "Clinical-risk evidence is current; reusing local artifact." -ForegroundColor Green
        return
    }
    if (-not (Test-KnowledgeEnrichmentReady)) {
        throw "Clinical-risk evidence requires READY Stage 24. Select item 24 first."
    }
    $helper = Get-ClinicalRiskEvidenceHelper
    $python = Get-Python
    New-Item -ItemType Directory -Force -Path $ClinicalRiskEvidenceDir | Out-Null
    $label = "Clinical risk evidence: Stage 24 + local ClinGen Actionability"
    Write-Step $label
    Write-Host "No OpenAI/API calls. Patient rsIDs, genes and genotypes remain local." -ForegroundColor Green
    Write-Host "Bare ClinGen likelihood scores are not converted into risk; only explicit quantitative estimates are accepted." -ForegroundColor DarkCyan
    Invoke-ObservedNative -Activity $label -Program $python -Arguments @(
        $helper, "build",
        "--root", $Root,
        "--patient-root", $PatientRoot
    ) -WatchDirectory $ClinicalRiskEvidenceDir -HeartbeatSeconds 20
    if (-not (Test-ClinicalRiskEvidenceReady)) {
        throw "Clinical-risk evidence was generated but failed its READY validator."
    }
}

function Test-GenomeReportReady([ValidateSet("commercial", "research")][string]$Kind) {
    if (-not (Test-KnowledgeEnrichmentReady)) { return $false }
    if (-not (Test-ClinicalRiskEvidenceReady)) { return $false }
    $reportDir = Join-Path $Results "reports"
    $manifestPath = Get-LatestGenomeReportManifest $Kind
    if (-not $manifestPath -or -not (Test-NonEmptyFile $manifestPath)) { return $false }
    try {
        $manifest = Get-Content -LiteralPath $manifestPath -Raw -Encoding UTF8 | ConvertFrom-Json
        $expectedMode = if ($Kind -eq "commercial") { "COMMERCIAL" } else { "RESEARCH_ONLY" }
        if (($manifest.schema -ne "dante-genome-report-v2") -or
            ([string]$manifest.mode -ne $expectedMode) -or
            ([string]$manifest.rules_schema -ne "dante-knowledge-rules-v8") -or
            ([string]$manifest.generation_mode -ne $ClinicalRiskGenerationMode) -or
            ([string]$manifest.publication_contract -ne $ClinicalRiskPublicationContract) -or
            ([string]$manifest.clinical_risk_evidence_schema -ne $ClinicalRiskEvidenceSchema) -or
            ([string]$manifest.table_of_contents.schema -ne $ClinicalRiskTocSchema) -or
            (-not [bool]$manifest.table_of_contents.clickable_internal_links) -or
            (-not [bool]$manifest.table_of_contents.pdf_outline_bookmarks) -or
            (-not [bool]$manifest.table_of_contents.heading_backlinks_to_contents) -or
            ([string]$manifest.table_of_contents.backlink_destination -ne "dante_toc_home") -or
            ([int]$manifest.report_version -lt 1) -or
            ([bool]$manifest.versioning_policy.overwrite_existing_outputs) -or
            ([bool]$manifest.ai_used) -or
            (-not [bool]$manifest.all_stage24_findings_included)) { return $false }

        # Do not compare the localized visual backlink label here. Windows PowerShell 5.1
        # can decode UTF-8 scripts without BOM through the active ANSI code page, making
        # non-ASCII source literals unsuitable as machine-readable READY contract values.
        # Backlink semantics are machine-validated by the boolean + ASCII destination;
        # the Python PDF regression separately verifies the visual localized label.

        if ([string]$manifest.source_files.state.sha256 -ne (Get-Sha256 $KnowledgeEnrichmentState)) { return $false }
        if ([string]$manifest.source_files.open_findings.sha256 -ne (Get-Sha256 $KnowledgeOpenFindings)) { return $false }
        if (-not ($manifest.source_files.PSObject.Properties.Name -contains "clinical_risk_evidence")) { return $false }
        if (-not ($manifest.source_files.PSObject.Properties.Name -contains "clinical_risk_state")) { return $false }
        if ([string]$manifest.source_files.clinical_risk_evidence.sha256 -ne (Get-Sha256 $ClinicalRiskEvidenceJson)) { return $false }
        if ([string]$manifest.source_files.clinical_risk_state.sha256 -ne (Get-Sha256 $ClinicalRiskEvidenceState)) { return $false }

        $pdf = Join-Path $reportDir ([string]$manifest.output.name)
        $referencePdf = Join-Path $reportDir ([string]$manifest.reference_appendix.name)
        if (-not (Test-NonEmptyFile $pdf) -or -not (Test-NonEmptyFile $referencePdf)) { return $false }
        if ([string]$manifest.output.sha256 -ne (Get-Sha256 $pdf)) { return $false }
        if ([string]$manifest.reference_appendix.sha256 -ne (Get-Sha256 $referencePdf)) { return $false }

        if ($Kind -eq "commercial") {
            if ([bool]$manifest.research_layer_opened) { return $false }
            if ($manifest.source_files.PSObject.Properties.Name -contains "research_overlays") { return $false }
            $usage = Join-Path $KnowledgeEnrichmentDir "Dante_knowledge_usage_manifest.json"
            if ([string]$manifest.source_files.usage.sha256 -ne (Get-Sha256 $usage)) { return $false }
        } else {
            if (-not [bool]$manifest.research_layer_opened) { return $false }
            $researchFindings = Join-Path $KnowledgeEnrichmentDir "Dante_knowledge_findings_research.tsv"
            $usage = Join-Path $KnowledgeEnrichmentDir "Dante_knowledge_usage_research_manifest.json"
            if ([string]$manifest.source_files.research_overlays.sha256 -ne (Get-Sha256 $researchFindings)) { return $false }
            if ([string]$manifest.source_files.usage.sha256 -ne (Get-Sha256 $usage)) { return $false }
        }
        return $true
    } catch {
        return $false
    }
}

function Get-GenomeReportMenuStatus([ValidateSet("commercial", "research")][string]$Kind) {
    if (-not (Test-KnowledgeEnrichmentReady)) { return "BLOCKED: ITEM 24" }
    if (Test-GenomeReportReady $Kind) {
        return $(if ($Kind -eq "commercial") { "READY: COMMERCIAL" } else { "READY: RESEARCH ONLY" })
    }
    if (-not (Test-ClinicalRiskEvidenceReady)) {
        return "READY TO RUN: RISK EVIDENCE PENDING"
    }
    return "READY TO RUN"
}

function Run-GenomeReport([ValidateSet("commercial", "research")][string]$Kind) {
    Ensure-ClinicalRiskEvidence
    Write-Host "Publication input: Stage 24 + deterministic local clinical-risk evidence." -ForegroundColor DarkCyan
    Write-Host "PDF navigation: clickable contents + outline bookmarks + heading return links." -ForegroundColor DarkCyan
    Write-Host "PDF generation itself remains NO AI / NO API." -ForegroundColor Green
    & $script:DanteOriginalRunGenomeReport $Kind
}
