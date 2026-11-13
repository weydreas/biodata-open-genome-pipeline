# Deterministic Stage-24 publication runtime for Dante v9.
# Dot-sourced after the seeded orchestrator has defined its functions and paths.

function Get-LatestGenomeReportManifest([ValidateSet("commercial", "research")][string]$Kind) {
    $reportDir = Join-Path $Results "reports"
    if (-not (Test-Path -LiteralPath $reportDir -PathType Container)) { return $null }
    $baseName = if ($Kind -eq "commercial") { "Dante_genome_report_COMMERCIAL" } else { "Dante_genome_report_RESEARCH_ONLY" }
    $pattern = '^' + [regex]::Escape($baseName) + '(?:_v(?<version>[0-9]+))?\.manifest\.json$'
    $candidates = @()
    foreach ($item in (Get-ChildItem -LiteralPath $reportDir -File -ErrorAction SilentlyContinue)) {
        if ($item.Name -match $pattern) {
            $versionText = [string]$Matches['version']
            $version = if ($versionText) { [int]$versionText } else { 1 }
            $candidates += [pscustomobject]@{ Version = $version; Path = $item.FullName }
        }
    }
    if ($candidates.Count -eq 0) { return $null }
    return [string](($candidates | Sort-Object Version -Descending | Select-Object -First 1).Path)
}

function Test-GenomeReportReady([ValidateSet("commercial", "research")][string]$Kind) {
    if (-not (Test-KnowledgeEnrichmentReady)) { return $false }
    $reportDir = Join-Path $Results "reports"
    $manifestPath = Get-LatestGenomeReportManifest $Kind
    if (-not $manifestPath -or -not (Test-NonEmptyFile $manifestPath)) { return $false }
    try {
        $manifest = Get-Content -LiteralPath $manifestPath -Raw -Encoding UTF8 | ConvertFrom-Json
        $expectedMode = if ($Kind -eq "commercial") { "COMMERCIAL" } else { "RESEARCH_ONLY" }
        if (($manifest.schema -ne "dante-genome-report-v2") -or
            ([string]$manifest.mode -ne $expectedMode) -or
            ([string]$manifest.rules_schema -ne "dante-knowledge-rules-v8") -or
            ([string]$manifest.generation_mode -ne "DETERMINISTIC_STAGE24_ONLY") -or
            ([string]$manifest.publication_contract -ne "dante-report-publication-v10") -or
            ([int]$manifest.report_version -lt 1) -or
            ([bool]$manifest.versioning_policy.overwrite_existing_outputs) -or
            ([bool]$manifest.ai_used) -or
            (-not [bool]$manifest.all_stage24_findings_included)) { return $false }
        if ([string]$manifest.source_files.state.sha256 -ne (Get-Sha256 $KnowledgeEnrichmentState)) { return $false }
        if ([string]$manifest.source_files.open_findings.sha256 -ne (Get-Sha256 $KnowledgeOpenFindings)) { return $false }

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
    return "READY TO RUN"
}

function Ensure-GenomeReportRuntime {
    $python = Get-Python
    $probe = Invoke-NativeCapture $python @("-c", "import reportlab; print(reportlab.Version)")
    if (($probe.ExitCode -eq 0) -and ($probe.Output -match '^\s*4\.2\.5\s*$')) { return }
    Write-Host "Installing pinned PDF runtime reportlab 4.2.5..." -ForegroundColor DarkCyan
    Invoke-NativeWithRetry $python @("-m", "pip", "install", "--disable-pip-version-check", "--user", "reportlab==4.2.5")
    $verify = Invoke-NativeCapture $python @("-c", "import reportlab; print(reportlab.Version)")
    if (($verify.ExitCode -ne 0) -or ($verify.Output -notmatch '^\s*4\.2\.5\s*$')) {
        throw "reportlab 4.2.5 installation completed but the pinned PDF runtime is not importable."
    }
}

function Run-GenomeReport([ValidateSet("commercial", "research")][string]$Kind) {
    if (-not (Test-KnowledgeEnrichmentReady)) { throw "Final report requires READY Stage 24. Select item 24 first." }
    $generator = Join-Path $CodeRoot "genome_report.py"
    if (-not (Test-Path -LiteralPath $generator -PathType Leaf)) { throw "Genome report generator is missing: $generator" }
    $python = Get-Python
    Ensure-GenomeReportRuntime
    $args = @($generator, $Kind, "--root", $Root, "--patient-root", $PatientRoot)
    $label = if ($Kind -eq "commercial") { "Stage 25 deterministic commercial/open PDF" } else { "Stage 26 deterministic RESEARCH ONLY PDF" }
    Write-Step $label
    Write-Host "Deterministic Stage-24 renderer: no OpenAI/API calls and no top-N variant limit." -ForegroundColor Green
    Write-Host "Existing report files are never overwritten; the next free _vN package is created automatically." -ForegroundColor DarkCyan
    if ($Kind -eq "commercial") {
        Write-Host "Commercial publication gate: research-only Stage-24 files are not opened." -ForegroundColor Green
    } else {
        Write-Host "RESEARCH ONLY. The generator requires explicit RESEARCH acknowledgement before opening research-only findings." -ForegroundColor Yellow
    }
    Invoke-ObservedNative -Activity $label -Program $python -Arguments $args -WatchDirectory (Join-Path $Results "reports") -HeartbeatSeconds 20
    if (-not (Test-GenomeReportReady $Kind)) { throw "Final report was generated but failed its publication/READY validator." }
    $manifestPath = Get-LatestGenomeReportManifest $Kind
    $manifest = Get-Content -LiteralPath $manifestPath -Raw -Encoding UTF8 | ConvertFrom-Json
    $reportDir = Join-Path $Results "reports"
    Write-Host "Created: $(Join-Path $reportDir ([string]$manifest.output.name))" -ForegroundColor Green
    Write-Host "Reference: $(Join-Path $reportDir ([string]$manifest.reference_appendix.name))" -ForegroundColor Green
    Write-Host "Manifest: $manifestPath" -ForegroundColor DarkGray
}

function Select-Mode {
    Show-AutomaticCheck
    Write-Host ""
    Write-Host "Dante WGS local pipeline v9.0.0" -ForegroundColor Cyan
    Write-Host ("Patient: {0}" -f $Patient) -ForegroundColor Green
    Write-Host ("Workspace: {0}" -f $PatientRoot) -ForegroundColor DarkGray
    Write-Host ""
    $vepCacheStatus = Get-VepCacheMenuStatus
    $rsidStatus = Get-RsidExportMenuStatus
    $annotSvDataStatus = Get-AnnotSVMenuStatus
    $annotSvOutputStatus = Get-AnnotSVOutputMenuStatus
    $stage1Status = Get-CompleteStageMenuStatus -Statuses @($vepCacheStatus, $rsidStatus, $annotSvDataStatus, $annotSvOutputStatus)
    $fastaStatus = if (Test-FastaReady) { "READY" } else { "MISSING" }
    $clinVarStatus = if (Test-ClinVarReady) { "READY" } else { "MISSING" }
    $normalizedStatus = Get-NormalizedLayersMenuStatus
    $fullVepStatus = Get-FullVepMenuStatus
    $catalogStatus = Get-UnifiedCatalogMenuStatus
    $shortlistStatus = Get-ShortlistMenuStatus
    $stage2Status = Get-CompleteStageMenuStatus -Statuses @($fastaStatus, $clinVarStatus, $normalizedStatus, $fullVepStatus, $catalogStatus, $shortlistStatus)
    $referenceTargetsStatus = Get-ReferenceTargetsMenuStatus
    $referenceSourceMergeStatus = Get-ReferenceSourceMergeMenuStatus
    $referenceBamStatus = Get-ReferenceBamMenuStatus
    $referenceCallsStatus = Get-ReferenceCallsMenuStatus
    $referenceExportStatus = Get-ReferenceExportMenuStatus
    $stage3Status = Get-CompleteStageMenuStatus -Statuses @($referenceTargetsStatus, $referenceSourceMergeStatus, $referenceBamStatus, $referenceCallsStatus, $referenceExportStatus)
    $telomereStatus = Get-TelomereMenuStatus
    $allLocalStatus = Get-CompleteStageMenuStatus -Statuses @($stage1Status, $stage2Status)
    $openReferenceStatus = Get-OpenReferenceMenuStatus
    $researchReferenceStatus = Get-ResearchReferenceMenuStatus
    $knowledgeEnrichmentStatus = Get-KnowledgeEnrichmentMenuStatus
    $commercialReportStatus = Get-GenomeReportMenuStatus "commercial"
    $researchReportStatus = Get-GenomeReportMenuStatus "research"
    $genomeBrowserStatus = Get-GenomeBrowserMenuStatus

    Write-Host "STAGE 1 - MINIMAL EXPORT AND STRUCTURAL-VARIANT ANNOTATION" -ForegroundColor DarkCyan
    Write-MenuItem "1" "Install/resume VEP 116 GRCh37 cache" $vepCacheStatus
    Write-MenuItem "2" "Create Dante_23andMe_style.txt.gz" $rsidStatus
    Write-MenuItem "3" "Install AnnotSV public annotations" $annotSvDataStatus
    Write-MenuItem "4" "Annotate CNV and SV with AnnotSV" $annotSvOutputStatus
    Write-MenuItem "5" "Run complete stage 1" $stage1Status
    Write-Host ""
    Write-Host "STAGE 2 - FULL DETERMINISTIC SMALL-VARIANT ANNOTATION" -ForegroundColor DarkCyan
    Write-MenuItem "6" "Install GRCh37 reference FASTA" $fastaStatus
    Write-MenuItem "7" "Install/update current ClinVar GRCh37" $clinVarStatus
    Write-MenuItem "8" "Prepare normalized SNP, INDEL and residual layers" $normalizedStatus
    Write-MenuItem "9" "Run full VEP annotation: SNP + INDEL + residual" $fullVepStatus
    Write-MenuItem "10" "Build unified annotated small-variant catalog" $catalogStatus
    Write-MenuItem "11" "Build deterministic shortlist + genome-wide VUS/rare coding sweep" $shortlistStatus
    Write-MenuItem "12" "Run complete stage 2 (local, no API)" $stage2Status
    Write-Host ""
    Write-Host "STAGE 3 - REFERENCE-GENOME COMPLETION: ARRAY SOURCES + FASTQ" -ForegroundColor DarkCyan
    Write-MenuItem "13" "Build GRCh37 reference-genotype target panel" $referenceTargetsStatus
    Write-MenuItem "14" "Merge known array genotypes (23andMe/Genotek/Atlas; no FASTQ)" $referenceSourceMergeStatus
    Write-MenuItem "15" "Align R1 (single-end) or R1/R2 (paired-end)" $referenceBamStatus
    Write-MenuItem "16" "Call strict reference and alternate genotypes" $referenceCallsStatus
    Write-MenuItem "17" "Build tiered reference-complete export from FASTQ evidence" $referenceExportStatus
    Write-MenuItem "18" "Run complete stage 3 (array + FASTQ, local)" $stage3Status
    Write-Host ""
    Write-Host "STAGE 4 - TELOMERE CALCULATION FROM WGS FASTQ" -ForegroundColor DarkCyan
    Write-MenuItem "19" "Estimate mean telomere length and publish TXT report" $telomereStatus
    Write-Host ""
    Write-Host "STAGE 5 - REFERENCE KNOWLEDGE + FULL-GENOME ENRICHMENT + REPORTS" -ForegroundColor DarkCyan
    Write-MenuItem "22" "Install/update open reference knowledge bundle" $openReferenceStatus
    Write-MenuItem "23" "Install research only knowledge bundle (NO COMMERCIAL USE ALLOWED)" $researchReferenceStatus
    Write-MenuItem "24" "Enrich full genomic outputs with reference knowledge (NO OPENAI)" $knowledgeEnrichmentStatus
    Write-MenuItem "25" "Generate CLIENT PDF report (OPEN / COMMERCIAL-SAFE / NO AI)" $commercialReportStatus
    Write-MenuItem "26" "Generate RESEARCH ONLY PDF report (NON-COMMERCIAL / NO AI)" $researchReportStatus
    Write-Host "     LEGAL SPLIT: item 25 never opens research-only findings; item 26 requires acknowledgement; both are deterministic Stage-24 PDFs." -ForegroundColor DarkGray
    Write-Host ""
    Write-Host "STAGE 6 - LOCAL GENOME BROWSER" -ForegroundColor DarkCyan
    Write-MenuItem "27" "Launch local IGV.js browser: GRCh37 + BAM reads + calls + ClinVar" $genomeBrowserStatus
    Write-Host ""
    Write-Host "AUTOMATION" -ForegroundColor DarkCyan
    Write-MenuItem "21" "Install and run stages 1 + 2 (local, no API)" $allLocalStatus
    Write-MenuItem "0" "Exit"
    $choice = Read-Host "Select"
    switch ($choice) {
        "1" { return "setup-vep-cache" }
        "2" { return "rsid" }
        "3" { return "setup-annotsv" }
        "4" { return "annotsv" }
        "5" { return "stage1" }
        "6" { return "setup-fasta" }
        "7" { return "setup-clinvar" }
        "8" { return "prepare-small" }
        "9" { return "vep-full" }
        "10" { return "catalog" }
        "11" { return "shortlist" }
        "12" { return "stage2" }
        "13" { return "reference-targets" }
        "14" { return "reference-source-merge" }
        "15" { return "reference-align" }
        "16" { return "reference-call" }
        "17" { return "reference-export" }
        "18" { return "stage3" }
        "19" { return "telomere" }
        "21" { return "local-all" }
        "22" { return "open-reference" }
        "23" { return "research-reference" }
        "24" { return "knowledge-enrichment" }
        "25" { return "commercial-report" }
        "26" { return "research-report" }
        "27" { return "genome-browser" }
        "0" { return "exit" }
        default { throw "Unknown menu choice: $choice" }
    }
}
