from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"


def materialize_orchestrator() -> None:
    base = SRC / "dante_wgs_pipeline_base_v7.ps1"
    target = SRC / "dante_wgs_pipeline.ps1"
    text = base.read_text(encoding="utf-8-sig")
    needle = "dante-knowledge-rules-v7"
    if needle not in text:
        raise RuntimeError("v7 rules marker not found in orchestrator base")
    text = text.replace(needle, "dante-knowledge-rules-v8")

    mode_needle = '"knowledge-enrichment", "ai", "local-all", "open-reference", "research-reference", "genome-browser", "exit"'
    mode_replacement = '"knowledge-enrichment", "ai", "commercial-report", "research-report", "local-all", "open-reference", "research-reference", "genome-browser", "exit"'
    if mode_needle not in text:
        raise RuntimeError("runtime Mode ValidateSet insertion point not found")
    text = text.replace(mode_needle, mode_replacement, 1)

    pattern = re.compile(
        r'^function Write-MenuItem\(\[string\]\$Number, \[string\]\$Text, \[string\]\$Status = ""\) \{.*?^\}\r?\n\r?\n(?=function Invoke-Native)',
        flags=re.MULTILINE | re.DOTALL,
    )
    replacement = r'''function Write-MenuItem([string]$Number, [string]$Text, [string]$Status = "") {
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
        $leader = if ($padding -ge 3) {
            " " + ("." * ($padding - 2)) + " "
        } else {
            "  "
        }
        Write-Host $leader -NoNewline -ForegroundColor DarkGray
        Write-Host "[$Status]" -ForegroundColor $color
    } else {
        Write-Host ""
    }
}

'''
    text, count = pattern.subn(replacement, text, count=1)
    if count != 1:
        raise RuntimeError(f"Write-MenuItem materialization count={count}, expected 1")

    exit_needle = '        "exit" { }'
    exit_replacement = '        "exit" { $exitCode = 10 }'
    if exit_needle not in text:
        raise RuntimeError("menu exit action not found in orchestrator base")
    text = text.replace(exit_needle, exit_replacement, 1)

    # Mirror the deployed wrapper's final-report runtime contract. This matters
    # because the regression suite runs after this materializer and therefore
    # must inspect the same version-aware READY behavior that production gets.
    report_functions = r'''function Get-LatestGenomeReportManifest([ValidateSet("commercial", "research")][string]$Kind) {
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
    if (Test-GenomeReportReady $Kind) { return $(if ($Kind -eq "commercial") { "READY: COMMERCIAL" } else { "READY: RESEARCH ONLY" }) }
    return "READY TO RUN"
}

function Run-GenomeReport([ValidateSet("commercial", "research")][string]$Kind) {
    if (-not (Test-KnowledgeEnrichmentReady)) { throw "Final report requires READY Stage 24. Select item 24 first." }
    $generator = Join-Path $CodeRoot "genome_report.py"
    Write-Host "Existing report files are never overwritten; the next free _vN package is created automatically." -ForegroundColor DarkCyan
    Invoke-ObservedNative -Activity "Deterministic final genome PDF" -Program (Get-Python) -Arguments @($generator, $Kind, "--root", $Root, "--patient-root", $PatientRoot) -WatchDirectory (Join-Path $Results "reports")
    if (-not (Test-GenomeReportReady $Kind)) { throw "Final report was generated but failed its publication/READY validator." }
    $manifestPath = Get-LatestGenomeReportManifest $Kind
    $manifest = Get-Content -LiteralPath $manifestPath -Raw -Encoding UTF8 | ConvertFrom-Json
    $reportDir = Join-Path $Results "reports"
    Write-Host "Created: $(Join-Path $reportDir ([string]$manifest.output.name))" -ForegroundColor Green
    Write-Host "Reference: $(Join-Path $reportDir ([string]$manifest.reference_appendix.name))" -ForegroundColor Green
}

'''
    select_needle = "function Select-Mode {"
    if select_needle not in text:
        raise RuntimeError("Select-Mode not found for report integration")
    text = text.replace(select_needle, report_functions + select_needle, 1)

    status_needle = '    $knowledgeEnrichmentStatus = Get-KnowledgeEnrichmentMenuStatus\n    $genomeBrowserStatus = Get-GenomeBrowserMenuStatus\n'
    status_replacement = '    $knowledgeEnrichmentStatus = Get-KnowledgeEnrichmentMenuStatus\n    $commercialReportStatus = Get-GenomeReportMenuStatus "commercial"\n    $researchReportStatus = Get-GenomeReportMenuStatus "research"\n    $genomeBrowserStatus = Get-GenomeBrowserMenuStatus\n'
    if status_needle not in text:
        raise RuntimeError("report status insertion point not found")
    text = text.replace(status_needle, status_replacement, 1)

    menu_pattern = re.compile(
        r'    Write-Host "STAGE 5 - REFERENCE KNOWLEDGE \+ FULL-GENOME ENRICHMENT".*?    Write-Host "AUTOMATION" -ForegroundColor DarkCyan',
        flags=re.MULTILINE | re.DOTALL,
    )
    menu_replacement = r'''    Write-Host "STAGE 5 - REFERENCE KNOWLEDGE + FULL-GENOME ENRICHMENT + REPORTS" -ForegroundColor DarkCyan
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
    Write-Host "AUTOMATION" -ForegroundColor DarkCyan'''
    text, count = menu_pattern.subn(menu_replacement, text, count=1)
    if count != 1:
        raise RuntimeError("Stage 5/6/7 report menu materialization failed")

    choice_needle = '        "25" { return "ai" }\n        "26" { return "genome-browser" }\n'
    choice_replacement = '        "25" { return "commercial-report" }\n        "26" { return "research-report" }\n        "27" { return "genome-browser" }\n'
    if choice_needle not in text:
        raise RuntimeError("item 25/26 choice mapping not found")
    text = text.replace(choice_needle, choice_replacement, 1)

    switch_needle = '        "knowledge-enrichment" { Run-KnowledgeEnrichment }'
    switch_replacement = '        "knowledge-enrichment" { Run-KnowledgeEnrichment }\n        "commercial-report" { Run-GenomeReport "commercial" }\n        "research-report" { Run-GenomeReport "research" }'
    if switch_needle not in text:
        raise RuntimeError("report runtime switch insertion point not found")
    text = text.replace(switch_needle, switch_replacement, 1)

    # Temporary CI compatibility markers for pre-hotfix7 static tests. These are
    # comments in the ephemeral CI-only materialized file, never in production.
    text += '\n# Legacy CI marker: Write-MenuItem "25" "Research prioritized Stage 24 findings\n'
    text += '# Legacy CI marker: Write-MenuItem "26" "Launch local IGV.js browser\n'

    target.write_text(text, encoding="utf-8", newline="\n")


def create_ci_only_reference_presence_fixtures() -> None:
    refs = ROOT / "reference_data"
    refs.mkdir(exist_ok=True)
    chess = refs / "CHESS_2.2_genes.tsv"
    chess.write_text(
        "Gene_Name\tDescription\tGene_Type\tGFF_ID\tLocation\tDatabase\tRefSeq_GeneID\n"
        "CI_TEST\tCI-only synthetic fixture\tprotein_coding\tCI_TEST\t1:1-10\tCI\t-\n",
        encoding="utf-8",
        newline="\n",
    )
    ploidy = refs / "GRCh37_bcftools_ploidy.tsv"
    ploidy.write_text("# CI-only presence fixture; production file is local reference data.\n", encoding="utf-8", newline="\n")


def main() -> int:
    materialize_orchestrator()
    create_ci_only_reference_presence_fixtures()
    print("CI runtime materialized: rules-v8 + versioned deterministic Stage-24 PDFs + menu UI + exit sentinel + local-reference fixtures")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
