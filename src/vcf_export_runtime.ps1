# Generic single-sample VCF -> 23andMe-compatible export profile (v9.1.0).
# This overlay is sourced LAST, after WGS/array/TXT overlays, and preserves every
# previous implementation as a fallback. It activates only for one generic VCF
# that the conservative patient_input_profile currently leaves as UNKNOWN.

$script:DantePreVcfShowAutomaticCheck = (Get-Command Show-AutomaticCheck -CommandType Function).ScriptBlock
$script:DantePreVcfSelectMode = (Get-Command Select-Mode -CommandType Function).ScriptBlock
$script:DantePreVcfRunRsidExport = (Get-Command Run-RsidExport -CommandType Function).ScriptBlock
$script:DantePreVcfGetRsidExportMenuStatus = (Get-Command Get-RsidExportMenuStatus -CommandType Function).ScriptBlock

function Test-IsGenericSingleVcfProfile {
    $profile = Get-DanteInputProfile
    if ([string]$profile.kind -ne "UNKNOWN") { return $false }
    if ($null -eq $profile.source_file) { return $false }
    if ([int]$profile.inspection.sample_count -ne 1) { return $false }
    $name = [string]$profile.source_file.name
    if ($name -notmatch '(?i)\.vcf(?:\.gz)?$') { return $false }
    if (@("GRCh37", "GRCh38", "UNKNOWN") -notcontains [string]$profile.assembly) { return $false }
    return $true
}

function Get-GenericVcfSourceFile {
    $profile = Get-DanteInputProfile
    if (-not (Test-IsGenericSingleVcfProfile)) {
        throw "The active patient is not a supported generic single-sample VCF profile."
    }
    $name = [string]$profile.source_file.name
    $path = Join-Path $PatientInput $name
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "Generic VCF source is missing: $path"
    }
    return (Get-Item -LiteralPath $path)
}

function Get-Vcf23meConverterPath {
    $path = Join-Path $CodeRoot "vcf_23me_converter.py"
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "VCF -> 23andMe converter is missing: $path"
    }
    return $path
}

function Get-Vcf23meStatePath {
    return (Join-Path $Results "Dante_23andMe_v5.state.json")
}

function Test-Vcf23meExportReady {
    if (-not (Test-IsGenericSingleVcfProfile)) { return $false }
    $source = Get-GenericVcfSourceFile
    $converter = Get-Vcf23meConverterPath
    $output = Join-Path $Results "Dante_23andMe_v5.txt"
    $audit = Join-Path $Qc "Dante_23andMe_v5_unresolved.tsv"
    $report = Join-Path $Qc "Dante_23andMe_v5_report.txt"
    $statePath = Get-Vcf23meStatePath
    foreach ($required in @($output, $audit, $report, $statePath)) {
        if (-not (Test-NonEmptyFile $required)) { return $false }
    }
    try {
        $state = Get-Content -LiteralPath $statePath -Raw -Encoding UTF8 | ConvertFrom-Json
        return (
            ([string]$state.schema -eq "dante-vcf-23me-export-v1") -and
            ([string]$state.status -eq "complete") -and
            ([string]$state.pipeline_release -eq "9.1.0") -and
            ([string]$state.source_sha256 -eq (Get-Sha256 $source.FullName)) -and
            ([string]$state.converter_sha256 -eq (Get-Sha256 $converter)) -and
            ([string]$state.output_sha256 -eq (Get-Sha256 $output)) -and
            ([string]$state.audit_sha256 -eq (Get-Sha256 $audit))
        )
    } catch {
        return $false
    }
}

function Get-RsidExportMenuStatus {
    if (-not (Test-IsGenericSingleVcfProfile)) {
        return (& $script:DantePreVcfGetRsidExportMenuStatus)
    }
    $profile = Get-DanteInputProfile
    if ([string]$profile.assembly -eq "UNKNOWN") { return "ASSEMBLY UNKNOWN" }
    $paths = @(
        (Join-Path $Results "Dante_23andMe_v5.txt"),
        (Join-Path $Results "Dante_23andMe_v5.txt.part"),
        (Get-Vcf23meStatePath)
    )
    return (Get-ArtifactMenuStatus -Ready (Test-Vcf23meExportReady) -Paths $paths)
}

function Run-RsidExport {
    if (-not (Test-IsGenericSingleVcfProfile)) {
        & $script:DantePreVcfRunRsidExport
        return
    }

    $profile = Get-DanteInputProfile
    if ([string]$profile.assembly -eq "UNKNOWN") {
        throw "Generic VCF assembly is unknown. Add ##reference=GRCh37/GRCh38 to the VCF or run the converter CLI with --assembly explicitly. No coordinate guess was made."
    }
    $source = Get-GenericVcfSourceFile
    $converter = Get-Vcf23meConverterPath
    $python = Get-Python
    New-Item -ItemType Directory -Force -Path $Results, $Qc | Out-Null

    if ((Test-Vcf23meExportReady) -and -not $ForceDownload) {
        Write-Host "VCF -> 23andMe v5-compatible export is current for the same source; reusing it." -ForegroundColor Green
        return
    }

    $output = Join-Path $Results "Dante_23andMe_v5.txt"
    $audit = Join-Path $Qc "Dante_23andMe_v5_unresolved.tsv"
    $report = Join-Path $Qc "Dante_23andMe_v5_report.txt"
    Write-Step "Convert single-sample VCF to 23andMe v5-compatible SNP TXT"
    Write-Host "All valid CALLED canonical SNVs are retained. Missing rsID never causes silent loss: DANTE_* marker IDs are used." -ForegroundColor Cyan
    if ([string]$profile.assembly -eq "GRCh38") {
        Write-Host "GRCh38 coordinates will remain GRCh38. This is column-compatible, not a fake GRCh37 liftover." -ForegroundColor Yellow
    }

    $arguments = @(
        $converter, "convert", $source.FullName,
        "--output", $output,
        "--audit", $audit,
        "--report", $report,
        "--assembly", ([string]$profile.assembly)
    )
    if ($env:DANTE_VCF_23ME_NO_RSID_RESOLUTION -eq "1") {
        $arguments += "--no-rsid-resolution"
        Write-Host "NCBI/dbSNP rsID enrichment disabled by DANTE_VCF_23ME_NO_RSID_RESOLUTION=1; deterministic fallback markers remain lossless." -ForegroundColor Yellow
    }
    Invoke-ObservedNative -Activity "VCF -> 23andMe v5-compatible export" -Program $python -Arguments $arguments -WatchFile "$output.part"
    foreach ($required in @($output, $audit, $report)) {
        if (-not (Test-NonEmptyFile $required)) { throw "VCF export required artifact was not created: $required" }
    }

    $state = [ordered]@{
        schema = "dante-vcf-23me-export-v1"
        status = "complete"
        pipeline_release = "9.1.0"
        profile = "GENERIC_SINGLE_SAMPLE_VCF_EXPORT"
        source = $source.Name
        assembly = [string]$profile.assembly
        source_sha256 = Get-Sha256 $source.FullName
        converter_sha256 = Get-Sha256 $converter
        output_sha256 = Get-Sha256 $output
        audit_sha256 = Get-Sha256 $audit
        report_sha256 = Get-Sha256 $report
        whole_genome_negative_evidence = $false
        coordinate_policy = "source assembly preserved; no silent liftover"
        called_snv_loss_policy = "every representable called canonical SNV exported; missing rsID uses deterministic DANTE marker"
    }
    Set-AtomicText (Get-Vcf23meStatePath) ($state | ConvertTo-Json -Depth 5)
    if (-not (Test-Vcf23meExportReady)) {
        throw "VCF -> 23andMe export completed but failed the READY validator."
    }
    Write-Host "Created: $output" -ForegroundColor Green
    Write-Host "Audit:   $audit" -ForegroundColor Green
    Write-Host "Report:  $report" -ForegroundColor Green
}

function Show-AutomaticCheck {
    if (-not (Test-IsGenericSingleVcfProfile)) {
        & $script:DantePreVcfShowAutomaticCheck
        return
    }
    $profile = Get-DanteInputProfile
    $source = Get-GenericVcfSourceFile
    Write-Step "Automatic startup check"
    Write-CheckStatus "Code directory" (Test-Path -LiteralPath $CodeRoot -PathType Container) $CodeRoot
    Write-CheckStatus "Active patient" $true "$Patient | $PatientRoot"
    Write-CheckStatus "Input adapter" $true "generic single-sample VCF export"
    Write-CheckStatus "VCF" $true "$($source.Name) | $([Math]::Round($source.Length / 1MB, 2)) MiB"
    Write-CheckStatus "Assembly" ([string]$profile.assembly -in @("GRCh37", "GRCh38")) ([string]$profile.assembly)
    $pythonCommand = Get-Command py -ErrorAction SilentlyContinue
    if (-not $pythonCommand) { $pythonCommand = Get-Command python -ErrorAction SilentlyContinue }
    Write-CheckStatus "Python 3" ([bool]$pythonCommand) $(if ($pythonCommand) { $pythonCommand.Source } else { "not found" })
    Write-Host "Docker is not required for VCF -> 23andMe TXT conversion." -ForegroundColor DarkGray
}

function Select-Mode {
    if (-not (Test-IsGenericSingleVcfProfile)) {
        return (& $script:DantePreVcfSelectMode)
    }

    Show-AutomaticCheck
    $profile = Get-DanteInputProfile
    Write-Host ""
    Write-Host ("Dante genomic pipeline v{0}" -f (Get-DanteDisplayVersion)) -ForegroundColor Cyan
    Write-Host ("Patient: {0}" -f $Patient) -ForegroundColor Green
    Write-Host ("Workspace: {0}" -f $PatientRoot) -ForegroundColor DarkGray
    Write-Host ("Input: SINGLE-SAMPLE VCF | {0} | {1}" -f $profile.source, $profile.assembly) -ForegroundColor Cyan
    Write-Host "Scope: loci present in this VCF only; absence elsewhere is UNKNOWN." -ForegroundColor Yellow
    Write-Host ""
    Write-Host "VCF CONVERSION" -ForegroundColor DarkCyan
    Write-MenuItem "2" "Convert VCF -> 23andMe v5-compatible SNP TXT" (Get-RsidExportMenuStatus)
    Write-MenuItem "0" "Exit"

    $choice = Read-Host "Select"
    switch ($choice) {
        "2" { return "rsid" }
        "0" { return "exit" }
        default { throw "This action is not available for the active generic VCF profile: $choice" }
    }
}
