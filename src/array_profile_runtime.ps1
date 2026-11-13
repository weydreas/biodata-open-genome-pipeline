# Genotyping-array patient runtime profile for Dante v9.
# Dot-sourced only after the seeded orchestrator has defined its functions/paths.

$script:DanteInputProfileCache = $null
$script:DanteOriginalShowAutomaticCheck = (Get-Command Show-AutomaticCheck -CommandType Function).ScriptBlock
$script:DanteOriginalSelectMode = (Get-Command Select-Mode -CommandType Function).ScriptBlock
$script:DanteOriginalRunRsidExport = (Get-Command Run-RsidExport -CommandType Function).ScriptBlock
$script:DanteOriginalGetRsidExportMenuStatus = (Get-Command Get-RsidExportMenuStatus -CommandType Function).ScriptBlock
$script:DanteOriginalPrepareNormalizedLayers = (Get-Command Prepare-NormalizedLayers -CommandType Function).ScriptBlock
$script:DanteOriginalTestNormalizedLayersReady = (Get-Command Test-NormalizedLayersReady -CommandType Function).ScriptBlock
$script:DanteOriginalGetNormalizedLayersMenuStatus = (Get-Command Get-NormalizedLayersMenuStatus -CommandType Function).ScriptBlock
$script:DanteOriginalRunFullVepAnnotation = (Get-Command Run-FullVepAnnotation -CommandType Function).ScriptBlock
$script:DanteOriginalRunStage1 = (Get-Command Run-Stage1 -CommandType Function).ScriptBlock
$script:DanteOriginalGetAnnotSVOutputMenuStatus = (Get-Command Get-AnnotSVOutputMenuStatus -CommandType Function).ScriptBlock
$script:DanteOriginalRunKnowledgeEnrichment = (Get-Command Run-KnowledgeEnrichment -CommandType Function).ScriptBlock

function Get-DanteInputProfile {
    if ($null -ne $script:DanteInputProfileCache) { return $script:DanteInputProfileCache }
    $helper = Join-Path $CodeRoot "patient_input_profile.py"
    if (-not (Test-Path -LiteralPath $helper -PathType Leaf)) {
        throw "Patient input profile detector is missing: $helper"
    }
    $python = Get-Python
    $probe = Invoke-NativeCapture $python @($helper, "detect", "--patient-root", $PatientRoot, "--json")
    if ($probe.ExitCode -ne 0) {
        throw "Patient input profile detection failed: $($probe.Output)"
    }
    try {
        $profile = $probe.Output | ConvertFrom-Json
    } catch {
        throw "Patient input profile detector returned invalid JSON."
    }
    if ([string]$profile.schema -ne "dante-patient-input-profile-v1") {
        throw "Unsupported patient input profile schema: $($profile.schema)"
    }
    $script:DanteInputProfileCache = $profile
    return $profile
}

function Test-IsArrayProfile {
    $profile = Get-DanteInputProfile
    return ([string]$profile.kind -eq "GENOTYPING_ARRAY_VCF")
}

function Get-ArraySourceFile {
    $profile = Get-DanteInputProfile
    if ([string]$profile.kind -ne "GENOTYPING_ARRAY_VCF") {
        throw "The active patient is not a supported genotyping-array VCF profile."
    }
    if ([string]$profile.assembly -ne "GRCh37") {
        throw "Genotyping-array adapter currently requires hg19/GRCh37 coordinates; detected: $($profile.assembly)"
    }
    $name = [string]$profile.source_file.name
    if (-not $name) { throw "Array profile has no source VCF identity." }
    $path = Join-Path $PatientInput $name
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "Array source VCF is missing: $path"
    }
    return (Get-Item -LiteralPath $path)
}

function Get-ArrayAdapterPath {
    $path = Join-Path $CodeRoot "array_vcf_tools.py"
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "Genotyping-array adapter is missing: $path"
    }
    return $path
}

function Get-ArrayRsidStatePath {
    return (Join-Path $Results "Dante_23andMe_style.array.state.json")
}

function Test-ArrayRsidExportReady {
    $source = Get-ArraySourceFile
    $adapter = Get-ArrayAdapterPath
    $export = Join-Path $Results "Dante_23andMe_style.txt.gz"
    $report = Join-Path $Qc "Dante_23andMe_style_report.txt"
    $statePath = Get-ArrayRsidStatePath
    foreach ($required in @($export, $report, $statePath)) {
        if (-not (Test-NonEmptyFile $required)) { return $false }
    }
    try {
        $state = Get-Content -LiteralPath $statePath -Raw -Encoding UTF8 | ConvertFrom-Json
        return (
            ([string]$state.schema -eq "dante-array-rsid-export-v1") -and
            ([string]$state.status -eq "complete") -and
            ([string]$state.source_sha256 -eq (Get-Sha256 $source.FullName)) -and
            ([string]$state.adapter_sha256 -eq (Get-Sha256 $adapter)) -and
            ([string]$state.output_sha256 -eq (Get-Sha256 $export))
        )
    } catch {
        return $false
    }
}

function Get-RsidExportMenuStatus {
    if (-not (Test-IsArrayProfile)) {
        return (& $script:DanteOriginalGetRsidExportMenuStatus)
    }
    $paths = @(
        (Join-Path $Results "Dante_23andMe_style.txt.gz"),
        (Join-Path $Results "Dante_23andMe_style.txt.gz.part"),
        (Get-ArrayRsidStatePath)
    )
    return (Get-ArtifactMenuStatus -Ready (Test-ArrayRsidExportReady) -Paths $paths)
}

function Run-RsidExport {
    if (-not (Test-IsArrayProfile)) {
        & $script:DanteOriginalRunRsidExport
        return
    }
    $source = Get-ArraySourceFile
    $adapter = Get-ArrayAdapterPath
    $python = Get-Python
    New-Item -ItemType Directory -Force -Path $Results, $Qc | Out-Null
    if ((Test-ArrayRsidExportReady) -and -not $ForceDownload) {
        Write-Host "Array genotype export is current for the same source VCF; reusing it." -ForegroundColor Green
        return
    }
    $export = Join-Path $Results "Dante_23andMe_style.txt.gz"
    $report = Join-Path $Qc "Dante_23andMe_style_report.txt"
    Write-Step "Create array-aware 23andMe-style genotype export"
    Write-Host "Measured 0/0 loci are preserved; absent loci remain UNKNOWN. No VEP lookup is required." -ForegroundColor Cyan
    Invoke-ObservedNative -Activity "Build Genotek/array genotype export" -Program $python -Arguments @(
        $adapter, "export-23andme", $source.FullName,
        "--output", $export,
        "--report", $report
    ) -WatchFile "$export.part"
    if (-not (Test-NonEmptyFile $export)) { throw "Array genotype export was not created." }
    $state = [ordered]@{
        schema = "dante-array-rsid-export-v1"
        status = "complete"
        profile = "GENOTYPING_ARRAY_VCF"
        source = $source.Name
        source_sha256 = Get-Sha256 $source.FullName
        adapter_sha256 = Get-Sha256 $adapter
        output_sha256 = Get-Sha256 $export
        whole_genome_negative_evidence = $false
    }
    Set-AtomicText (Get-ArrayRsidStatePath) ($state | ConvertTo-Json -Depth 5)
    Write-Host "Created: $export" -ForegroundColor Green
}

function Get-ArrayNormalizedStatePath {
    return (Join-Path $Work "Dante_array_normalized.state.json")
}

function Test-ArrayNormalizedLayersReady {
    if (-not (Test-FastaReady)) { return $false }
    $source = Get-ArraySourceFile
    $adapter = Get-ArrayAdapterPath
    $statePath = Get-ArrayNormalizedStatePath
    if (-not (Test-NonEmptyFile $statePath)) { return $false }
    foreach ($name in @("SNP", "INDEL", "RESIDUAL")) {
        $path = Join-Path $Work "Dante_$($name)_normalized.vcf.gz"
        if (-not (Test-NonEmptyFile $path) -or -not (Test-NonEmptyFile "$path.tbi")) { return $false }
    }
    try {
        $state = Get-Content -LiteralPath $statePath -Raw -Encoding UTF8 | ConvertFrom-Json
        $fasta = Join-Path $VepFastaData $FastaName
        return (
            ([string]$state.schema -eq "dante-array-normalized-v1") -and
            ([string]$state.status -eq "complete") -and
            ([string]$state.source_sha256 -eq (Get-Sha256 $source.FullName)) -and
            ([string]$state.adapter_sha256 -eq (Get-Sha256 $adapter)) -and
            ([string]$state.fasta_sha256 -eq (Get-Sha256 $fasta))
        )
    } catch {
        return $false
    }
}

function Test-NormalizedLayersReady {
    if (Test-IsArrayProfile) { return (Test-ArrayNormalizedLayersReady) }
    return (& $script:DanteOriginalTestNormalizedLayersReady)
}

function Get-NormalizedLayersMenuStatus {
    if (-not (Test-IsArrayProfile)) {
        return (& $script:DanteOriginalGetNormalizedLayersMenuStatus)
    }
    $paths = @(
        (Join-Path $Work "Dante_SNP_normalized.vcf.gz"),
        (Join-Path $Work "Dante_INDEL_normalized.vcf.gz"),
        (Join-Path $Work "Dante_RESIDUAL_normalized.vcf.gz"),
        (Get-ArrayNormalizedStatePath)
    )
    return (Get-ArtifactMenuStatus -Ready (Test-ArrayNormalizedLayersReady) -Paths $paths)
}

function Invoke-ArrayBcftoolsNormalize([string]$InputName, [string]$OutputName, [string]$Label, [string]$BcftoolsRuntimeImage) {
    $input = Join-Path $Work $InputName
    if (-not (Test-NonEmptyFile $input)) { throw "Prepared array $Label layer is missing: $input" }
    $output = Join-Path $Work $OutputName
    $partName = ($OutputName -replace '\.vcf\.gz$', '.part.vcf.gz')
    $part = Join-Path $Work $partName
    foreach ($stale in @($part, "$part.tbi")) {
        if (Test-Path -LiteralPath $stale) { Remove-Item -LiteralPath $stale -Force }
    }
    Write-Step "Normalize array $Label carried variants"
    Invoke-ObservedNative -Activity "bcftools normalize array $Label" -Program "docker" -Arguments @(
        "run", "--rm",
        "--mount", "type=bind,source=$PatientRoot,target=/patient",
        "--mount", "type=bind,source=$VepData,target=/data",
        "--entrypoint", "bcftools", $BcftoolsRuntimeImage,
        "norm", "--check-ref", "x",
        "--fasta-ref", "/data/fasta/$FastaName",
        "--multiallelics", "-any",
        "--threads", "$script:BcftoolsThreads",
        "--output-type", "z",
        "--output", "/patient/Dante_work/$partName",
        "/patient/Dante_work/$InputName"
    ) -WatchFile $part
    Invoke-ObservedNative -Activity "Index array $Label normalized layer" -Program "docker" -Arguments @(
        "run", "--rm",
        "--mount", "type=bind,source=$PatientRoot,target=/patient",
        "--entrypoint", "bcftools", $BcftoolsRuntimeImage,
        "index", "--tbi", "--force", "/patient/Dante_work/$partName"
    ) -WatchFile "$part.tbi"
    Move-Item -LiteralPath $part -Destination $output -Force
    Move-Item -LiteralPath "$part.tbi" -Destination "$output.tbi" -Force
}

function Prepare-NormalizedLayers {
    if (-not (Test-IsArrayProfile)) {
        & $script:DanteOriginalPrepareNormalizedLayers
        return
    }
    if (-not (Test-FastaReady)) { throw "GRCh37 FASTA is not installed. Select item 6 first." }
    Assert-Docker
    $source = Get-ArraySourceFile
    $adapter = Get-ArrayAdapterPath
    $python = Get-Python
    New-Item -ItemType Directory -Force -Path $Work, $Qc | Out-Null
    if ((Test-ArrayNormalizedLayersReady) -and -not $ForceDownload) {
        Write-Host "Array carried-variant layers are current; reusing them." -ForegroundColor Green
        return
    }

    $preparedSnp = Join-Path $Work "Dante_ARRAY_SNP_carried.vcf.gz"
    $preparedIndel = Join-Path $Work "Dante_ARRAY_INDEL_carried.vcf.gz"
    $preparedResidual = Join-Path $Work "Dante_ARRAY_RESIDUAL_carried.vcf.gz"
    $prepReport = Join-Path $Qc "Dante_ARRAY_preparation_report.txt"
    Write-Step "Prepare carried variants from genotyping-array VCF"
    Write-Host "0/0 and missing loci are excluded from variant annotation; source | phasing separators are not trusted." -ForegroundColor Cyan
    Invoke-ObservedNative -Activity "Prepare array carried-variant layers" -Program $python -Arguments @(
        $adapter, "prepare-layers", $source.FullName,
        "--snp-out", $preparedSnp,
        "--indel-out", $preparedIndel,
        "--residual-out", $preparedResidual,
        "--report", $prepReport
    ) -WatchFile "$preparedSnp.part"

    $bcftoolsRuntimeImage = Get-BcftoolsRuntimeImage
    Assert-BcftoolsCompatibility $bcftoolsRuntimeImage
    Invoke-ArrayBcftoolsNormalize "Dante_ARRAY_SNP_carried.vcf.gz" "Dante_SNP_normalized.vcf.gz" "SNP" $bcftoolsRuntimeImage
    Invoke-ArrayBcftoolsNormalize "Dante_ARRAY_INDEL_carried.vcf.gz" "Dante_INDEL_normalized.vcf.gz" "INDEL" $bcftoolsRuntimeImage
    Invoke-ArrayBcftoolsNormalize "Dante_ARRAY_RESIDUAL_carried.vcf.gz" "Dante_RESIDUAL_normalized.vcf.gz" "RESIDUAL" $bcftoolsRuntimeImage

    $fasta = Join-Path $VepFastaData $FastaName
    $state = [ordered]@{
        schema = "dante-array-normalized-v1"
        status = "complete"
        profile = "GENOTYPING_ARRAY_VCF"
        source = $source.Name
        source_sha256 = Get-Sha256 $source.FullName
        adapter_sha256 = Get-Sha256 $adapter
        fasta_sha256 = Get-Sha256 $fasta
        reference_mismatch_policy = "bcftools norm --check-ref x; mismatching REF records excluded, never strand-flipped"
        whole_genome_negative_evidence = $false
    }
    Set-AtomicText (Get-ArrayNormalizedStatePath) ($state | ConvertTo-Json -Depth 5)
    if (-not (Test-ArrayNormalizedLayersReady)) {
        throw "Array normalized layers were created but failed their READY validator."
    }
}

function Invoke-ArrayFullVepLayer([string]$Layer, [string]$InputName) {
    $input = Join-Path $Work $InputName
    if (-not (Test-NonEmptyFile $input)) { throw "Missing prepared array layer: $input" }
    $upper = $Layer.ToUpper()
    $outputName = "Dante_$($upper)_full.vep.vcf.gz"
    $output = Join-Path $FullVepResults $outputName
    $state = "$output.state.sha256"
    $cacheInfo = Join-Path $VepData "homo_sapiens\116_GRCh37\info.txt"
    $fingerprint = "VEP-116-array-GRCh37|$(Get-Sha256 $cacheInfo)|$(Get-Sha256 $input)"
    if ((Test-NonEmptyFile $output) -and (Test-NonEmptyFile $state) -and
        ((Get-Content -LiteralPath $state -Raw).Trim() -eq $fingerprint) -and -not $ForceDownload) {
        Write-Host "$upper array VEP output is current; reusing it."
    } else {
        $partName = "Dante_$($upper)_full.part.vep.vcf.gz"
        $part = Join-Path $FullVepResults $partName
        if (Test-Path -LiteralPath $part) { Remove-Item -LiteralPath $part -Force }
        Write-Step "Full VEP annotation: array $upper"
        Invoke-ObservedNative -Activity "Full VEP annotation: array $upper" -Program "docker" -Arguments @(
            "run", "--rm",
            "--mount", "type=bind,source=$PatientRoot,target=/patient",
            "--mount", "type=bind,source=$VepData,target=/data",
            $VepImage, "vep", "--cache", "--offline", "--cache_version", $VepCacheVersion,
            "--species", "homo_sapiens", "--assembly", "GRCh37", "--dir_cache", "/data",
            "--fasta", "/data/fasta/$FastaName", "--format", "vcf", "--vcf",
            "--everything", "--allele_number", "--show_ref_allele", "--uploaded_allele",
            "--transcript_version", "--check_ref", "--exclude_null_alleles",
            "--fork", "$Threads", "--buffer_size", "5000", "--no_stats",
            "--compress_output", "bgzip", "--force_overwrite",
            "--input_file", "/patient/Dante_work/$InputName",
            "--output_file", "/patient/Dante_results/full_vep/$partName",
            "--warning_file", "/patient/Dante_qc/Dante_$($upper)_full_vep.warnings.txt"
        ) -WatchFile $part
        if (-not (Test-GzipFile $part)) { throw "$upper array VEP output failed gzip validation." }
        Move-Item -LiteralPath $part -Destination $output -Force
        Set-AtomicText $state $fingerprint
    }

    $python = Get-Python
    $summary = Join-Path $FullVepResults "Dante_$($upper)_VEP_summary.tsv.gz"
    $rsidVcf = Join-Path $FullVepResults "Dante_$($upper)_rsid.vcf.gz"
    Invoke-ObservedNative -Activity "Postprocess array VEP layer: $upper" -Program $python -Arguments @(
        (Join-Path $CodeRoot "vep_postprocess.py"), $output,
        "--rsid-vcf", $rsidVcf, "--summary", $summary,
        "--report", (Join-Path $Qc "Dante_$($upper)_VEP_report.txt"),
        "--layer", $Layer
    ) -WatchFile "$summary.part"
}

function Run-FullVepAnnotation {
    if (-not (Test-IsArrayProfile)) {
        & $script:DanteOriginalRunFullVepAnnotation
        return
    }
    if (-not (Test-VepCacheReady)) { throw "VEP 116 GRCh37 cache is missing. Select item 1 first." }
    if (-not (Test-FastaReady)) { throw "GRCh37 FASTA is missing. Select item 6 first." }
    if (-not (Test-ArrayNormalizedLayersReady)) { throw "Array carried-variant layers are not READY. Select item 8 first." }
    Assert-Docker
    New-Item -ItemType Directory -Force -Path $Results, $FullVepResults, $Work, $Qc | Out-Null
    Invoke-ArrayFullVepLayer "snp" "Dante_SNP_normalized.vcf.gz"
    Invoke-ArrayFullVepLayer "indel" "Dante_INDEL_normalized.vcf.gz"
    Invoke-ArrayFullVepLayer "residual" "Dante_RESIDUAL_normalized.vcf.gz"
    Write-Host "Array VEP layers are ready in $FullVepResults" -ForegroundColor Green
}

function Run-Stage1 {
    if (Test-IsArrayProfile) {
        Run-RsidExport
        return
    }
    & $script:DanteOriginalRunStage1
}

function Ensure-ArrayStructuralSentinels {
    if (-not (Test-IsArrayProfile)) { return }
    New-Item -ItemType Directory -Force -Path $Results, $Qc | Out-Null
    $header = "AnnotSV_ID`tSV_chrom`tSV_start`tSV_end`tSV_type"
    foreach ($kind in @("CNV", "SV")) {
        $output = Join-Path $Results "Dante_$($kind)_AnnotSV.tsv"
        $report = Join-Path $Qc "Dante_$($kind)_AnnotSV_input_report.txt"
        Set-AtomicText $output $header
        Set-AtomicText $report @"
DANTE ARRAY STRUCTURAL INPUT SENTINEL
Selected rows: 0
Scope: NOT APPLICABLE. The active genotyping-array VCF does not provide validated $kind calls.
This file exists only so Stage 24 can represent zero structural inputs without claiming they were tested.
"@
    }
}

function Get-AnnotSVOutputMenuStatus {
    if (Test-IsArrayProfile) { return "READY" }
    return (& $script:DanteOriginalGetAnnotSVOutputMenuStatus)
}

function Run-KnowledgeEnrichment {
    if (-not (Test-IsArrayProfile)) {
        & $script:DanteOriginalRunKnowledgeEnrichment
        return
    }
    if (-not (Test-UnifiedCatalogReady)) {
        throw "Array Stage 24 requires the current small-variant catalog. Select item 10/12 first."
    }
    Ensure-ArrayStructuralSentinels
    $oldKind = [Environment]::GetEnvironmentVariable("DANTE_INPUT_PROFILE_KIND", "Process")
    try {
        [Environment]::SetEnvironmentVariable("DANTE_INPUT_PROFILE_KIND", "GENOTYPING_ARRAY_VCF", "Process")
        & $script:DanteOriginalRunKnowledgeEnrichment
    } finally {
        [Environment]::SetEnvironmentVariable("DANTE_INPUT_PROFILE_KIND", $oldKind, "Process")
    }
    $python = Get-Python
    Invoke-Native $python @(
        (Join-Path $CodeRoot "patient_input_profile.py"),
        "annotate-stage24", "--patient-root", $PatientRoot
    )
    if (-not (Test-KnowledgeEnrichmentReady)) {
        throw "Array Stage 24 completed but input-scope annotation failed the READY validator."
    }
}

function Show-AutomaticCheck {
    if (-not (Test-IsArrayProfile)) {
        & $script:DanteOriginalShowAutomaticCheck
        return
    }
    $profile = Get-DanteInputProfile
    $source = Get-ArraySourceFile
    Write-Step "Automatic startup check"
    Write-CheckStatus "Code directory" (Test-Path -LiteralPath $CodeRoot -PathType Container) $CodeRoot
    Write-CheckStatus "Active patient" $true "$Patient | $PatientRoot"
    Write-CheckStatus "Input profile" $true "$($profile.kind) | $($profile.source) | $($profile.assembly)"
    Write-CheckStatus "Array VCF" $true "$($source.Name) | $([Math]::Round($source.Length / 1MB, 2)) MiB"
    Write-Host "Scope: $($profile.scope)" -ForegroundColor Yellow
    Write-Host "Not available from this source: WGS coverage, FASTQ evidence, CNV/SV calls, telomere estimate, BAM browser." -ForegroundColor DarkGray
    $dockerCli = [bool](Get-Command docker -ErrorAction SilentlyContinue)
    Write-CheckStatus "Docker CLI" $dockerCli $(if ($dockerCli) { "available" } else { "not found" })
    $pythonCommand = Get-Command py -ErrorAction SilentlyContinue
    if (-not $pythonCommand) { $pythonCommand = Get-Command python -ErrorAction SilentlyContinue }
    Write-CheckStatus "Python 3" ([bool]$pythonCommand) $(if ($pythonCommand) { $pythonCommand.Source } else { "not found" })
}

function Select-Mode {
    if (-not (Test-IsArrayProfile)) {
        return (& $script:DanteOriginalSelectMode)
    }

    Show-AutomaticCheck
    $profile = Get-DanteInputProfile
    Write-Host ""
    Write-Host ("Dante genomic pipeline v{0}" -f (Get-DanteDisplayVersion)) -ForegroundColor Cyan
    Write-Host ("Patient: {0}" -f $Patient) -ForegroundColor Green
    Write-Host ("Workspace: {0}" -f $PatientRoot) -ForegroundColor DarkGray
    Write-Host ("Input: {0} | {1} | {2}" -f $profile.kind, $profile.source, $profile.assembly) -ForegroundColor Cyan
    Write-Host "Scope: interrogated array loci only; absence outside the panel is UNKNOWN." -ForegroundColor Yellow
    Write-Host ""

    $rsidStatus = Get-RsidExportMenuStatus
    $vepCacheStatus = Get-VepCacheMenuStatus
    $fastaStatus = if (Test-FastaReady) { "READY" } else { "MISSING" }
    $clinVarStatus = if (Test-ClinVarReady) { "READY" } else { "MISSING" }
    $normalizedStatus = Get-NormalizedLayersMenuStatus
    $fullVepStatus = Get-FullVepMenuStatus
    $catalogStatus = Get-UnifiedCatalogMenuStatus
    $shortlistStatus = Get-ShortlistMenuStatus
    $stage2Status = Get-CompleteStageMenuStatus -Statuses @(
        $vepCacheStatus, $fastaStatus, $clinVarStatus, $normalizedStatus,
        $fullVepStatus, $catalogStatus, $shortlistStatus
    )
    $openReferenceStatus = Get-OpenReferenceMenuStatus
    $researchReferenceStatus = Get-ResearchReferenceMenuStatus
    $knowledgeEnrichmentStatus = Get-KnowledgeEnrichmentMenuStatus
    $commercialReportStatus = Get-GenomeReportMenuStatus "commercial"
    $researchReportStatus = Get-GenomeReportMenuStatus "research"
    $allLocalStatus = Get-CompleteStageMenuStatus -Statuses @($rsidStatus, $stage2Status)

    Write-Host "ARRAY GENOTYPE EXPORT" -ForegroundColor DarkCyan
    Write-MenuItem "2" "Create measured-loci 23andMe-style genotype export" $rsidStatus
    Write-Host ""
    Write-Host "SMALL-VARIANT ANNOTATION - ARRAY ALT-CARRIERS ONLY" -ForegroundColor DarkCyan
    Write-MenuItem "1" "Install/resume VEP 116 GRCh37 cache" $vepCacheStatus
    Write-MenuItem "6" "Install GRCh37 reference FASTA" $fastaStatus
    Write-MenuItem "7" "Install/update current ClinVar GRCh37" $clinVarStatus
    Write-MenuItem "8" "Prepare normalized array ALT-carrier SNP/INDEL layers" $normalizedStatus
    Write-MenuItem "9" "Run full VEP annotation of measured ALT-carriers" $fullVepStatus
    Write-MenuItem "10" "Build unified annotated measured-variant catalog" $catalogStatus
    Write-MenuItem "11" "Build deterministic shortlist within measured loci" $shortlistStatus
    Write-MenuItem "12" "Run complete array small-variant annotation" $stage2Status
    Write-Host ""
    Write-Host "REFERENCE KNOWLEDGE + REPORTS" -ForegroundColor DarkCyan
    Write-MenuItem "22" "Install/update open reference knowledge bundle" $openReferenceStatus
    Write-MenuItem "23" "Install research only knowledge bundle (NO COMMERCIAL USE ALLOWED)" $researchReferenceStatus
    Write-MenuItem "24" "Enrich measured array variants with reference knowledge (NO OPENAI)" $knowledgeEnrichmentStatus
    Write-MenuItem "25" "Generate CLIENT PDF report (ARRAY SCOPE / OPEN / NO AI)" $commercialReportStatus
    Write-MenuItem "26" "Generate RESEARCH ONLY PDF report (ARRAY SCOPE / NON-COMMERCIAL)" $researchReportStatus
    Write-Host ""
    Write-Host "AUTOMATION" -ForegroundColor DarkCyan
    Write-MenuItem "21" "Run all applicable array processing (local, no API)" $allLocalStatus
    Write-MenuItem "0" "Exit"

    $choice = Read-Host "Select"
    switch ($choice) {
        "1" { return "setup-vep-cache" }
        "2" { return "rsid" }
        "6" { return "setup-fasta" }
        "7" { return "setup-clinvar" }
        "8" { return "prepare-small" }
        "9" { return "vep-full" }
        "10" { return "catalog" }
        "11" { return "shortlist" }
        "12" { return "stage2" }
        "21" { return "local-all" }
        "22" { return "open-reference" }
        "23" { return "research-reference" }
        "24" { return "knowledge-enrichment" }
        "25" { return "commercial-report" }
        "26" { return "research-report" }
        "0" { return "exit" }
        default { throw "This action is not available for the active genotyping-array profile: $choice" }
    }
}
