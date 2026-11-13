# Array normalization QC overlay for Dante v9.
# Loaded after array_profile_runtime.ps1. It does not change genotype semantics;
# it makes bcftools REF-mismatch exclusions explicit and auditable.

$script:DanteArrayQcOriginalPrepareNormalizedLayers = (Get-Command Prepare-NormalizedLayers -CommandType Function).ScriptBlock
$script:DanteArrayQcOriginalGetNormalizedLayersMenuStatus = (Get-Command Get-NormalizedLayersMenuStatus -CommandType Function).ScriptBlock

function Get-ArrayNormalizationQcHelperPath {
    $path = Join-Path $CodeRoot "array_normalization_qc.py"
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "Array normalization QC helper is missing: $path"
    }
    return $path
}

function Test-ArrayNormalizationQcReady {
    if (-not (Test-IsArrayProfile)) { return $true }
    $statePath = Get-ArrayNormalizedStatePath
    if (-not (Test-NonEmptyFile $statePath)) { return $false }
    try {
        $state = Get-Content -LiteralPath $statePath -Raw -Encoding UTF8 | ConvertFrom-Json
        if (-not ($state.PSObject.Properties.Name -contains "normalization_qc")) { return $false }
        $qc = $state.normalization_qc
        $helper = Get-ArrayNormalizationQcHelperPath
        return (
            ([string]$qc.schema -eq "dante-array-normalization-qc-v1") -and
            ([string]$qc.check_ref_policy -eq "x") -and
            ([string]$qc.helper_sha256 -eq (Get-Sha256 $helper))
        )
    } catch {
        return $false
    }
}

function Get-NormalizedLayersMenuStatus {
    if (-not (Test-IsArrayProfile)) {
        return (& $script:DanteArrayQcOriginalGetNormalizedLayersMenuStatus)
    }
    $baseStatus = & $script:DanteArrayQcOriginalGetNormalizedLayersMenuStatus
    if (($baseStatus -match '^READY') -and -not (Test-ArrayNormalizationQcReady)) {
        return "READY: QC PENDING"
    }
    return $baseStatus
}

function Read-ArrayBcftoolsNormQc([string]$LogPath, [string]$Label, [string]$OutputPath, [string]$BcftoolsRuntimeImage) {
    $helper = Get-ArrayNormalizationQcHelperPath
    $python = Get-Python
    $probe = Invoke-NativeCapture $python @($helper, "parse-log", $LogPath, "--json")
    if ($probe.ExitCode -ne 0) {
        throw "Could not parse bcftools normalization QC for ${Label}: $($probe.Output)"
    }
    try {
        $summary = $probe.Output | ConvertFrom-Json
    } catch {
        throw "bcftools normalization QC parser returned invalid JSON for $Label."
    }

    $outputName = [System.IO.Path]::GetFileName($OutputPath)
    $countProbe = Invoke-NativeCapture "docker" @(
        "run", "--rm",
        "--mount", "type=bind,source=$PatientRoot,target=/patient",
        "--entrypoint", "bcftools", $BcftoolsRuntimeImage,
        "index", "--nrecords", "/patient/Dante_work/$outputName"
    )
    if ($countProbe.ExitCode -ne 0) {
        throw "Could not count normalized $Label records: $($countProbe.Output)"
    }
    $countMatch = [regex]::Match([string]$countProbe.Output, '(?m)^\s*(\d+)\s*$')
    if (-not $countMatch.Success) {
        throw "Unexpected bcftools index --nrecords output for ${Label}: $($countProbe.Output)"
    }

    return [pscustomobject][ordered]@{
        layer = $Label.ToUpper()
        input_records = [int64]$summary.input_records
        output_records = [int64]$countMatch.Groups[1].Value
        split_records = [int64]$summary.split_records
        joined_records = [int64]$summary.joined_records
        realigned_records = [int64]$summary.realigned_records
        reference_mismatch_records_excluded = [int64]$summary.reference_mismatch_records_excluded
        reference_mismatch_counter_source = [string]$summary.reference_mismatch_counter_source
        duplicate_records_removed = [int64]$summary.duplicate_records_removed
        skipped_records = [int64]$summary.skipped_records
        summary_format = [string]$summary.summary_format
    }
}

function Invoke-ArrayBcftoolsNormalizeWithQc(
    [string]$InputName,
    [string]$OutputName,
    [string]$Label,
    [string]$BcftoolsRuntimeImage,
    [switch]$ProbeOnly
) {
    $input = Join-Path $Work $InputName
    if (-not (Test-NonEmptyFile $input)) { throw "Prepared array $Label layer is missing: $input" }
    $output = Join-Path $Work $OutputName
    if ($ProbeOnly -and -not (Test-NonEmptyFile $output)) {
        throw "Normalized array $Label layer is missing for QC backfill: $output"
    }

    $partName = ($OutputName -replace '\.vcf\.gz$', '.part.vcf.gz')
    $part = Join-Path $Work $partName
    $logName = "Dante_ARRAY_$($Label.ToUpper())_bcftools_norm.log"
    $workLog = Join-Path $Work $logName
    $qcLog = Join-Path $Qc $logName
    foreach ($stale in @($workLog)) {
        if (Test-Path -LiteralPath $stale) { Remove-Item -LiteralPath $stale -Force }
    }
    if (-not $ProbeOnly) {
        foreach ($stale in @($part, "$part.tbi")) {
            if (Test-Path -LiteralPath $stale) { Remove-Item -LiteralPath $stale -Force }
        }
    }

    $containerOutput = if ($ProbeOnly) { "/dev/null" } else { "/patient/Dante_work/$partName" }
    $outputType = if ($ProbeOnly) { "u" } else { "z" }
    $containerLog = "/patient/Dante_work/$logName"
    $activity = if ($ProbeOnly) { "Audit bcftools REF exclusions for array $Label" } else { "bcftools normalize array $Label" }
    $step = if ($ProbeOnly) { "Backfill normalization QC for array $Label" } else { "Normalize array $Label carried variants" }
    Write-Step $step

    $command = "bcftools norm --check-ref x --fasta-ref '/data/fasta/$FastaName' --multiallelics -any --threads $script:BcftoolsThreads --output-type $outputType --output '$containerOutput' '/patient/Dante_work/$InputName' 2> '$containerLog'; status=`$?; cat '$containerLog' >&2; exit `$status"
    $watch = if ($ProbeOnly) { $workLog } else { $part }
    Invoke-ObservedNative -Activity $activity -Program "docker" -Arguments @(
        "run", "--rm",
        "--mount", "type=bind,source=$PatientRoot,target=/patient",
        "--mount", "type=bind,source=$VepData,target=/data",
        "--entrypoint", "/bin/sh", $BcftoolsRuntimeImage,
        "-c", $command
    ) -WatchFile $watch

    if (-not (Test-NonEmptyFile $workLog)) {
        throw "bcftools normalization QC log was not created for $Label."
    }

    if (-not $ProbeOnly) {
        Invoke-ObservedNative -Activity "Index array $Label normalized layer" -Program "docker" -Arguments @(
            "run", "--rm",
            "--mount", "type=bind,source=$PatientRoot,target=/patient",
            "--entrypoint", "bcftools", $BcftoolsRuntimeImage,
            "index", "--tbi", "--force", "/patient/Dante_work/$partName"
        ) -WatchFile "$part.tbi"
        Move-Item -LiteralPath $part -Destination $output -Force
        Move-Item -LiteralPath "$part.tbi" -Destination "$output.tbi" -Force
    }

    New-Item -ItemType Directory -Force -Path $Qc | Out-Null
    Copy-Item -LiteralPath $workLog -Destination $qcLog -Force
    $result = Read-ArrayBcftoolsNormQc $qcLog $Label $output $BcftoolsRuntimeImage
    Remove-Item -LiteralPath $workLog -Force -ErrorAction SilentlyContinue
    return $result
}

function New-ArrayNormalizationQcPayload([object[]]$Layers, [string]$BcftoolsRuntimeImage) {
    $helper = Get-ArrayNormalizationQcHelperPath
    $mismatchTotal = [int64](($Layers | Measure-Object -Property reference_mismatch_records_excluded -Sum).Sum)
    $duplicateRemovedTotal = [int64](($Layers | Measure-Object -Property duplicate_records_removed -Sum).Sum)
    $skippedTotal = [int64](($Layers | Measure-Object -Property skipped_records -Sum).Sum)
    $inputTotal = [int64](($Layers | Measure-Object -Property input_records -Sum).Sum)
    $outputTotal = [int64](($Layers | Measure-Object -Property output_records -Sum).Sum)
    $layerMap = [ordered]@{}
    foreach ($layer in $Layers) {
        $layerMap[[string]$layer.layer] = [ordered]@{
            input_records = [int64]$layer.input_records
            output_records = [int64]$layer.output_records
            split_records = [int64]$layer.split_records
            joined_records = [int64]$layer.joined_records
            realigned_records = [int64]$layer.realigned_records
            reference_mismatch_records_excluded = [int64]$layer.reference_mismatch_records_excluded
            reference_mismatch_counter_source = [string]$layer.reference_mismatch_counter_source
            duplicate_records_removed = [int64]$layer.duplicate_records_removed
            skipped_records = [int64]$layer.skipped_records
            summary_format = [string]$layer.summary_format
        }
    }
    return [ordered]@{
        schema = "dante-array-normalization-qc-v1"
        helper_sha256 = Get-Sha256 $helper
        bcftools_image = $BcftoolsRuntimeImage
        check_ref_policy = "x"
        interpretation = "bcftools norm --check-ref x excludes REF/reference mismatches. In bcftools 1.24 the explicit mismatch_removed counter is authoritative; legacy supported summaries use skipped. Variants are never strand-flipped."
        input_records_total = $inputTotal
        output_records_total = $outputTotal
        reference_mismatch_records_excluded_total = $mismatchTotal
        duplicate_records_removed_total = $duplicateRemovedTotal
        skipped_records_total = $skippedTotal
        layers = $layerMap
    }
}

function Add-ArrayNormalizationQcToReport([string]$ReportPath, [object]$Payload) {
    $lines = @(
        "",
        "BCFTOOLS NORMALIZATION QC",
        ("=" * 78),
        "Policy: bcftools norm --check-ref x; REF/reference mismatches are excluded, never strand-flipped.",
        "For bcftools 1.24, mismatch_removed is the authoritative REF-mismatch exclusion counter; dup_removed and skipped are reported separately.",
        "",
        "normalization_input_records_total: $($Payload.input_records_total)",
        "normalization_output_records_total: $($Payload.output_records_total)",
        "reference_mismatch_records_excluded_total: $($Payload.reference_mismatch_records_excluded_total)",
        "duplicate_records_removed_total: $($Payload.duplicate_records_removed_total)",
        "skipped_records_total: $($Payload.skipped_records_total)"
    )
    foreach ($name in @("SNP", "INDEL", "RESIDUAL")) {
        $layer = $Payload.layers.$name
        $prefix = $name.ToLower()
        $lines += ""
        $lines += "$($prefix)_input_records: $($layer.input_records)"
        $lines += "$($prefix)_output_records: $($layer.output_records)"
        $lines += "$($prefix)_split_records: $($layer.split_records)"
        $lines += "$($prefix)_joined_records: $($layer.joined_records)"
        $lines += "$($prefix)_realigned_records: $($layer.realigned_records)"
        $lines += "$($prefix)_reference_mismatch_records_excluded: $($layer.reference_mismatch_records_excluded)"
        $lines += "$($prefix)_reference_mismatch_counter_source: $($layer.reference_mismatch_counter_source)"
        $lines += "$($prefix)_duplicate_records_removed: $($layer.duplicate_records_removed)"
        $lines += "$($prefix)_skipped_records: $($layer.skipped_records)"
    }
    Add-Content -LiteralPath $ReportPath -Value $lines -Encoding UTF8
}

function Set-ArrayNormalizationQcState([string]$StatePath, [object]$Payload) {
    if (-not (Test-NonEmptyFile $StatePath)) { throw "Array normalized state is missing: $StatePath" }
    $state = Get-Content -LiteralPath $StatePath -Raw -Encoding UTF8 | ConvertFrom-Json
    $state | Add-Member -NotePropertyName normalization_qc -NotePropertyValue $Payload -Force
    Set-AtomicText $StatePath ($state | ConvertTo-Json -Depth 10)
}

function Prepare-NormalizedLayers {
    if (-not (Test-IsArrayProfile)) {
        & $script:DanteArrayQcOriginalPrepareNormalizedLayers
        return
    }
    if (-not (Test-FastaReady)) { throw "GRCh37 FASTA is not installed. Select item 6 first." }
    Assert-Docker
    $source = Get-ArraySourceFile
    $adapter = Get-ArrayAdapterPath
    $python = Get-Python
    New-Item -ItemType Directory -Force -Path $Work, $Qc | Out-Null

    $preparedSnp = Join-Path $Work "Dante_ARRAY_SNP_carried.vcf.gz"
    $preparedIndel = Join-Path $Work "Dante_ARRAY_INDEL_carried.vcf.gz"
    $preparedResidual = Join-Path $Work "Dante_ARRAY_RESIDUAL_carried.vcf.gz"
    $prepReport = Join-Path $Qc "Dante_ARRAY_preparation_report.txt"
    $baseReady = Test-ArrayNormalizedLayersReady

    if ($baseReady -and (Test-ArrayNormalizationQcReady) -and -not $ForceDownload) {
        Write-Host "Array carried-variant layers and normalization QC are current; reusing them." -ForegroundColor Green
        return
    }

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

    if ($baseReady -and -not $ForceDownload) {
        Write-Host "Normalized VCFs are already current; running QC-only replay without changing them." -ForegroundColor Cyan
        $layers = @(
            (Invoke-ArrayBcftoolsNormalizeWithQc "Dante_ARRAY_SNP_carried.vcf.gz" "Dante_SNP_normalized.vcf.gz" "SNP" $bcftoolsRuntimeImage -ProbeOnly),
            (Invoke-ArrayBcftoolsNormalizeWithQc "Dante_ARRAY_INDEL_carried.vcf.gz" "Dante_INDEL_normalized.vcf.gz" "INDEL" $bcftoolsRuntimeImage -ProbeOnly),
            (Invoke-ArrayBcftoolsNormalizeWithQc "Dante_ARRAY_RESIDUAL_carried.vcf.gz" "Dante_RESIDUAL_normalized.vcf.gz" "RESIDUAL" $bcftoolsRuntimeImage -ProbeOnly)
        )
        $payload = New-ArrayNormalizationQcPayload $layers $bcftoolsRuntimeImage
        Add-ArrayNormalizationQcToReport $prepReport $payload
        Set-ArrayNormalizationQcState (Get-ArrayNormalizedStatePath) $payload
        Write-Host "Normalization QC backfilled without changing normalized variant outputs." -ForegroundColor Green
        Write-Host "REF-mismatch exclusions: $($payload.reference_mismatch_records_excluded_total)" -ForegroundColor Green
        return
    }

    $layers = @(
        (Invoke-ArrayBcftoolsNormalizeWithQc "Dante_ARRAY_SNP_carried.vcf.gz" "Dante_SNP_normalized.vcf.gz" "SNP" $bcftoolsRuntimeImage),
        (Invoke-ArrayBcftoolsNormalizeWithQc "Dante_ARRAY_INDEL_carried.vcf.gz" "Dante_INDEL_normalized.vcf.gz" "INDEL" $bcftoolsRuntimeImage),
        (Invoke-ArrayBcftoolsNormalizeWithQc "Dante_ARRAY_RESIDUAL_carried.vcf.gz" "Dante_RESIDUAL_normalized.vcf.gz" "RESIDUAL" $bcftoolsRuntimeImage)
    )
    $payload = New-ArrayNormalizationQcPayload $layers $bcftoolsRuntimeImage
    Add-ArrayNormalizationQcToReport $prepReport $payload

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
        normalization_qc = $payload
        whole_genome_negative_evidence = $false
    }
    Set-AtomicText (Get-ArrayNormalizedStatePath) ($state | ConvertTo-Json -Depth 10)
    if (-not (Test-ArrayNormalizedLayersReady)) {
        throw "Array normalized layers were created but failed their READY validator."
    }
    if (-not (Test-ArrayNormalizationQcReady)) {
        throw "Array normalized layers were created but normalization QC failed its READY validator."
    }
    Write-Host "REF-mismatch exclusions: $($payload.reference_mismatch_records_excluded_total)" -ForegroundColor Green
}
