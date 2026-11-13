# 23andMe-style genotype-table patient runtime for Dante v9.
# Loaded after the WGS, array-VCF and array-QC overlays.  Every override below
# delegates unchanged for non-TXT profiles, so this is a third independent path.

$script:DanteGenotypeTxtProfileCache = $null
$script:DanteGenotypeTxtOriginalGetDanteInputProfile = (Get-Command Get-DanteInputProfile -CommandType Function).ScriptBlock
$script:DanteGenotypeTxtOriginalShowAutomaticCheck = (Get-Command Show-AutomaticCheck -CommandType Function).ScriptBlock
$script:DanteGenotypeTxtOriginalSelectMode = (Get-Command Select-Mode -CommandType Function).ScriptBlock
$script:DanteGenotypeTxtOriginalRunRsidExport = (Get-Command Run-RsidExport -CommandType Function).ScriptBlock
$script:DanteGenotypeTxtOriginalGetRsidExportMenuStatus = (Get-Command Get-RsidExportMenuStatus -CommandType Function).ScriptBlock
$script:DanteGenotypeTxtOriginalPrepareNormalizedLayers = (Get-Command Prepare-NormalizedLayers -CommandType Function).ScriptBlock
$script:DanteGenotypeTxtOriginalTestNormalizedLayersReady = (Get-Command Test-NormalizedLayersReady -CommandType Function).ScriptBlock
$script:DanteGenotypeTxtOriginalGetNormalizedLayersMenuStatus = (Get-Command Get-NormalizedLayersMenuStatus -CommandType Function).ScriptBlock
$script:DanteGenotypeTxtOriginalRunFullVepAnnotation = (Get-Command Run-FullVepAnnotation -CommandType Function).ScriptBlock
$script:DanteGenotypeTxtOriginalRunStage1 = (Get-Command Run-Stage1 -CommandType Function).ScriptBlock
$script:DanteGenotypeTxtOriginalGetAnnotSVOutputMenuStatus = (Get-Command Get-AnnotSVOutputMenuStatus -CommandType Function).ScriptBlock
$script:DanteGenotypeTxtOriginalRunKnowledgeEnrichment = (Get-Command Run-KnowledgeEnrichment -CommandType Function).ScriptBlock

function Get-GenotypeTxtProfileHelperPath {
    return (Join-Path $CodeRoot "genotype_txt_profile.py")
}

function Get-DanteInputProfile {
    if ($null -ne $script:DanteGenotypeTxtProfileCache) { return $script:DanteGenotypeTxtProfileCache }

    $baseProfile = & $script:DanteGenotypeTxtOriginalGetDanteInputProfile
    if ([string]$baseProfile.kind -ne "UNKNOWN") {
        $script:DanteGenotypeTxtProfileCache = $baseProfile
        return $baseProfile
    }

    $helper = Get-GenotypeTxtProfileHelperPath
    if (-not (Test-Path -LiteralPath $helper -PathType Leaf)) {
        return $baseProfile
    }
    $python = Get-Python
    $probe = Invoke-NativeCapture $python @(
        $helper, "detect", "--patient-root", $PatientRoot, "--json", "--save"
    )
    if ($probe.ExitCode -ne 0) {
        throw "Genotype-table TXT profile detection failed: $($probe.Output)"
    }
    try {
        $profile = $probe.Output | ConvertFrom-Json
    } catch {
        throw "Genotype-table TXT detector returned invalid JSON."
    }
    if ([string]$profile.kind -eq "GENOTYPE_TABLE_TXT") {
        if ([string]$profile.schema -ne "dante-patient-input-profile-v1") {
            throw "Unsupported genotype-table TXT profile schema: $($profile.schema)"
        }
        $script:DanteGenotypeTxtProfileCache = $profile
        return $profile
    }
    return $baseProfile
}

function Test-IsGenotypeTxtProfile {
    $profile = Get-DanteInputProfile
    return ([string]$profile.kind -eq "GENOTYPE_TABLE_TXT")
}

function Get-GenotypeTxtSourceFile {
    $profile = Get-DanteInputProfile
    if ([string]$profile.kind -ne "GENOTYPE_TABLE_TXT") {
        throw "The active patient is not a supported genotype-table TXT profile."
    }
    if ([string]$profile.assembly -ne "GRCh37") {
        throw "Genotype-table TXT adapter currently requires GRCh37/build 37 coordinates; detected: $($profile.assembly)"
    }
    $name = [string]$profile.source_file.name
    if (-not $name) { throw "Genotype-table TXT profile has no source-file identity." }
    $path = Join-Path $PatientInput $name
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "Genotype-table TXT source is missing: $path"
    }
    return (Get-Item -LiteralPath $path)
}

function Get-GenotypeTxtAdapterPath {
    $path = Join-Path $CodeRoot "genotype_txt_tools.py"
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "Genotype-table TXT adapter is missing: $path"
    }
    return $path
}

function Get-GenotypeTxtRsidStatePath {
    return (Join-Path $Results "Dante_23andMe_style.txt.state.json")
}

function Test-GenotypeTxtRsidExportReady {
    $source = Get-GenotypeTxtSourceFile
    $adapter = Get-GenotypeTxtAdapterPath
    $export = Join-Path $Results "Dante_23andMe_style.txt.gz"
    $report = Join-Path $Qc "Dante_23andMe_style_report.txt"
    $statePath = Get-GenotypeTxtRsidStatePath
    foreach ($required in @($export, $report, $statePath)) {
        if (-not (Test-NonEmptyFile $required)) { return $false }
    }
    try {
        $state = Get-Content -LiteralPath $statePath -Raw -Encoding UTF8 | ConvertFrom-Json
        return (
            ([string]$state.schema -eq "dante-genotype-txt-export-v1") -and
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
    if (-not (Test-IsGenotypeTxtProfile)) {
        return (& $script:DanteGenotypeTxtOriginalGetRsidExportMenuStatus)
    }
    $paths = @(
        (Join-Path $Results "Dante_23andMe_style.txt.gz"),
        (Join-Path $Results "Dante_23andMe_style.txt.gz.part"),
        (Get-GenotypeTxtRsidStatePath)
    )
    return (Get-ArtifactMenuStatus -Ready (Test-GenotypeTxtRsidExportReady) -Paths $paths)
}

function Run-RsidExport {
    if (-not (Test-IsGenotypeTxtProfile)) {
        & $script:DanteGenotypeTxtOriginalRunRsidExport
        return
    }
    $source = Get-GenotypeTxtSourceFile
    $adapter = Get-GenotypeTxtAdapterPath
    $python = Get-Python
    New-Item -ItemType Directory -Force -Path $Results, $Qc | Out-Null
    if ((Test-GenotypeTxtRsidExportReady) -and -not $ForceDownload) {
        Write-Host "Genotype-table direct export is current; reusing it." -ForegroundColor Green
        return
    }
    $export = Join-Path $Results "Dante_23andMe_style.txt.gz"
    $report = Join-Path $Qc "Dante_23andMe_style_report.txt"
    Write-Step "Create measured-loci genotype-table export"
    Write-Host "All valid source rows are preserved, including -- no-calls and symbolic I/D genotypes." -ForegroundColor Cyan
    Invoke-ObservedNative -Activity "Build genotype-table direct export" -Program $python -Arguments @(
        $adapter, "export-23andme", $source.FullName,
        "--output", $export,
        "--report", $report
    ) -WatchFile "$export.part"
    if (-not (Test-NonEmptyFile $export)) { throw "Genotype-table direct export was not created." }
    $state = [ordered]@{
        schema = "dante-genotype-txt-export-v1"
        status = "complete"
        profile = "GENOTYPE_TABLE_TXT"
        source = $source.Name
        source_sha256 = Get-Sha256 $source.FullName
        adapter_sha256 = Get-Sha256 $adapter
        output_sha256 = Get-Sha256 $export
        whole_genome_negative_evidence = $false
        sequence_resolved_indels = $false
    }
    Set-AtomicText (Get-GenotypeTxtRsidStatePath) ($state | ConvertTo-Json -Depth 5)
    Write-Host "Created: $export" -ForegroundColor Green
}

function Get-GenotypeTxtFastaPath {
    return (Join-Path $VepFastaData $FastaName)
}

function Ensure-GenotypeTxtFastaIndex {
    if (-not (Test-FastaReady)) { throw "GRCh37 FASTA is not installed. Select item 6 first." }
    $fasta = Get-GenotypeTxtFastaPath
    $fai = "$fasta.fai"
    if ((Test-NonEmptyFile $fai) -and ((Get-Item -LiteralPath $fai).LastWriteTimeUtc -ge (Get-Item -LiteralPath $fasta).LastWriteTimeUtc)) {
        return $fai
    }
    Assert-Docker
    $samtoolsRuntimeImage = if (Get-Command Get-SamtoolsRuntimeImage -ErrorAction SilentlyContinue) {
        Get-SamtoolsRuntimeImage
    } else {
        $SamtoolsImage
    }
    Write-Step "Index GRCh37 FASTA for genotype-table REF resolution"
    Invoke-ObservedNative -Activity "samtools faidx GRCh37 reference" -Program "docker" -Arguments @(
        "run", "--rm",
        "--mount", "type=bind,source=$VepFastaData,target=/data",
        "--entrypoint", "samtools", $samtoolsRuntimeImage,
        "faidx", "/data/$FastaName"
    ) -WatchFile $fai
    if (-not (Test-NonEmptyFile $fai)) { throw "samtools did not create the GRCh37 FASTA index: $fai" }
    return $fai
}

function Get-GenotypeTxtNormalizedStatePath {
    return (Join-Path $Work "Dante_genotype_txt_normalized.state.json")
}

function Test-GenotypeTxtNormalizedLayersReady {
    if (-not (Test-FastaReady)) { return $false }
    $source = Get-GenotypeTxtSourceFile
    $adapter = Get-GenotypeTxtAdapterPath
    $fasta = Get-GenotypeTxtFastaPath
    $fai = "$fasta.fai"
    $statePath = Get-GenotypeTxtNormalizedStatePath
    if (-not (Test-NonEmptyFile $fai) -or -not (Test-NonEmptyFile $statePath)) { return $false }
    foreach ($name in @("SNP", "INDEL", "RESIDUAL")) {
        $path = Join-Path $Work "Dante_$($name)_normalized.vcf.gz"
        if (-not (Test-NonEmptyFile $path) -or -not (Test-NonEmptyFile "$path.tbi")) { return $false }
    }
    try {
        $state = Get-Content -LiteralPath $statePath -Raw -Encoding UTF8 | ConvertFrom-Json
        return (
            ([string]$state.schema -eq "dante-genotype-txt-normalized-v1") -and
            ([string]$state.status -eq "complete") -and
            ([string]$state.source_sha256 -eq (Get-Sha256 $source.FullName)) -and
            ([string]$state.adapter_sha256 -eq (Get-Sha256 $adapter)) -and
            ([string]$state.fasta_sha256 -eq (Get-Sha256 $fasta)) -and
            ([string]$state.fai_sha256 -eq (Get-Sha256 $fai))
        )
    } catch {
        return $false
    }
}

function Test-NormalizedLayersReady {
    if (Test-IsGenotypeTxtProfile) { return (Test-GenotypeTxtNormalizedLayersReady) }
    return (& $script:DanteGenotypeTxtOriginalTestNormalizedLayersReady)
}

function Get-NormalizedLayersMenuStatus {
    if (-not (Test-IsGenotypeTxtProfile)) {
        return (& $script:DanteGenotypeTxtOriginalGetNormalizedLayersMenuStatus)
    }
    $paths = @(
        (Join-Path $Work "Dante_SNP_normalized.vcf.gz"),
        (Join-Path $Work "Dante_INDEL_normalized.vcf.gz"),
        (Join-Path $Work "Dante_RESIDUAL_normalized.vcf.gz"),
        (Get-GenotypeTxtNormalizedStatePath)
    )
    return (Get-ArtifactMenuStatus -Ready (Test-GenotypeTxtNormalizedLayersReady) -Paths $paths)
}

function Invoke-GenotypeTxtBcftoolsNormalize([string]$InputName, [string]$OutputName, [string]$Label, [string]$BcftoolsRuntimeImage) {
    $input = Join-Path $Work $InputName
    if (-not (Test-NonEmptyFile $input)) { throw "Prepared genotype-table $Label layer is missing: $input" }
    $output = Join-Path $Work $OutputName
    $partName = ($OutputName -replace '\.vcf\.gz$', '.part.vcf.gz')
    $part = Join-Path $Work $partName
    foreach ($stale in @($part, "$part.tbi")) {
        if (Test-Path -LiteralPath $stale) { Remove-Item -LiteralPath $stale -Force }
    }
    Write-Step "Normalize genotype-table $Label layer"
    Invoke-ObservedNative -Activity "bcftools normalize genotype-table $Label" -Program "docker" -Arguments @(
        "run", "--rm",
        "--mount", "type=bind,source=$PatientRoot,target=/patient",
        "--mount", "type=bind,source=$VepData,target=/data",
        "--entrypoint", "bcftools", $BcftoolsRuntimeImage,
        "norm", "--check-ref", "e",
        "--fasta-ref", "/data/fasta/$FastaName",
        "--multiallelics", "-any",
        "--threads", "$script:BcftoolsThreads",
        "--output-type", "z",
        "--output", "/patient/Dante_work/$partName",
        "/patient/Dante_work/$InputName"
    ) -WatchFile $part
    Invoke-ObservedNative -Activity "Index genotype-table $Label normalized layer" -Program "docker" -Arguments @(
        "run", "--rm",
        "--mount", "type=bind,source=$PatientRoot,target=/patient",
        "--entrypoint", "bcftools", $BcftoolsRuntimeImage,
        "index", "--tbi", "--force", "/patient/Dante_work/$partName"
    ) -WatchFile "$part.tbi"
    Move-Item -LiteralPath $part -Destination $output -Force
    Move-Item -LiteralPath "$part.tbi" -Destination "$output.tbi" -Force
}

function Prepare-NormalizedLayers {
    if (-not (Test-IsGenotypeTxtProfile)) {
        & $script:DanteGenotypeTxtOriginalPrepareNormalizedLayers
        return
    }
    if (-not (Test-FastaReady)) { throw "GRCh37 FASTA is not installed. Select item 6 first." }
    Assert-Docker
    $source = Get-GenotypeTxtSourceFile
    $adapter = Get-GenotypeTxtAdapterPath
    $python = Get-Python
    New-Item -ItemType Directory -Force -Path $Work, $Qc | Out-Null
    if ((Test-GenotypeTxtNormalizedLayersReady) -and -not $ForceDownload) {
        Write-Host "Genotype-table reference-resolved layers are current; reusing them." -ForegroundColor Green
        return
    }

    $fasta = Get-GenotypeTxtFastaPath
    $fai = Ensure-GenotypeTxtFastaIndex
    $preparedSnp = Join-Path $Work "Dante_TXT_SNP_carried.vcf.gz"
    $preparedIndel = Join-Path $Work "Dante_TXT_INDEL_carried.vcf.gz"
    $preparedResidual = Join-Path $Work "Dante_TXT_RESIDUAL_carried.vcf.gz"
    $prepReport = Join-Path $Qc "Dante_GENOTYPE_TXT_preparation_report.txt"
    Write-Step "Resolve genotype-table calls against GRCh37"
    Write-Host "Only sequence-resolved A/C/G/T ALT-carriers enter VEP; -- and symbolic I/D calls remain outside sequence annotation." -ForegroundColor Cyan
    Invoke-ObservedNative -Activity "Prepare genotype-table reference-resolved layers" -Program $python -Arguments @(
        $adapter, "prepare-layers", $source.FullName,
        "--fasta", $fasta,
        "--fai", $fai,
        "--snp-out", $preparedSnp,
        "--indel-out", $preparedIndel,
        "--residual-out", $preparedResidual,
        "--report", $prepReport
    ) -WatchFile "$preparedSnp.part"

    $bcftoolsRuntimeImage = Get-BcftoolsRuntimeImage
    Assert-BcftoolsCompatibility $bcftoolsRuntimeImage
    Invoke-GenotypeTxtBcftoolsNormalize "Dante_TXT_SNP_carried.vcf.gz" "Dante_SNP_normalized.vcf.gz" "SNP" $bcftoolsRuntimeImage
    Invoke-GenotypeTxtBcftoolsNormalize "Dante_TXT_INDEL_carried.vcf.gz" "Dante_INDEL_normalized.vcf.gz" "INDEL" $bcftoolsRuntimeImage
    Invoke-GenotypeTxtBcftoolsNormalize "Dante_TXT_RESIDUAL_carried.vcf.gz" "Dante_RESIDUAL_normalized.vcf.gz" "RESIDUAL" $bcftoolsRuntimeImage

    $state = [ordered]@{
        schema = "dante-genotype-txt-normalized-v1"
        status = "complete"
        profile = "GENOTYPE_TABLE_TXT"
        source = $source.Name
        source_sha256 = Get-Sha256 $source.FullName
        adapter_sha256 = Get-Sha256 $adapter
        fasta_sha256 = Get-Sha256 $fasta
        fai_sha256 = Get-Sha256 $fai
        reference_resolution = "Exact GRCh37 FASTA base; bcftools --check-ref e fails closed on disagreement"
        symbolic_indels = "NOT SEQUENCE-RESOLVED; excluded from VEP"
        whole_genome_negative_evidence = $false
    }
    Set-AtomicText (Get-GenotypeTxtNormalizedStatePath) ($state | ConvertTo-Json -Depth 6)
    if (-not (Test-GenotypeTxtNormalizedLayersReady)) {
        throw "Genotype-table normalized layers were created but failed their READY validator."
    }
}

function Invoke-GenotypeTxtFullVepLayer([string]$Layer, [string]$InputName) {
    $input = Join-Path $Work $InputName
    if (-not (Test-NonEmptyFile $input)) { throw "Missing genotype-table layer: $input" }
    $upper = $Layer.ToUpper()
    $outputName = "Dante_$($upper)_full.vep.vcf.gz"
    $output = Join-Path $FullVepResults $outputName
    $state = "$output.state.sha256"
    $cacheInfo = Join-Path $VepData "homo_sapiens\116_GRCh37\info.txt"
    $fingerprint = "VEP-116-genotype-txt-GRCh37|$(Get-Sha256 $cacheInfo)|$(Get-Sha256 $input)"
    if ((Test-NonEmptyFile $output) -and (Test-NonEmptyFile $state) -and
        ((Get-Content -LiteralPath $state -Raw).Trim() -eq $fingerprint) -and -not $ForceDownload) {
        Write-Host "$upper genotype-table VEP output is current; reusing it."
    } else {
        $partName = "Dante_$($upper)_full.part.vep.vcf.gz"
        $part = Join-Path $FullVepResults $partName
        if (Test-Path -LiteralPath $part) { Remove-Item -LiteralPath $part -Force }
        Write-Step "Full VEP annotation: genotype-table $upper"
        Invoke-ObservedNative -Activity "Full VEP annotation: genotype-table $upper" -Program "docker" -Arguments @(
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
        if (-not (Test-GzipFile $part)) { throw "$upper genotype-table VEP output failed gzip validation." }
        Move-Item -LiteralPath $part -Destination $output -Force
        Set-AtomicText $state $fingerprint
    }

    $python = Get-Python
    $summary = Join-Path $FullVepResults "Dante_$($upper)_VEP_summary.tsv.gz"
    $rsidVcf = Join-Path $FullVepResults "Dante_$($upper)_rsid.vcf.gz"
    Invoke-ObservedNative -Activity "Postprocess genotype-table VEP layer: $upper" -Program $python -Arguments @(
        (Join-Path $CodeRoot "vep_postprocess.py"), $output,
        "--rsid-vcf", $rsidVcf, "--summary", $summary,
        "--report", (Join-Path $Qc "Dante_$($upper)_VEP_report.txt"),
        "--layer", $Layer
    ) -WatchFile "$summary.part"
}

function Run-FullVepAnnotation {
    if (-not (Test-IsGenotypeTxtProfile)) {
        & $script:DanteGenotypeTxtOriginalRunFullVepAnnotation
        return
    }
    if (-not (Test-VepCacheReady)) { throw "VEP 116 GRCh37 cache is missing. Select item 1 first." }
    if (-not (Test-FastaReady)) { throw "GRCh37 FASTA is missing. Select item 6 first." }
    if (-not (Test-GenotypeTxtNormalizedLayersReady)) { throw "Genotype-table small-variant layers are not READY. Select item 8 first." }
    Assert-Docker
    New-Item -ItemType Directory -Force -Path $Results, $FullVepResults, $Work, $Qc | Out-Null
    Invoke-GenotypeTxtFullVepLayer "snp" "Dante_SNP_normalized.vcf.gz"
    Invoke-GenotypeTxtFullVepLayer "indel" "Dante_INDEL_normalized.vcf.gz"
    Invoke-GenotypeTxtFullVepLayer "residual" "Dante_RESIDUAL_normalized.vcf.gz"
    Write-Host "Genotype-table VEP layers are ready in $FullVepResults" -ForegroundColor Green
}

function Run-Stage1 {
    if (Test-IsGenotypeTxtProfile) {
        Run-RsidExport
        return
    }
    & $script:DanteGenotypeTxtOriginalRunStage1
}

function Ensure-GenotypeTxtStructuralSentinels {
    if (-not (Test-IsGenotypeTxtProfile)) { return }
    New-Item -ItemType Directory -Force -Path $Results, $Qc | Out-Null
    $header = "AnnotSV_ID`tSV_chrom`tSV_start`tSV_end`tSV_type"
    foreach ($kind in @("CNV", "SV")) {
        $output = Join-Path $Results "Dante_$($kind)_AnnotSV.tsv"
        $report = Join-Path $Qc "Dante_$($kind)_AnnotSV_input_report.txt"
        Set-AtomicText $output $header
        Set-AtomicText $report @"
DANTE GENOTYPE-TABLE STRUCTURAL INPUT SENTINEL
Selected rows: 0
Scope: NOT APPLICABLE. The active genotype-table TXT source does not provide validated $kind calls.
This file exists only so Stage 24 can represent zero structural inputs without claiming they were tested.
"@
    }
}

function Get-AnnotSVOutputMenuStatus {
    if (Test-IsGenotypeTxtProfile) { return "READY" }
    return (& $script:DanteGenotypeTxtOriginalGetAnnotSVOutputMenuStatus)
}

function Run-KnowledgeEnrichment {
    if (-not (Test-IsGenotypeTxtProfile)) {
        & $script:DanteGenotypeTxtOriginalRunKnowledgeEnrichment
        return
    }
    if (-not (Test-UnifiedCatalogReady)) {
        throw "Genotype-table Stage 24 requires the current small-variant catalog. Select item 10/12 first."
    }
    Ensure-GenotypeTxtStructuralSentinels
    $oldKind = [Environment]::GetEnvironmentVariable("DANTE_INPUT_PROFILE_KIND", "Process")
    try {
        [Environment]::SetEnvironmentVariable("DANTE_INPUT_PROFILE_KIND", "GENOTYPE_TABLE_TXT", "Process")
        & $script:DanteGenotypeTxtOriginalRunKnowledgeEnrichment
    } finally {
        [Environment]::SetEnvironmentVariable("DANTE_INPUT_PROFILE_KIND", $oldKind, "Process")
    }
    $python = Get-Python
    Invoke-Native $python @(
        (Get-GenotypeTxtProfileHelperPath),
        "annotate-stage24", "--patient-root", $PatientRoot
    )
    if (-not (Test-KnowledgeEnrichmentReady)) {
        throw "Genotype-table Stage 24 completed but input-scope annotation failed the READY validator."
    }
}

function Show-AutomaticCheck {
    if (-not (Test-IsGenotypeTxtProfile)) {
        & $script:DanteGenotypeTxtOriginalShowAutomaticCheck
        return
    }
    $profile = Get-DanteInputProfile
    $source = Get-GenotypeTxtSourceFile
    Write-Step "Automatic startup check"
    Write-CheckStatus "Code directory" (Test-Path -LiteralPath $CodeRoot -PathType Container) $CodeRoot
    Write-CheckStatus "Active patient" $true "$Patient | $PatientRoot"
    Write-CheckStatus "Input profile" $true "$($profile.kind) | $($profile.source) | $($profile.assembly)"
    Write-CheckStatus "Genotype TXT" $true "$($source.Name) | $([Math]::Round($source.Length / 1MB, 2)) MiB"
    Write-Host "Scope: $($profile.scope)" -ForegroundColor Yellow
    Write-Host "Sequence annotation: A/C/G/T calls only after exact GRCh37 REF resolution. Symbolic I/D indels remain unannotated." -ForegroundColor Yellow
    Write-Host "Not available from this source: WGS coverage, FASTQ evidence, CNV/SV calls, telomere estimate, BAM browser." -ForegroundColor DarkGray
    $dockerCli = [bool](Get-Command docker -ErrorAction SilentlyContinue)
    Write-CheckStatus "Docker CLI" $dockerCli $(if ($dockerCli) { "available" } else { "not found" })
    $pythonCommand = Get-Command py -ErrorAction SilentlyContinue
    if (-not $pythonCommand) { $pythonCommand = Get-Command python -ErrorAction SilentlyContinue }
    Write-CheckStatus "Python 3" ([bool]$pythonCommand) $(if ($pythonCommand) { $pythonCommand.Source } else { "not found" })
}

function Select-Mode {
    if (-not (Test-IsGenotypeTxtProfile)) {
        return (& $script:DanteGenotypeTxtOriginalSelectMode)
    }

    Show-AutomaticCheck
    $profile = Get-DanteInputProfile
    Write-Host ""
    Write-Host ("Dante genomic pipeline v{0}" -f (Get-DanteDisplayVersion)) -ForegroundColor Cyan
    Write-Host ("Patient: {0}" -f $Patient) -ForegroundColor Green
    Write-Host ("Workspace: {0}" -f $PatientRoot) -ForegroundColor DarkGray
    Write-Host ("Input: {0} | {1} | {2}" -f $profile.kind, $profile.source, $profile.assembly) -ForegroundColor Cyan
    Write-Host "Scope: interrogated genotype-table loci only; no-calls and absence outside the table are UNKNOWN." -ForegroundColor Yellow
    Write-Host "Symbolic I/D indels are retained in export/QC but cannot enter VEP without sequence alleles." -ForegroundColor Yellow
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

    Write-Host "GENOTYPE-TABLE EXPORT" -ForegroundColor DarkCyan
    Write-MenuItem "2" "Create measured-loci 23andMe-style genotype export" $rsidStatus
    Write-Host ""
    Write-Host "SMALL-VARIANT ANNOTATION - REFERENCE-RESOLVED SNP CARRIERS ONLY" -ForegroundColor DarkCyan
    Write-MenuItem "1" "Install/resume VEP 116 GRCh37 cache" $vepCacheStatus
    Write-MenuItem "6" "Install GRCh37 reference FASTA" $fastaStatus
    Write-MenuItem "7" "Install/update current ClinVar GRCh37" $clinVarStatus
    Write-MenuItem "8" "Resolve TXT genotypes against GRCh37 and normalize SNP carriers" $normalizedStatus
    Write-MenuItem "9" "Run full VEP annotation of sequence-resolved measured SNPs" $fullVepStatus
    Write-MenuItem "10" "Build unified annotated measured-variant catalog" $catalogStatus
    Write-MenuItem "11" "Build deterministic shortlist within sequence-resolved measured loci" $shortlistStatus
    Write-MenuItem "12" "Run complete genotype-table small-variant annotation" $stage2Status
    Write-Host ""
    Write-Host "REFERENCE KNOWLEDGE + REPORTS" -ForegroundColor DarkCyan
    Write-MenuItem "22" "Install/update open reference knowledge bundle" $openReferenceStatus
    Write-MenuItem "23" "Install research only knowledge bundle (NO COMMERCIAL USE ALLOWED)" $researchReferenceStatus
    Write-MenuItem "24" "Enrich measured genotype-table variants with reference knowledge (NO OPENAI)" $knowledgeEnrichmentStatus
    Write-MenuItem "25" "Generate CLIENT PDF report (TXT SCOPE / OPEN / NO AI)" $commercialReportStatus
    Write-MenuItem "26" "Generate RESEARCH ONLY PDF report (TXT SCOPE / NON-COMMERCIAL)" $researchReportStatus
    Write-Host ""
    Write-Host "AUTOMATION" -ForegroundColor DarkCyan
    Write-MenuItem "21" "Run all applicable genotype-table processing (local, no API)" $allLocalStatus
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
        default { throw "This action is not available for the active genotype-table TXT profile: $choice" }
    }
}
