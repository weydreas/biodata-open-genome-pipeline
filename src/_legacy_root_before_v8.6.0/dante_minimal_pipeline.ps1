param(
    [ValidateSet("menu", "check", "setup-rsid", "setup-annotsv", "rsid", "annotsv", "all")]
    [string]$Mode = "menu",
    [int]$Threads = 0,
    [switch]$ForceDownload,
    [switch]$KeepIntermediate,
    [switch]$NoPause
)

Set-StrictMode -Version 2.0
$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$AnnotationData = Join-Path $Root "annotation_data"
$VepData = Join-Path $AnnotationData "vep_rsid"
$AnnotSVData = Join-Path $AnnotationData "annotsv"
$Results = Join-Path $Root "Dante_results"
$Work = Join-Path $Root "Dante_work"
$Qc = Join-Path $Root "Dante_qc"

$VepImage = "ensemblorg/ensembl-vep:release_116.1"
$VepCacheVersion = "116"
$VepCacheArchiveName = "homo_sapiens_vep_116_GRCh37.tar.gz"
$VepCacheUrl = "https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/$VepCacheArchiveName"
$AnnotSVImage = "quay.io/biocontainers/annotsv:3.5.10--hdfd78af_0"

if ($Threads -le 0) {
    $Threads = [Math]::Max(2, [Math]::Min(4, [Environment]::ProcessorCount - 1))
}

function Write-Step([string]$Text) {
    Write-Host ""
    Write-Host ("=" * 72) -ForegroundColor DarkCyan
    Write-Host $Text -ForegroundColor Cyan
    Write-Host ("=" * 72) -ForegroundColor DarkCyan
}

function Invoke-Native([string]$Program, [string[]]$Arguments) {
    # Native programs write normal progress to stderr.  With the script-wide
    # ErrorActionPreference=Stop, PowerShell 7 can otherwise turn that stream
    # into a terminating PowerShell error before LASTEXITCODE is inspected.
    $savedErrorActionPreference = $ErrorActionPreference
    $nativeExitCode = -1
    try {
        $ErrorActionPreference = "Continue"
        & $Program @Arguments 2>&1 | ForEach-Object { Write-Host $_ }
        $nativeExitCode = $LASTEXITCODE
    } finally {
        $ErrorActionPreference = $savedErrorActionPreference
    }
    if ($nativeExitCode -ne 0) {
        throw "Command failed with exit code $nativeExitCode`: $Program"
    }
}

function Invoke-NativeWithRetry(
    [string]$Program,
    [string[]]$Arguments,
    [int]$MaxAttempts = 5,
    [int]$InitialDelaySeconds = 10
) {
    $delaySeconds = $InitialDelaySeconds
    for ($attempt = 1; $attempt -le $MaxAttempts; $attempt++) {
        try {
            Invoke-Native $Program $Arguments
            return
        } catch {
            if ($attempt -ge $MaxAttempts) {
                throw
            }
            Write-Host ""
            Write-Host "Network operation failed (attempt $attempt of $MaxAttempts)." -ForegroundColor Yellow
            Write-Host "Completed data will be reused when possible; retrying in $delaySeconds seconds..." -ForegroundColor Yellow
            Start-Sleep -Seconds $delaySeconds
            $delaySeconds = [Math]::Min(60, $delaySeconds * 2)
        }
    }
}

function Get-Python {
    if (Get-Command py -ErrorAction SilentlyContinue) { return "py" }
    if (Get-Command python -ErrorAction SilentlyContinue) { return "python" }
    throw "Python 3 was not found. Install it and make py.exe or python.exe available in PATH."
}

function Get-Curl {
    $curl = Get-Command curl.exe -ErrorAction SilentlyContinue
    if ($curl) { return $curl.Source }
    throw "curl.exe was not found. It is included with supported Windows 10/11 installations."
}

function Format-ByteSize([long]$Bytes) {
    if ($Bytes -ge 1GB) { return ("{0:N2} GiB" -f ($Bytes / 1GB)) }
    if ($Bytes -ge 1MB) { return ("{0:N1} MiB" -f ($Bytes / 1MB)) }
    if ($Bytes -ge 1KB) { return ("{0:N1} KiB" -f ($Bytes / 1KB)) }
    return "$Bytes bytes"
}

function Invoke-ResumableDownload(
    [string]$Url,
    [string]$Destination,
    [int]$MaxConsecutiveNoProgress = 20
) {
    # curl's generic --retry policy does not reliably retry error 56
    # (connection reset).  Restart curl explicitly so every attempt rereads
    # the current file length and requests only the missing byte range.
    $curl = Get-Curl
    $attempt = 0
    $consecutiveNoProgress = 0

    while ($true) {
        $attempt++
        $before = 0L
        if (Test-Path -LiteralPath $Destination -PathType Leaf) {
            $before = (Get-Item -LiteralPath $Destination).Length
        }

        Write-Host "Download attempt $attempt; already saved: $(Format-ByteSize $before)"
        try {
            Invoke-Native $curl @(
                "--fail", "--location", "--http1.1",
                "--connect-timeout", "60",
                "--speed-limit", "1024", "--speed-time", "300",
                "--continue-at", "-",
                "--output", $Destination,
                $Url
            )
            $after = (Get-Item -LiteralPath $Destination).Length
            Write-Host "Download command completed; saved: $(Format-ByteSize $after)"
            return
        } catch {
            $after = 0L
            if (Test-Path -LiteralPath $Destination -PathType Leaf) {
                $after = (Get-Item -LiteralPath $Destination).Length
            }
            $gained = [Math]::Max(0L, $after - $before)

            if ($gained -gt 0) {
                $consecutiveNoProgress = 0
            } else {
                $consecutiveNoProgress++
            }

            Write-Host ""
            Write-Host "Connection interrupted. The partial file is preserved." -ForegroundColor Yellow
            Write-Host "Saved: $(Format-ByteSize $after); added this attempt: $(Format-ByteSize $gained)." -ForegroundColor Yellow

            if ($consecutiveNoProgress -ge $MaxConsecutiveNoProgress) {
                throw "The download made no progress for $MaxConsecutiveNoProgress consecutive attempts. The partial file is intact; check the network and run item 2 again."
            }

            $delaySeconds = [Math]::Min(120, 10 + (10 * $consecutiveNoProgress))
            Write-Host "Automatically resuming in $delaySeconds seconds. Press Ctrl+C only if you want to pause." -ForegroundColor Yellow
            Start-Sleep -Seconds $delaySeconds
        }
    }
}

function Assert-Docker {
    if (-not (Get-Command docker -ErrorAction SilentlyContinue)) {
        throw "Docker was not found. Install and start Docker Desktop for Windows."
    }
    & docker info *> $null
    if ($LASTEXITCODE -ne 0) {
        throw "Docker Desktop is installed, but its engine is not running."
    }
}

function Find-OneDanteFile([string]$Role, [string]$Pattern) {
    $found = @(Get-ChildItem -LiteralPath $Root -File | Where-Object { $_.Name -match $Pattern })
    if ($found.Count -eq 0) {
        throw "Missing $Role VCF beside the scripts (expected pattern: $Pattern)."
    }
    if ($found.Count -gt 1) {
        throw "More than one $Role VCF found: $($found.Name -join ', ')"
    }
    return $found[0]
}

function Get-DanteInputs {
    return @{
        snp = Find-OneDanteFile "filtered SNP" '\.filtered\.snp\.vcf(\.gz)?$'
        cnv = Find-OneDanteFile "CNV" '\.cnv\.vcf(\.gz)?$'
        sv  = Find-OneDanteFile "SV" '\.sv\.vcf(\.gz)?$'
    }
}

function Test-VepCacheReady([string]$CacheDir) {
    $info = Join-Path $CacheDir "info.txt"
    $chr1 = Join-Path $CacheDir "1"
    return (
        (Test-Path -LiteralPath $CacheDir -PathType Container) -and
        (Test-Path -LiteralPath $info -PathType Leaf) -and
        (Test-Path -LiteralPath $chr1 -PathType Container)
    )
}

function Test-VepCacheArchive([string]$ArchiveName) {
    $archivePath = Join-Path $VepData $ArchiveName
    if (-not (Test-Path -LiteralPath $archivePath -PathType Leaf)) {
        return $false
    }
    $savedErrorActionPreference = $ErrorActionPreference
    try {
        $ErrorActionPreference = "Continue"
        & docker @(
            "run", "--rm",
            "--mount", "type=bind,source=$VepData,target=/data",
            $VepImage,
            "bash", "-lc", "gzip -t /data/$ArchiveName"
        ) *> $null
        return ($LASTEXITCODE -eq 0)
    } finally {
        $ErrorActionPreference = $savedErrorActionPreference
    }
}

function Show-Check {
    Write-Step "Environment and input check"
    $inputs = Get-DanteInputs
    Assert-Docker
    $python = Get-Python
    Write-Host "Root: $Root"
    Write-Host "Python: $python"
    Write-Host "Docker: available and running"
    Write-Host "VEP workers: $Threads"
    foreach ($role in @("snp", "cnv", "sv")) {
        Write-Host (("{0,7}: {1}" -f $role, $inputs[$role].Name))
    }
    Write-Host "Only these three Dante VCF files are required by the minimal pipeline."
}

function Install-RsidData {
    Write-Step "Install Ensembl VEP 116 GRCh37 cache for rsID matching"
    Assert-Docker
    New-Item -ItemType Directory -Force -Path $VepData | Out-Null
    $cacheDir = Join-Path $VepData "homo_sapiens\116_GRCh37"
    $cacheArchive = Join-Path $VepData $VepCacheArchiveName
    if ((Test-VepCacheReady $cacheDir) -and -not $ForceDownload) {
        Write-Host "VEP cache already exists: $cacheDir"
        return
    }
    Invoke-NativeWithRetry docker @("pull", $VepImage)

    if ($ForceDownload) {
        if (Test-Path -LiteralPath $cacheArchive) {
            Remove-Item -LiteralPath $cacheArchive -Force
        }
        if (Test-Path -LiteralPath $cacheDir) {
            Remove-Item -LiteralPath $cacheDir -Recurse -Force
        }
    }

    # VEP 115/116 Docker images currently lack HTML::TableExtract, so
    # INSTALL.pl --USE_HTTPS_PROTO silently skips the cache. Download the
    # official indexed cache directly, which Ensembl documents as supported.
    if (-not (Test-VepCacheArchive $VepCacheArchiveName)) {
        Write-Host "Downloading official Ensembl indexed cache. Network resets will be resumed automatically."
        Invoke-ResumableDownload $VepCacheUrl $cacheArchive
    } else {
        Write-Host "A complete downloaded cache archive was found; reusing it."
    }

    if (-not (Test-VepCacheArchive $VepCacheArchiveName)) {
        if (Test-Path -LiteralPath $cacheArchive) {
            Remove-Item -LiteralPath $cacheArchive -Force
        }
        throw "The downloaded VEP cache archive failed its gzip integrity check. Run item 2 again to redownload it."
    }

    if (Test-Path -LiteralPath $cacheDir) {
        Remove-Item -LiteralPath $cacheDir -Recurse -Force
    }
    Write-Host "Archive integrity is valid. Extracting the GRCh37 cache..."
    Invoke-Native docker @(
        "run", "--rm",
        "--mount", "type=bind,source=$VepData,target=/data",
        $VepImage,
        "bash", "-lc", "tar -xzf /data/$VepCacheArchiveName -C /data"
    )
    if (-not (Test-VepCacheReady $cacheDir)) {
        throw "The official archive was extracted, but the expected GRCh37 cache structure was not found: $cacheDir"
    }
    if (Test-Path -LiteralPath $cacheArchive) {
        Remove-Item -LiteralPath $cacheArchive -Force
    }
    Write-Host "The rsID cache is ready. A reference FASTA was intentionally not downloaded."
}

function Install-AnnotSVData {
    Write-Step "Install AnnotSV 3.5.10 public annotations"
    Assert-Docker
    New-Item -ItemType Directory -Force -Path $AnnotSVData | Out-Null
    $bundle = Join-Path $AnnotSVData "AnnotSV_annotations"
    if ((Test-Path -LiteralPath $bundle) -and -not $ForceDownload) {
        Write-Host "AnnotSV annotations already exist: $bundle"
        return
    }
    Invoke-NativeWithRetry docker @("pull", $AnnotSVImage)
    Invoke-Native docker @(
        "run", "--rm",
        "--mount", "type=bind,source=$AnnotSVData,target=/data",
        "--workdir", "/data",
        $AnnotSVImage,
        "INSTALL_annotations.sh"
    )
    if (-not (Test-Path -LiteralPath $bundle)) {
        throw "AnnotSV installer finished, but AnnotSV_annotations was not found in $AnnotSVData"
    }
    Write-Host "AnnotSV public annotations are ready."
}

function Run-RsidExport {
    Show-Check
    $cacheDir = Join-Path $VepData "homo_sapiens\116_GRCh37"
    if (-not (Test-VepCacheReady $cacheDir)) {
        throw "VEP GRCh37 cache is not installed. Select menu item 2 first."
    }
    New-Item -ItemType Directory -Force -Path $Results, $Work, $Qc | Out-Null
    $inputs = Get-DanteInputs
    $python = Get-Python
    $vepOutputName = "Dante_SNP_rsid_lookup.vep.vcf.gz"
    $vepOutput = Join-Path $Work $vepOutputName
    $warningName = "Dante_SNP_rsid_lookup.warnings.txt"

    Write-Step "Allele-specific rsID matching for Dante SNPs"
    Invoke-Native docker @(
        "run", "--rm",
        "--mount", "type=bind,source=$Root,target=/work",
        "--mount", "type=bind,source=$VepData,target=/data",
        $VepImage,
        "vep",
        "--cache", "--offline",
        "--cache_version", $VepCacheVersion,
        "--species", "homo_sapiens",
        "--assembly", "GRCh37",
        "--dir_cache", "/data",
        "--format", "vcf",
        "--vcf",
        "--check_existing",
        "--exclude_null_alleles",
        "--allele_number",
        "--fields", "Allele,Consequence,Existing_variation,ALLELE_NUM",
        "--fork", "$Threads",
        "--buffer_size", "10000",
        "--no_stats",
        "--compress_output", "bgzip",
        "--force_overwrite",
        "--input_file", "/work/$($inputs.snp.Name)",
        "--output_file", "/work/Dante_work/$vepOutputName",
        "--warning_file", "/work/Dante_qc/$warningName"
    )

    Write-Step "Create strict four-column 23andMe-style export"
    $export = Join-Path $Results "Dante_23andMe_style.txt.gz"
    $report = Join-Path $Qc "Dante_23andMe_style_report.txt"
    Invoke-Native $python @(
        (Join-Path $Root "vep_to_23andme.py"),
        $vepOutput,
        "--output", $export,
        "--report", $report
    )
    if (-not (Test-Path -LiteralPath $export)) {
        throw "The 23andMe-style output was not created."
    }
    if (-not $KeepIntermediate) {
        Remove-Item -LiteralPath $vepOutput -Force
    }
    Write-Host "Created: $export"
}

function Invoke-AnnotSVOne(
    [string]$Python,
    [System.IO.FileInfo]$InputFile,
    [string]$Kind
) {
    $upper = $Kind.ToUpper()
    $preparedName = "Dante_$($upper)_PASS_canonical.vcf"
    $prepared = Join-Path $Work $preparedName
    $prepareReport = Join-Path $Qc "Dante_$($upper)_AnnotSV_input_report.txt"
    Invoke-Native $Python @(
        (Join-Path $Root "prepare_annotsv_input.py"),
        $InputFile.FullName,
        $prepared,
        "--kind", $Kind,
        "--report", $prepareReport
    )

    $annotatedName = "Dante_$($upper)_AnnotSV.tsv"
    Write-Step "AnnotSV annotation: $upper"
    Invoke-Native docker @(
        "run", "--rm",
        "--mount", "type=bind,source=$Root,target=/work",
        $AnnotSVImage,
        "AnnotSV",
        "-SVinputFile", "/work/Dante_work/$preparedName",
        "-annotationsDir", "/work/annotation_data/annotsv/AnnotSV_annotations",
        "-genomeBuild", "GRCh37",
        "-outputDir", "/work/Dante_results",
        "-outputFile", $annotatedName,
        "-overwrite", "1",
        "-includeCI", "1",
        "-SVinputInfo", "1",
        "-annotationMode", "both",
        "-tx", "RefSeq",
        "-vcf", "0"
    )
    $annotated = Join-Path $Results $annotatedName
    if (-not (Test-Path -LiteralPath $annotated)) {
        throw "AnnotSV did not create $annotatedName"
    }
    if (-not $KeepIntermediate) {
        Remove-Item -LiteralPath $prepared -Force
    }
    Write-Host "Created: $annotated"
}

function Run-StructuralAnnotation {
    Show-Check
    $bundle = Join-Path $AnnotSVData "AnnotSV_annotations"
    if (-not (Test-Path -LiteralPath $bundle)) {
        throw "AnnotSV data are not installed. Select menu item 4 first."
    }
    New-Item -ItemType Directory -Force -Path $Results, $Work, $Qc | Out-Null
    $inputs = Get-DanteInputs
    $python = Get-Python
    Invoke-AnnotSVOne $python $inputs.cnv "cnv"
    Invoke-AnnotSVOne $python $inputs.sv "sv"
}

function Select-Mode {
    Write-Host "Dante WGS minimal local export" -ForegroundColor Cyan
    Write-Host "1. Check Docker, Python and the three required Dante VCF files"
    Write-Host "2. Download/install the GRCh37 rsID cache"
    Write-Host "3. Create Dante_23andMe_style.txt.gz"
    Write-Host "4. Download/install AnnotSV public annotations"
    Write-Host "5. Create Dante_CNV_AnnotSV.tsv and Dante_SV_AnnotSV.tsv"
    Write-Host "6. Install everything and create all three outputs"
    Write-Host "0. Exit"
    $choice = Read-Host "Select"
    switch ($choice) {
        "1" { return "check" }
        "2" { return "setup-rsid" }
        "3" { return "rsid" }
        "4" { return "setup-annotsv" }
        "5" { return "annotsv" }
        "6" { return "all" }
        "0" { return "exit" }
        default { throw "Unknown menu choice: $choice" }
    }
}

$exitCode = 0
try {
    if ($Mode -eq "menu") { $Mode = Select-Mode }
    switch ($Mode) {
        "exit" { }
        "check" { Show-Check }
        "setup-rsid" { Install-RsidData }
        "setup-annotsv" { Install-AnnotSVData }
        "rsid" { Run-RsidExport }
        "annotsv" { Run-StructuralAnnotation }
        "all" {
            Install-RsidData
            Install-AnnotSVData
            Run-RsidExport
            Run-StructuralAnnotation
        }
    }
    if ($Mode -ne "exit") {
        Write-Step "Done"
        Write-Host "Final outputs: $Results"
        Write-Host "Quality-control reports: $Qc"
        Write-Host "Original Dante VCF files were not modified."
        Write-Host "No OpenAI API or literature-search component was used."
    }
} catch {
    $exitCode = 1
    Write-Host ""
    Write-Host "ERROR: $($_.Exception.Message)" -ForegroundColor Red
    Write-Host $_.ScriptStackTrace -ForegroundColor DarkRed
}

if (-not $NoPause) {
    [void](Read-Host "Press Enter to close")
}
exit $exitCode
