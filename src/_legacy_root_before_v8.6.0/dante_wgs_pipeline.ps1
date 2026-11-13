param(
    [ValidateSet(
        "menu", "setup-vep-cache", "rsid", "setup-annotsv", "annotsv", "stage1",
        "setup-fasta", "setup-clinvar", "prepare-small", "vep-full", "catalog",
        "shortlist", "stage2", "reference-targets", "reference-source-merge", "reference-align",
        "reference-call", "reference-export", "stage3", "telomere", "ai", "local-all", "open-reference", "research-reference", "genome-browser", "exit"
    )]
    [string]$Mode = "menu",
    [string]$Patient = "",
    [int]$Threads = 0,
    [ValidateSet("auto", "max", "balanced")]
    [string]$PerformanceProfile = "auto",
    [switch]$ForceDownload,
    [switch]$KeepIntermediate,
    [switch]$NoPause
)

Set-StrictMode -Version 2.0
$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$PatientsRoot = Join-Path $Root "patients"
$ActivePatientState = Join-Path $Root ".active_patient.json"

if (-not $Patient) {
    if (-not (Test-Path -LiteralPath $ActivePatientState -PathType Leaf)) {
        throw "No active patient is selected. Launch start.cmd and choose a patient first, or pass -Patient explicitly."
    }
    try {
        $activePatient = Get-Content -LiteralPath $ActivePatientState -Raw -Encoding UTF8 | ConvertFrom-Json
        if ([string]$activePatient.schema -ne "dante-active-patient-v1") {
            throw "Unsupported active-patient state schema."
        }
        $Patient = [string]$activePatient.name
    } catch {
        throw "Active patient state is invalid: $ActivePatientState"
    }
}
$Patient = $Patient.Trim()
if (-not $Patient) { throw "Patient name is empty." }
if ($Patient -match '[<>:"/\\|?*]') { throw "Patient name contains an invalid Windows path character." }

$PatientRoot = Join-Path $PatientsRoot $Patient
$PatientInput = Join-Path $PatientRoot "input"
$PatientEnvFile = Join-Path $PatientRoot ".env"
if (-not (Test-Path -LiteralPath $PatientRoot -PathType Container)) {
    throw "Patient workspace does not exist: $PatientRoot. Launch start.cmd and create/select the patient first."
}
New-Item -ItemType Directory -Force -Path $PatientInput | Out-Null

$AnnotationData = Join-Path $Root "annotation_data"
$VepData = Join-Path $AnnotationData "vep_rsid"
$VepFastaData = Join-Path $VepData "fasta"
$ClinVarData = Join-Path $AnnotationData "clinvar_grch37"
$AnnotSVData = Join-Path $AnnotationData "annotsv"
$BwaIndexData = Join-Path $AnnotationData "bwa_grch37"
$ReferencePanelSources = Join-Path $PatientInput "reference_panel_sources"
$Results = Join-Path $PatientRoot "Dante_results"
$FullVepResults = Join-Path $Results "full_vep"
$Work = Join-Path $PatientRoot "Dante_work"
$Qc = Join-Path $PatientRoot "Dante_qc"
$EnvFile = Join-Path $Root ".env"
$ActivityRunner = Join-Path $Root "dante_activity.py"
$OpenReferenceManager = Join-Path $Root "open_reference_manager.py"
$OpenReferenceData = Join-Path $Root "reference_data\open_reference"
$OpenReferenceManifest = Join-Path $OpenReferenceData "open_reference_manifest.json"
$ResearchReferenceManager = Join-Path $Root "research_reference_manager.py"
$ResearchReferenceData = Join-Path $Root "reference_data\research_reference"
$ResearchReferenceManifest = Join-Path $ResearchReferenceData "research_reference_manifest.json"
$ReportPublicationGuard = Join-Path $Root "report_publication_guard.py"
$GenomeBrowserApp = Join-Path $Root "genome_browser.py"
$GenomeBrowserRequirements = Join-Path $Root "requirements-browser.txt"

$VepImage = "ensemblorg/ensembl-vep:release_116.1"
$VepCacheVersion = "116"
$VepCacheArchiveName = "homo_sapiens_vep_116_GRCh37.tar.gz"
$VepCacheUrl = "https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/$VepCacheArchiveName"
$FastaArchiveName = "Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
$FastaName = "Homo_sapiens.GRCh37.dna.primary_assembly.fa"
$FastaUrl = "https://ftp.ensembl.org/pub/grch37/release-116/fasta/homo_sapiens/dna/$FastaArchiveName"
$ClinVarUrl = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz"
$ClinVarIndexUrl = "$ClinVarUrl.tbi"
$ClinVarMd5Url = "$ClinVarUrl.md5"
$AnnotSVImage = "quay.io/biocontainers/annotsv:3.5.10--hdfd78af_0"
$AnnotSVHumanVersion = "3.5"
$AnnotSVHumanArchiveName = "Annotations_Human_$AnnotSVHumanVersion.tar.gz"
$AnnotSVHumanArchiveUrl = "https://www.lbgi.fr/~geoffroy/Annotations/$AnnotSVHumanArchiveName"
$AnnotSVCompleteMarkerName = ".dante_annotsv_3.5.10_grch37.complete"
$BcftoolsImage = "quay.io/biocontainers/bcftools:1.24--h487d631_1"
$BcftoolsFallbackImages = @("staphb/bcftools:1.24")
$BcftoolsToolContract = "bcftools-1.24"
$BwaImage = "quay.io/biocontainers/bwa:0.7.19--h577a1d6_1"
$BwaFallbackImages = @("staphb/bwa:0.7.19")
$BwaToolContract = "bwa-mem-0.7.19"
$SamtoolsImage = "quay.io/biocontainers/samtools:1.22.1--h96c455f_0"
$SamtoolsFallbackImages = @("staphb/samtools:1.22.1")
$SamtoolsToolContract = "samtools-1.22.1"
$script:ResolvedBcftoolsImage = $null
$script:ResolvedBwaImage = $null
$script:ResolvedSamtoolsImage = $null
$ReferenceContigTargetsName = "GRCh37_reference_contigs.targets.tsv"
$ReferenceContigAuditName = "Dante_reference_contig_audit.txt"
$ReferenceTargetsName = "Dante_reference_targets.tsv"
$ReferenceTargetsBedName = "Dante_reference_targets.bed"
$ReferenceTargetsReportName = "Dante_reference_targets_report.txt"
$ReferenceCoordinateConflictsName = "Dante_reference_target_coordinate_conflicts.tsv"
$ReferenceSourceGenotypeConflictsName = "Dante_reference_source_genotype_conflicts.tsv"
$ReferenceProvenanceName = "Dante_reference_genotype_provenance.tsv"
$ReferenceFastqEvidenceName = "Dante_reference_FASTQ_evidence.tsv.gz"
$ReferenceBamName = "Dante_reference_targets.bam"
$ReferenceCallsName = "Dante_reference_calls.vcf.gz"
$ReferenceCallsMetadataName = "Dante_reference_calls.metadata.json"
$ReferenceCallsUncallableName = "Dante_reference_calls_uncallable.tsv"
$ReferenceExportName = "Dante_23andMe_reference_complete.txt.gz"
$TelomereReportName = "Dante_telomere_report.txt"
$TelomereMetricsName = "Dante_telomere_metrics.json"
$TelomereStateName = "Dante_telomere_report.state.json"
$ReferenceExportMinTier = "B"
$ReferenceTierAMinDp = 10
$ReferenceTierAMinQuality = 30
$ReferenceTierAMinAlleleDepth = 3
$ReferenceTierAHomFraction = 0.90
$ReferenceTierAHetMinFraction = 0.20
$ReferenceTierAHetMaxFraction = 0.80
$ReferenceTierBMinDp = 5
$ReferenceTierBMinQuality = 10
$ReferenceTierBHomMinSupport = 4
$ReferenceTierBHomFraction = 0.80
$ReferenceTierBHetMinAlleleDepth = 2
$ReferenceTierBHetMinFraction = 0.20
$ReferenceTierBHetMaxFraction = 0.80
$ReferenceTierCMinDp = 5
$ReferenceTierCMinQuality = 5
$script:VepCacheRemoteLength = [long]0
$script:HashCache = @{}

$script:ThreadsExplicit = ($Threads -gt 0)
$script:RequestedThreads = $Threads
$script:LogicalCpuCount = [Math]::Max(1, [Environment]::ProcessorCount)
$script:HostMemoryBytes = [long]0
try {
    $computerSystem = Get-CimInstance Win32_ComputerSystem -ErrorAction Stop
    $script:HostMemoryBytes = [long]$computerSystem.TotalPhysicalMemory
} catch {
    $script:HostMemoryBytes = [long]0
}
if ($Threads -le 0) {
    $cpuRatio = if ($PerformanceProfile -eq "balanced") { 0.50 } else { 0.75 }
    $Threads = [Math]::Max(2, [Math]::Min($script:LogicalCpuCount, [Math]::Floor($script:LogicalCpuCount * $cpuRatio)))
}
$script:GeneralThreads = [int]$Threads
$script:BwaThreads = [int]$Threads
$script:SamtoolsSortThreads = [Math]::Max(1, $script:LogicalCpuCount - $script:BwaThreads)
$script:SamtoolsIndexThreads = [int]$Threads
$script:BcftoolsThreads = [Math]::Max(2, [Math]::Floor($script:LogicalCpuCount / 2))
$script:SamtoolsSortMemory = "512M"
$script:DockerCpuCount = 0
$script:DockerMemoryBytes = [long]0
$script:ResourceProfilePrinted = $false

function Write-Step([string]$Text) {
    Write-Host ""
    Write-Host ("=" * 78) -ForegroundColor DarkCyan
    Write-Host $Text -ForegroundColor Cyan
    Write-Host ("=" * 78) -ForegroundColor DarkCyan
}

function Write-CheckStatus([string]$Name, [bool]$Ok, [string]$Detail) {
    if ($Ok) {
        Write-Host "[OK]     " -NoNewline -ForegroundColor Green
    } else {
        Write-Host "[NOT OK] " -NoNewline -ForegroundColor Red
    }
    Write-Host ("{0,-22} {1}" -f $Name, $Detail)
}

function Write-MenuItem([string]$Number, [string]$Text, [string]$Status = "") {
    Write-Host ("{0,2}. {1}" -f $Number, $Text) -NoNewline
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
        Write-Host "  [$Status]" -ForegroundColor $color
    } else {
        Write-Host ""
    }
}

function Invoke-Native([string]$Program, [string[]]$Arguments) {
    # Native tools legitimately write progress to stderr. PowerShell 7 can
    # otherwise turn that progress into a terminating PowerShell exception.
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
        throw "Command failed with exit code ${nativeExitCode}: $Program"
    }
}

function Invoke-NativeCapture([string]$Program, [string[]]$Arguments) {
    # Capture combined stdout/stderr without treating a non-zero native exit code
    # as a PowerShell exception.  This is used for version/help probes where
    # tools such as `bwa` legitimately return non-zero after printing usage.
    $savedErrorActionPreference = $ErrorActionPreference
    $nativeExitCode = -1
    $captured = @()
    try {
        $ErrorActionPreference = "Continue"
        $captured = @(& $Program @Arguments 2>&1 | ForEach-Object { [string]$_ })
        $nativeExitCode = $LASTEXITCODE
    } finally {
        $ErrorActionPreference = $savedErrorActionPreference
    }
    return [pscustomobject]@{
        ExitCode = [int]$nativeExitCode
        Output = ($captured -join "`n")
    }
}


function Set-DockerResourceProfile {
    $probe = Invoke-NativeCapture "docker" @("info", "--format", "{{.NCPU}}|{{.MemTotal}}")
    $dockerCpu = $script:LogicalCpuCount
    $dockerMem = [long]0
    if (($probe.ExitCode -eq 0) -and ($probe.Output -match '^\s*(\d+)\|(\d+)\s*$')) {
        $dockerCpu = [Math]::Max(1, [int]$Matches[1])
        $dockerMem = [long]$Matches[2]
    }
    $usableCpu = [Math]::Max(1, [Math]::Min($script:LogicalCpuCount, $dockerCpu))
    $pipelineCpu = $usableCpu
    if ((-not $script:ThreadsExplicit) -and ($PerformanceProfile -ne "balanced") -and ($usableCpu -ge 6)) {
        # Keep one logical CPU outside the heavy pipeline so Windows/Docker/WSL
        # scheduling and I/O do not become the bottleneck at 100% saturation.
        $pipelineCpu = $usableCpu - 1
    }
    if ($script:ThreadsExplicit) {
        $general = [Math]::Max(1, [Math]::Min([int]$script:RequestedThreads, $usableCpu))
        $bwa = $general
    } elseif ($PerformanceProfile -eq "max") {
        # Home-computer MAX profile: use all 12 logical CPUs for BWA on the
        # 6C/12T Ryzen 5 7500F. samtools is intentionally allowed to overlap
        # with four workers; Windows/WSL schedules the short overlap periods.
        $bwa = [Math]::Max(2, [Math]::Min(12, $usableCpu))
        $general = $bwa
    } else {
        $ratio = if ($PerformanceProfile -eq "balanced") { 0.50 } else { 0.75 }
        $general = [Math]::Max(2, [Math]::Min($pipelineCpu, [Math]::Floor($pipelineCpu * $ratio)))
        $bwa = $general
    }
    $sortBaseCpu = if ($PerformanceProfile -eq "max") { $usableCpu } else { $pipelineCpu }
    $sort = [Math]::Max(1, $sortBaseCpu - $bwa)
    if (($PerformanceProfile -eq "max") -and ($usableCpu -ge 12) -and (-not $script:ThreadsExplicit)) {
        # Intentional oversubscription for the 12-thread Home Computer:
        # BWA receives 12 workers while samtools sort receives four. BWA and sort
        # are bursty rather than continuously saturating all workers together.
        $sort = 4
    } elseif (($usableCpu -ge 4) -and ($sort -lt 2)) {
        $sort = 2
    }
    $sort = [Math]::Min(8, $sort)
    $index = [Math]::Max(1, [Math]::Min($usableCpu, $general))
    $bcf = [Math]::Max(1, [Math]::Min($general, [Math]::Max(2, [Math]::Floor($pipelineCpu / 2))))
    if (($PerformanceProfile -eq "max") -and ($usableCpu -ge 12) -and (-not $script:ThreadsExplicit)) {
        $bcf = 10
    }

    $sortMemory = "512M"
    if (($PerformanceProfile -eq "max") -and ($dockerMem -ge 40GB)) {
        # Large-memory MAX profile: keep sort fully in RAM where practical.
        $sortMemory = "2048M"
    } elseif (($PerformanceProfile -eq "max") -and ($dockerMem -ge 20GB)) {
        # Home Computer (23.5 GiB WSL RAM): 4 x 1536M = ~6 GiB sort budget.
        # This is intended to remove the small coordinate-sort spill files seen
        # with 4 x 768M while retaining comfortable headroom for BWA/WSL/Windows.
        $sortMemory = "1536M"
    } elseif ($dockerMem -ge 32GB) {
        $sortMemory = "1024M"
    } elseif ($dockerMem -ge 18GB) {
        $sortMemory = "768M"
    }

    $script:DockerCpuCount = $dockerCpu
    $script:DockerMemoryBytes = $dockerMem
    $script:GeneralThreads = [int]$general
    $script:Threads = [int]$general
    $script:BwaThreads = [int]$bwa
    $script:SamtoolsSortThreads = [int]$sort
    $script:SamtoolsIndexThreads = [int]$index
    $script:BcftoolsThreads = [int]$bcf
    $script:SamtoolsSortMemory = $sortMemory

    if (-not $script:ResourceProfilePrinted) {
        $hostRamGiB = if ($script:HostMemoryBytes -gt 0) { [Math]::Round($script:HostMemoryBytes / 1GB, 1) } else { 0 }
        $dockerRamGiB = if ($dockerMem -gt 0) { [Math]::Round($dockerMem / 1GB, 1) } else { 0 }
        Write-Host ""
        Write-Host "RESOURCE PROFILE" -ForegroundColor Cyan
        Write-Host ("  Mode:                 {0}" -f $PerformanceProfile.ToUpperInvariant())
        Write-Host ("  Host:                 {0} logical CPUs, {1} GiB RAM" -f $script:LogicalCpuCount, $(if ($hostRamGiB -gt 0) { $hostRamGiB } else { "unknown" }))
        Write-Host ("  Docker/WSL:           {0} CPUs, {1} GiB RAM" -f $dockerCpu, $(if ($dockerRamGiB -gt 0) { $dockerRamGiB } else { "unknown" }))
        Write-Host ("  BWA MEM:              {0} threads" -f $script:BwaThreads)
        Write-Host ("  samtools sort:        {0} threads, -m {1} per thread" -f $script:SamtoolsSortThreads, $script:SamtoolsSortMemory)
        Write-Host ("  samtools index/VEP:   {0} threads" -f $script:GeneralThreads)
        Write-Host ("  bcftools:             {0} threads" -f $script:BcftoolsThreads)
        if (($script:HostMemoryBytes -gt 0) -and ($dockerMem -gt 0) -and ($dockerMem -lt ($script:HostMemoryBytes * 0.65))) {
            Write-Host ("  NOTE: Docker/WSL sees only {0} of {1} GiB host RAM. More WSL memory can improve sort-heavy stages." -f $dockerRamGiB, $hostRamGiB) -ForegroundColor Yellow
        }
        $script:ResourceProfilePrinted = $true
    }
}

function Get-CommandProbeTail([string]$Text, [int]$MaxChars = 1200) {
    if (-not $Text) { return "<no output>" }
    if ($Text.Length -le $MaxChars) { return $Text }
    return $Text.Substring($Text.Length - $MaxChars)
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
            if ($attempt -ge $MaxAttempts) { throw }
            Write-Host ""
            Write-Host "Network operation failed (attempt $attempt of $MaxAttempts)." -ForegroundColor Yellow
            Write-Host "Completed data will be reused; retrying in $delaySeconds seconds..." -ForegroundColor Yellow
            Start-Sleep -Seconds $delaySeconds
            $delaySeconds = [Math]::Min(60, $delaySeconds * 2)
        }
    }
}

function Test-DockerImageLocal([string]$Image) {
    $savedErrorActionPreference = $ErrorActionPreference
    try {
        $ErrorActionPreference = "Continue"
        & docker image inspect $Image --format "{{.Id}}" *> $null
        return ($LASTEXITCODE -eq 0)
    } finally {
        $ErrorActionPreference = $savedErrorActionPreference
    }
}

function Get-ApprovedDockerImages([string]$Primary, [string[]]$Fallbacks = @()) {
    return @($Primary) + @($Fallbacks)
}

function Test-ApprovedDockerImage([string]$Image, [string]$Primary, [string[]]$Fallbacks = @()) {
    return ((Get-ApprovedDockerImages $Primary $Fallbacks) -contains $Image)
}

function Resolve-DockerImage(
    [string]$Primary,
    [string[]]$Fallbacks = @(),
    [string]$Label = "Docker image"
) {
    $candidates = @(Get-ApprovedDockerImages $Primary $Fallbacks)
    foreach ($candidate in $candidates) {
        if (Test-DockerImageLocal $candidate) {
            Write-Host "Using local approved $Label image: $candidate" -ForegroundColor Green
            return $candidate
        }
    }

    $failures = @()
    foreach ($candidate in $candidates) {
        Write-Host "Approved $Label image is not present locally; trying registry source: $candidate" -ForegroundColor Yellow
        try {
            Invoke-NativeWithRetry docker @("pull", $candidate) -MaxAttempts 2 -InitialDelaySeconds 5
            if (Test-DockerImageLocal $candidate) {
                if ($candidate -ne $Primary) {
                    Write-Host "Primary registry was unavailable; using approved fallback image: $candidate" -ForegroundColor Yellow
                } else {
                    Write-Host "Downloaded approved $Label image: $candidate" -ForegroundColor Green
                }
                return $candidate
            }
            $failures += "$candidate (pull returned but image is absent)"
        } catch {
            $failures += "$candidate ($($_.Exception.Message))"
            Write-Host "Registry source failed; trying the next approved source if available." -ForegroundColor Yellow
        }
    }
    throw ("No approved $Label image is available locally and all approved registry sources failed. " +
        "Tried: $($candidates -join ', '). Completed analytical outputs and FASTQ validation are preserved. " +
        "Registry failures: $($failures -join ' | ')")
}

function Get-BwaRuntimeImage {
    if ($script:ResolvedBwaImage -and (Test-DockerImageLocal $script:ResolvedBwaImage)) { return $script:ResolvedBwaImage }
    $script:ResolvedBwaImage = Resolve-DockerImage $BwaImage $BwaFallbackImages "BWA 0.7.19"
    return $script:ResolvedBwaImage
}

function Get-SamtoolsRuntimeImage {
    if ($script:ResolvedSamtoolsImage -and (Test-DockerImageLocal $script:ResolvedSamtoolsImage)) { return $script:ResolvedSamtoolsImage }
    $script:ResolvedSamtoolsImage = Resolve-DockerImage $SamtoolsImage $SamtoolsFallbackImages "samtools 1.22.1"
    return $script:ResolvedSamtoolsImage
}

function Get-BcftoolsRuntimeImage {
    if ($script:ResolvedBcftoolsImage -and (Test-DockerImageLocal $script:ResolvedBcftoolsImage)) { return $script:ResolvedBcftoolsImage }
    $script:ResolvedBcftoolsImage = Resolve-DockerImage $BcftoolsImage $BcftoolsFallbackImages "bcftools 1.24"
    return $script:ResolvedBcftoolsImage
}

function Write-Stage15Progress([double]$Percent, [string]$Detail) {
    $bounded = [Math]::Min(100.0, [Math]::Max(0.0, $Percent))
    Write-Host ("[OVERALL] Stage 15 | {0:N2}% | {1}" -f $bounded, $Detail) -ForegroundColor DarkCyan
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

function Invoke-ObservedNative(
    [string]$Activity,
    [string]$Program,
    [string[]]$Arguments,
    [string]$WatchFile = "",
    [string]$WatchDirectory = "",
    [string]$ProgressProfile = "",
    [string]$ContainerName = "",
    [int]$HeartbeatSeconds = 15
) {
    if (-not (Test-Path -LiteralPath $ActivityRunner -PathType Leaf)) {
        throw "Activity runner is missing: $ActivityRunner"
    }
    $python = Get-Python
    $runnerArguments = @(
        $ActivityRunner,
        "--label", $Activity,
        "--heartbeat-seconds", "$HeartbeatSeconds"
    )
    if ($WatchFile) { $runnerArguments += @("--watch-file", $WatchFile) }
    if ($WatchDirectory) { $runnerArguments += @("--watch-dir", $WatchDirectory) }
    if ($ProgressProfile) { $runnerArguments += @("--progress-profile", $ProgressProfile) }
    if ($ContainerName) { $runnerArguments += @("--container-name", $ContainerName) }
    $runnerArguments += "--"
    $runnerArguments += $Program
    $runnerArguments += $Arguments
    try {
        Invoke-Native $python $runnerArguments
    } catch {
        throw ("Observed activity failed: {0}; child program: {1}. {2}" -f $Activity, $Program, $_.Exception.Message)
    }
}

function Format-ByteSize([long]$Bytes) {
    if ($Bytes -ge 1073741824) { return ("{0:N2} GiB" -f ($Bytes / 1073741824)) }
    if ($Bytes -ge 1048576) { return ("{0:N1} MiB" -f ($Bytes / 1048576)) }
    if ($Bytes -ge 1024) { return ("{0:N1} KiB" -f ($Bytes / 1024)) }
    return "$Bytes bytes"
}

function Format-Duration([double]$Seconds) {
    if ([double]::IsNaN($Seconds) -or [double]::IsInfinity($Seconds) -or $Seconds -lt 0) {
        return "unknown"
    }
    $maximumSeconds = ([TimeSpan]::MaxValue).TotalSeconds
    $span = [TimeSpan]::FromSeconds([Math]::Min($Seconds, $maximumSeconds))
    if ($span.Days -gt 0) {
        return ("{0}d {1:00}:{2:00}:{3:00}" -f $span.Days, $span.Hours, $span.Minutes, $span.Seconds)
    }
    return ("{0:00}:{1:00}:{2:00}" -f $span.Hours, $span.Minutes, $span.Seconds)
}

function Get-DownloadProgressText(
    [long]$CurrentBytes,
    [long]$TotalBytes,
    [long]$AttemptStartBytes,
    [double]$ElapsedSeconds
) {
    $gained = [Math]::Max([long]0, $CurrentBytes - $AttemptStartBytes)
    $speed = if ($ElapsedSeconds -gt 0) { $gained / $ElapsedSeconds } else { [double]0 }
    $speedText = if ($speed -gt 0) { "$(Format-ByteSize ([long][Math]::Round($speed)))/s" } else { "--/s" }
    if ($TotalBytes -gt 0) {
        $percent = [Math]::Min(100.0, [Math]::Max(0.0, (100.0 * $CurrentBytes / $TotalBytes)))
        $remaining = [Math]::Max([long]0, $TotalBytes - $CurrentBytes)
        $eta = if ($speed -gt 0) { Format-Duration ($remaining / $speed) } else { "unknown" }
        return ("Overall {0,6:N2}% | {1} / {2} | {3} | ETA {4}" -f $percent, (Format-ByteSize $CurrentBytes), (Format-ByteSize $TotalBytes), $speedText, $eta)
    }
    return ("Saved {0} | +{1} this attempt | {2} | total size unknown" -f (Format-ByteSize $CurrentBytes), (Format-ByteSize $gained), $speedText)
}

function ConvertTo-NativeCommandLineArgument([string]$Argument) {
    if ($null -eq $Argument -or $Argument.Length -eq 0) { return '""' }
    if ($Argument.Contains('"')) { throw "A native command argument contains an unsupported quote character." }
    if ($Argument -match '\s') { return ('"' + $Argument + '"') }
    return $Argument
}

function Get-RemoteContentLength([string]$Url) {
    $curl = Get-Curl
    $saved = $ErrorActionPreference
    try {
        $ErrorActionPreference = "Continue"
        $headers = @(& $curl "--silent" "--show-error" "--location" "--head" "--http1.1" "--connect-timeout" "10" "--max-time" "20" $Url 2>$null)
        if ($LASTEXITCODE -ne 0) { return [long]0 }
    } finally {
        $ErrorActionPreference = $saved
    }
    $lengths = @()
    foreach ($line in $headers) {
        $match = [regex]::Match([string]$line, '(?i)^\s*content-length:\s*(\d+)')
        if ($match.Success) { $lengths += [long]$match.Groups[1].Value }
    }
    if ($lengths.Count -eq 0) { return [long]0 }
    return [long]($lengths[-1])
}

function Get-VepCacheRemoteLength {
    if ($script:VepCacheRemoteLength -gt 0) {
        return $script:VepCacheRemoteLength
    }
    try {
        $length = Get-RemoteContentLength $VepCacheUrl
    } catch {
        # Menu rendering must remain available even when curl is missing or the
        # network is temporarily unavailable. The installer will report the
        # actionable dependency error when its item is selected.
        return [long]0
    }
    if ($length -gt 0) {
        $script:VepCacheRemoteLength = [long]$length
    }
    return [long]$length
}

function Get-FileHashTracked(
    [string]$Path,
    [ValidateSet("SHA256", "MD5")]
    [string]$Algorithm = "SHA256",
    [string]$Label = ""
) {
    $item = Get-Item -LiteralPath $Path
    $cacheKey = "$Algorithm|$($item.FullName)|$($item.Length)|$($item.LastWriteTimeUtc.Ticks)"
    if ($script:HashCache.ContainsKey($cacheKey)) {
        return [string]$script:HashCache[$cacheKey]
    }
    if (-not $Label) { $Label = "$Algorithm $($item.Name)" }
    $hashAlgorithm = if ($Algorithm -eq "MD5") {
        [System.Security.Cryptography.MD5]::Create()
    } else {
        [System.Security.Cryptography.SHA256]::Create()
    }
    $stream = $null
    $lastWidth = 0
    $startedAt = [DateTime]::UtcNow
    try {
        $stream = [System.IO.File]::Open($item.FullName, [System.IO.FileMode]::Open, [System.IO.FileAccess]::Read, [System.IO.FileShare]::Read)
        $buffer = New-Object byte[] (8 * 1024 * 1024)
        $processed = [long]0
        $nextDisplay = [DateTime]::UtcNow
        while (($read = $stream.Read($buffer, 0, $buffer.Length)) -gt 0) {
            [void]$hashAlgorithm.TransformBlock($buffer, 0, $read, $buffer, 0)
            $processed += $read
            if ([DateTime]::UtcNow -ge $nextDisplay) {
                $elapsed = ([DateTime]::UtcNow - $startedAt).TotalSeconds
                $speed = if ($elapsed -gt 0) { $processed / $elapsed } else { [double]0 }
                $percent = if ($item.Length -gt 0) { 100.0 * $processed / $item.Length } else { 100.0 }
                $remaining = [Math]::Max([long]0, $item.Length - $processed)
                $eta = if ($speed -gt 0) { Format-Duration ($remaining / $speed) } else { "unknown" }
                $text = ("Hashing {0}: {1,6:N2}% | {2} / {3} | ETA {4}" -f $Label, $percent, (Format-ByteSize $processed), (Format-ByteSize $item.Length), $eta)
                Write-Host ("`r" + $text.PadRight($lastWidth)) -NoNewline
                $lastWidth = $text.Length
                $nextDisplay = [DateTime]::UtcNow.AddSeconds(1)
            }
        }
        [void]$hashAlgorithm.TransformFinalBlock($buffer, 0, 0)
        $hash = ([System.BitConverter]::ToString($hashAlgorithm.Hash) -replace '-', '').ToLowerInvariant()
        if ($lastWidth -gt 0) {
            $done = "Hashing ${Label}: 100.00% | $(Format-ByteSize $item.Length) / $(Format-ByteSize $item.Length)"
            Write-Host ("`r" + $done.PadRight($lastWidth))
        }
        $script:HashCache[$cacheKey] = $hash
        return $hash
    } finally {
        if ($stream) { $stream.Dispose() }
        $hashAlgorithm.Dispose()
    }
}

function Get-Sha256([string]$Path) {
    return Get-FileHashTracked -Path $Path -Algorithm "SHA256"
}

function Get-StringSha256([string]$Value) {
    $algorithm = [System.Security.Cryptography.SHA256]::Create()
    try {
        $bytes = [System.Text.Encoding]::UTF8.GetBytes($Value)
        return (([System.BitConverter]::ToString($algorithm.ComputeHash($bytes))) -replace '-', '').ToLowerInvariant()
    } finally {
        $algorithm.Dispose()
    }
}

function Get-PersistentFileSha256([string]$Path) {
    $item = Get-Item -LiteralPath $Path
    $identity = "$($item.FullName)|$($item.Length)|$($item.LastWriteTimeUtc.Ticks)"
    $identityKey = Get-StringSha256 ($item.FullName.ToLowerInvariant())
    $cacheDirectory = Join-Path $Work "input_hashes"
    $cachePath = Join-Path $cacheDirectory "$identityKey.json"
    if (Test-Path -LiteralPath $cachePath -PathType Leaf) {
        try {
            $cached = Get-Content -LiteralPath $cachePath -Raw -Encoding UTF8 | ConvertFrom-Json
            if (($cached.identity -eq $identity) -and ($cached.sha256 -match '^[0-9a-f]{64}$')) {
                return [string]$cached.sha256
            }
        } catch {
            # Invalid cache metadata is ignored; the source file remains untouched.
        }
    }
    New-Item -ItemType Directory -Force -Path $cacheDirectory | Out-Null
    $hash = Get-Sha256 $item.FullName
    $value = [ordered]@{
        identity = $identity
        sha256 = $hash
    } | ConvertTo-Json -Compress
    Set-AtomicText $cachePath $value
    return $hash
}

function Set-AtomicText([string]$Path, [string]$Value) {
    $part = "$Path.part"
    Set-Content -LiteralPath $part -Value $Value -Encoding UTF8 -NoNewline
    Move-Item -LiteralPath $part -Destination $Path -Force
}

function Import-DotEnv([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) { return }
    foreach ($raw in Get-Content -LiteralPath $Path -Encoding UTF8) {
        $line = $raw.Trim()
        if (-not $line -or $line.StartsWith("#")) { continue }
        if ($line.StartsWith("export ")) { $line = $line.Substring(7).Trim() }
        $separator = $line.IndexOf("=")
        if ($separator -lt 1) { continue }
        $name = $line.Substring(0, $separator).Trim()
        $value = $line.Substring($separator + 1).Trim()
        if (($value.StartsWith('"') -and $value.EndsWith('"')) -or ($value.StartsWith("'") -and $value.EndsWith("'"))) {
            $value = $value.Substring(1, $value.Length - 2)
        }
        if ($name -match '^[A-Za-z_][A-Za-z0-9_]*$') {
            [Environment]::SetEnvironmentVariable($name, $value, "Process")
        }
    }
}

function Get-DanteSampleSexMode {
    $value = [Environment]::GetEnvironmentVariable("DANTE_SAMPLE_SEX", "Process")
    if (-not $value) { return "AUTO" }
    switch ($value.Trim().ToUpperInvariant()) {
        "AUTO"   { return "AUTO" }
        "M"      { return "M" }
        "MALE"   { return "M" }
        "XY"     { return "M" }
        "F"      { return "F" }
        "FEMALE" { return "F" }
        "XX"     { return "F" }
        default { throw "DANTE_SAMPLE_SEX must be AUTO, M/MALE/XY, or F/FEMALE/XX." }
    }
}

function Test-ApiKeyConfigured {
    $key = [Environment]::GetEnvironmentVariable("OPENAI_API_KEY", "Process")
    return [bool]($key -and $key -notmatch '^(X+|YOUR_|REPLACE_|sk-XXXX)')
}

function Invoke-ResumableDownload(
    [string]$Url,
    [string]$Destination,
    [long]$RemoteLength = 0,
    [int]$MaxConsecutiveNoProgress = 20
) {
    # One curl process per attempt is intentional. A failed attempt keeps every
    # new byte, and the next process asks for the following byte range. curl's
    # own meter is intentionally hidden because after a Range request it treats
    # only the remaining response as 100%. The meter below uses the complete
    # remote Content-Length and therefore resumes at the true overall percent.
    $curl = Get-Curl
    if ($RemoteLength -le 0) { $RemoteLength = Get-RemoteContentLength $Url }
    $attempt = 0
    $consecutiveNoProgress = 0
    while ($true) {
        $before = [long]0
        if (Test-Path -LiteralPath $Destination -PathType Leaf) {
            $before = (Get-Item -LiteralPath $Destination).Length
        }
        if ($RemoteLength -gt 0) {
            if ($before -eq $RemoteLength) {
                Write-Host "The local file already has the exact server size: $(Format-ByteSize $before)."
                Write-Host "No network retry is needed; continuing with the caller's integrity validation."
                return
            }
            if ($before -gt $RemoteLength) {
                throw "The local download is larger than the server file ($before > $RemoteLength bytes). It was preserved and will not be appended to."
            }
        }
        $attempt++
        if ($RemoteLength -gt 0) {
            $startPercent = [Math]::Min(100.0, [Math]::Max(0.0, (100.0 * $before / $RemoteLength)))
            Write-Host ("Download attempt {0}; continuing at {1:N2}% ({2} of {3})." -f $attempt, $startPercent, (Format-ByteSize $before), (Format-ByteSize $RemoteLength))
        } else {
            Write-Host "Download attempt $attempt; already saved: $(Format-ByteSize $before); total size is temporarily unknown."
        }
        $curlArguments = @(
            "--fail", "--location", "--http1.1",
            "--connect-timeout", "60",
            "--speed-limit", "1024", "--speed-time", "300",
            "--continue-at", "-",
            "--silent", "--show-error",
            "--output", $Destination,
            $Url
        )
        $argumentLine = (($curlArguments | ForEach-Object { ConvertTo-NativeCommandLineArgument ([string]$_) }) -join " ")
        $process = $null
        $startedAt = [DateTime]::UtcNow
        $lastProgressWidth = 0
        $nativeExitCode = -1
        try {
            $process = Start-Process -FilePath $curl -ArgumentList $argumentLine -PassThru -NoNewWindow
            while (-not $process.WaitForExit(1000)) {
                $current = if (Test-Path -LiteralPath $Destination -PathType Leaf) { (Get-Item -LiteralPath $Destination).Length } else { [long]0 }
                $elapsed = ([DateTime]::UtcNow - $startedAt).TotalSeconds
                $progressText = Get-DownloadProgressText $current $RemoteLength $before $elapsed
                Write-Host ("`r" + $progressText.PadRight($lastProgressWidth)) -NoNewline
                $lastProgressWidth = $progressText.Length
            }
            $process.WaitForExit()
            $nativeExitCode = $process.ExitCode
            $current = if (Test-Path -LiteralPath $Destination -PathType Leaf) { (Get-Item -LiteralPath $Destination).Length } else { [long]0 }
            $elapsed = ([DateTime]::UtcNow - $startedAt).TotalSeconds
            $progressText = Get-DownloadProgressText $current $RemoteLength $before $elapsed
            Write-Host ("`r" + $progressText.PadRight($lastProgressWidth)) -NoNewline
            $lastProgressWidth = $progressText.Length
        } finally {
            if ($process -and -not $process.HasExited) {
                $process.Kill()
                $process.WaitForExit()
            }
            if ($lastProgressWidth -gt 0) { Write-Host "" }
        }
        $after = [long]0
        if (Test-Path -LiteralPath $Destination -PathType Leaf) {
            $after = (Get-Item -LiteralPath $Destination).Length
        }

        # File state is authoritative when the server supplied Content-Length.
        # curl can return code 56 after the final response bytes were already
        # flushed to disk. In that case retrying a zero-length range forever is
        # wrong: accept the transfer as size-complete and let the caller's
        # gzip/tar/hash validation decide whether the bytes are usable.
        if ($RemoteLength -gt 0) {
            if ($after -eq $RemoteLength) {
                if ($nativeExitCode -eq 0) {
                    Write-Host "Download command completed; saved the exact server size: $(Format-ByteSize $after)"
                } else {
                    Write-Host "curl exited with code $nativeExitCode after the exact server size was saved." -ForegroundColor Yellow
                    Write-Host "The transfer is size-complete; continuing with the caller's integrity validation." -ForegroundColor Yellow
                }
                return
            }
            if ($after -gt $RemoteLength) {
                throw "The local download became larger than the server file ($after > $RemoteLength bytes). It was preserved and will not be appended to."
            }
        } elseif ($nativeExitCode -eq 0) {
            Write-Host "Download command completed; saved: $(Format-ByteSize $after)"
            return
        }

        $gained = [Math]::Max([long]0, $after - $before)
        if ($gained -gt 0) { $consecutiveNoProgress = 0 } else { $consecutiveNoProgress++ }
        Write-Host ""
        if ($nativeExitCode -eq 0) {
            Write-Host "curl stopped before the expected byte count. The partial file is preserved." -ForegroundColor Yellow
        } else {
            Write-Host "Connection interrupted (curl exit code $nativeExitCode). The partial file is preserved." -ForegroundColor Yellow
        }
        Write-Host "Saved: $(Format-ByteSize $after); added this attempt: $(Format-ByteSize $gained)." -ForegroundColor Yellow
        if ($consecutiveNoProgress -ge $MaxConsecutiveNoProgress) {
            throw "The download made no progress for $MaxConsecutiveNoProgress consecutive attempts. The partial file is intact."
        }
        $delaySeconds = [Math]::Min(120, 10 + (10 * $consecutiveNoProgress))
        Write-Host "Automatically resuming in $delaySeconds seconds. Press Ctrl+C only to pause." -ForegroundColor Yellow
        Start-Sleep -Seconds $delaySeconds
    }
}

function Assert-Docker {
    if (-not (Get-Command docker -ErrorAction SilentlyContinue)) {
        throw "Docker was not found. Install and start Docker Desktop for Windows."
    }
    $saved = $ErrorActionPreference
    try {
        $ErrorActionPreference = "Continue"
        & docker info *> $null
        $code = $LASTEXITCODE
    } finally {
        $ErrorActionPreference = $saved
    }
    if ($code -ne 0) { throw "Docker Desktop is installed, but its engine is not running." }
    Set-DockerResourceProfile
}

function Remove-VepCacheTree([string]$CacheDir) {
    if (-not (Test-Path -LiteralPath $CacheDir -PathType Container)) { return }
    Invoke-ObservedNative -Activity "Remove incomplete VEP cache extraction" -Program "docker" -Arguments @(
        "run", "--rm", "--mount", "type=bind,source=$VepData,target=/data",
        $VepImage, "bash", "-lc", "rm -rf /data/homo_sapiens/116_GRCh37"
    ) -WatchDirectory $CacheDir
}

function Find-DanteCandidates([string]$Pattern) {
    return @(Get-ChildItem -LiteralPath $PatientInput -File | Where-Object { $_.Name -match $Pattern })
}

function Find-OneDanteFile([string]$Role, [string]$Pattern) {
    # PowerShell unwraps a function's pipeline output when there is exactly one
    # match (and returns $null for zero matches). Force an array at the call
    # site so .Count and [0] behave consistently for 0, 1, and many files.
    $found = @(Find-DanteCandidates $Pattern)
    if ($found.Count -eq 0) { throw "Missing $Role VCF in patient input folder $PatientInput (expected pattern: $Pattern)." }
    if ($found.Count -gt 1) { throw "More than one $Role VCF found: $($found.Name -join ', ')" }
    return $found[0]
}

function Get-DanteInputs([string[]]$Roles) {
    $definitions = @{
        raw   = @("raw", '\.raw\.vcf(\.gz)?$')
        snp   = @("filtered SNP", '\.filtered\.snp\.vcf(\.gz)?$')
        indel = @("filtered indel", '\.filtered\.indel\.vcf(\.gz)?$')
        cnv   = @("CNV", '\.cnv\.vcf(\.gz)?$')
        sv    = @("SV", '\.sv\.vcf(\.gz)?$')
    }
    $result = @{}
    foreach ($role in $Roles) {
        $result[$role] = Find-OneDanteFile $definitions[$role][0] $definitions[$role][1]
    }
    return $result
}

function Get-FastqPairInfo([bool]$FullGzipCheck = $false) {
    $handler = Join-Path $Root "fastq_pairing.py"
    if (-not (Test-Path -LiteralPath $handler -PathType Leaf)) {
        return [pscustomobject]@{ status = "HANDLER_MISSING"; r1 = $null; r2 = $null; details = @("fastq_pairing.py is missing") }
    }
    try { $python = Get-Python } catch {
        return [pscustomobject]@{ status = "PYTHON_MISSING"; r1 = $null; r2 = $null; details = @($_.Exception.Message) }
    }
    $arguments = @($handler, "--root", $PatientInput)
    if ($FullGzipCheck) { $arguments += "--full-gzip-check" }
    $saved = $ErrorActionPreference
    try {
        $ErrorActionPreference = "Continue"
        $output = @(& $python @arguments 2>&1)
        $code = $LASTEXITCODE
    } finally {
        $ErrorActionPreference = $saved
    }
    if ($code -ne 0) {
        return [pscustomobject]@{ status = "HANDLER_ERROR"; r1 = $null; r2 = $null; details = @($output -join "`n") }
    }
    try {
        return (($output -join "`n") | ConvertFrom-Json)
    } catch {
        return [pscustomobject]@{ status = "HANDLER_ERROR"; r1 = $null; r2 = $null; details = @("FASTQ handler returned invalid JSON") }
    }
}

function Get-FastqPairMenuStatus($Pair = $null) {
    $pair = $Pair
    if ($null -eq $pair) { $pair = Get-FastqPairInfo }
    switch ([string]$pair.status) {
        "READY" { return "READY" }
        "SINGLE_END_READY" { return "READY: SINGLE-END R1" }
        "R1_MISSING" { return "R1 MISSING" }
        "MISSING" { return "MISSING" }
        "MULTIPLE_PAIRS" { return "MULTIPLE PAIRS" }
        "PAIRING_ERROR" { return "PAIRING ERROR" }
        "INVALID" { return "INVALID" }
        default { return [string]$pair.status }
    }
}

function Test-FastqInputReady($Pair) {
    return (@("READY", "SINGLE_END_READY") -contains [string]$Pair.status)
}

function Test-RelocatedLegacyFastqIdentity([string]$CachedIdentity, $R1, $R2, [string]$ReadMode) {
    if (-not $CachedIdentity -or -not $CachedIdentity.StartsWith("fastq-input-v2|")) { return $false }
    $parts = @($CachedIdentity -split '\|')
    if ($ReadMode -eq "single") {
        if ($parts.Count -ne 6 -or $parts[1] -ne "single" -or $parts[5] -ne "NO_R2") { return $false }
        $cachedR1Name = Split-Path -Leaf ([string]$parts[2])
        return (
            ($cachedR1Name -eq $R1.Name) -and
            ([long]$parts[3] -eq [long]$R1.Length) -and
            ([long]$parts[4] -eq [long]$R1.LastWriteTimeUtc.Ticks)
        )
    }
    if ($ReadMode -eq "paired") {
        if ($parts.Count -ne 8 -or $parts[1] -ne "paired" -or $null -eq $R2) { return $false }
        $cachedR1Name = Split-Path -Leaf ([string]$parts[2])
        $cachedR2Name = Split-Path -Leaf ([string]$parts[5])
        return (
            ($cachedR1Name -eq $R1.Name) -and
            ([long]$parts[3] -eq [long]$R1.Length) -and
            ([long]$parts[4] -eq [long]$R1.LastWriteTimeUtc.Ticks) -and
            ($cachedR2Name -eq $R2.Name) -and
            ([long]$parts[6] -eq [long]$R2.Length) -and
            ([long]$parts[7] -eq [long]$R2.LastWriteTimeUtc.Ticks)
        )
    }
    return $false
}

function Confirm-FastqInputComplete([int]$ProgressStage = 15) {
    New-Item -ItemType Directory -Force -Path $Work | Out-Null
    $pair = Get-FastqPairInfo
    if (-not (Test-FastqInputReady $pair)) {
        $detail = @($pair.details) -join "; "
        throw "One valid R1 (single-end) or one matching R1/R2 pair is required in patient input folder $PatientInput. Status: $($pair.status). $detail"
    }
    $readMode = [string]$pair.read_mode
    if (@("single", "paired") -notcontains $readMode) { throw "The FASTQ handler returned an unsupported read mode: $readMode" }
    $r1 = Get-Item -LiteralPath ([string]$pair.r1)
    $r2 = $null
    if ($readMode -eq "paired") { $r2 = Get-Item -LiteralPath ([string]$pair.r2) }
    $identity = if ($readMode -eq "paired") {
        "fastq-input-v3|paired|$($r1.Name)|$($r1.Length)|$($r1.LastWriteTimeUtc.Ticks)|$($r2.Name)|$($r2.Length)|$($r2.LastWriteTimeUtc.Ticks)"
    } else {
        "fastq-input-v3|single|$($r1.Name)|$($r1.Length)|$($r1.LastWriteTimeUtc.Ticks)|NO_R2"
    }
    $state = Join-Path $Work "Dante_FASTQ_pair_validation.state"
    if (Test-NonEmptyFile $state) {
        try {
            $cached = Get-Content -LiteralPath $state -Raw -Encoding UTF8 | ConvertFrom-Json
            $cachedHashesReady = (
                ($cached.r1_content_sha256 -match '^[0-9a-f]{64}$') -and
                (($readMode -eq "single") -or ($cached.r2_content_sha256 -match '^[0-9a-f]{64}$'))
            )
            $cachedCountsReady = (
                ([long]$cached.r1_record_count -gt 0) -and
                (($readMode -eq "single") -or ([long]$cached.r2_record_count -gt 0))
            )
            $identityCompatible = (
                ($cached.identity -eq $identity) -or
                (Test-RelocatedLegacyFastqIdentity ([string]$cached.identity) $r1 $r2 $readMode)
            )
            if ($identityCompatible -and ([string]$cached.read_mode -eq $readMode) -and $cachedHashesReady -and $cachedCountsReady) {
                if ([string]$cached.identity -ne $identity) {
                    $cached.identity = $identity
                    Set-AtomicText $state ($cached | ConvertTo-Json -Depth 6)
                    Write-Host "Adopted relocated FASTQ validation state without rereading the FASTQ stream." -ForegroundColor DarkGray
                }
                $pair | Add-Member -NotePropertyName r1_content_sha256 -NotePropertyValue ([string]$cached.r1_content_sha256) -Force
                $pair | Add-Member -NotePropertyName r1_record_count -NotePropertyValue ([long]$cached.r1_record_count) -Force
                $r2Hash = if ($readMode -eq "paired") { [string]$cached.r2_content_sha256 } else { $null }
                $r2Records = if ($readMode -eq "paired") { [long]$cached.r2_record_count } else { [long]0 }
                $pair | Add-Member -NotePropertyName r2_content_sha256 -NotePropertyValue $r2Hash -Force
                $pair | Add-Member -NotePropertyName r2_record_count -NotePropertyValue $r2Records -Force
                Write-Host "The FASTQ gzip/record validation is current for $readMode-end input; reusing it."
                if ($ProgressStage -eq 15) {
                    Write-Stage15Progress 5.0 "FASTQ validation cached; exact read counts available"
                } else {
                    Write-Host ("[OVERALL] Stage {0} | 5.00% | FASTQ validation cached; exact read counts available" -f $ProgressStage) -ForegroundColor DarkCyan
                }
                return $pair
            } elseif ($identityCompatible -and $cachedHashesReady) {
                Write-Host "Existing FASTQ validation predates exact read-count progress metadata; upgrading it once." -ForegroundColor Yellow
            }
        } catch {
            # Invalid validation metadata is ignored and rebuilt from the immutable FASTQ inputs.
        }
    }
    Write-Step "Validate the complete FASTQ input streams"
    New-Item -ItemType Directory -Force -Path $Qc | Out-Null
    $validationReport = Join-Path $Qc "Dante_FASTQ_pair_validation.json"
    if (Test-Path -LiteralPath $validationReport) { Remove-Item -LiteralPath $validationReport -Force }
    Invoke-ObservedNative -Activity "Validate and fingerprint complete FASTQ input" -Program (Get-Python) -Arguments @(
        (Join-Path $Root "fastq_pairing.py"),
        "--root", $PatientInput,
        "--full-gzip-check",
        "--output-json", $validationReport,
        "--progress-stage", "$ProgressStage",
        "--progress-start", "0",
        "--progress-span", "5"
    ) -WatchFile "$validationReport.part"
    if (-not (Test-NonEmptyFile $validationReport)) {
        throw "The FASTQ integrity handler did not create its validation report."
    }
    $checked = Get-Content -LiteralPath $validationReport -Raw -Encoding UTF8 | ConvertFrom-Json
    if (-not (Test-FastqInputReady $checked)) {
        $detail = @($checked.details) -join "; "
        throw "The FASTQ integrity check failed. Status: $($checked.status). $detail"
    }
    $checkedMode = [string]$checked.read_mode
    if (($checked.r1_content_sha256 -notmatch '^[0-9a-f]{64}$') -or
        (($checkedMode -eq "paired") -and ($checked.r2_content_sha256 -notmatch '^[0-9a-f]{64}$'))) {
        throw "The FASTQ integrity checker did not return all required content fingerprints."
    }
    if ([long]$checked.r1_record_count -le 0 -or (($checkedMode -eq "paired") -and ([long]$checked.r2_record_count -le 0))) {
        throw "The FASTQ integrity checker did not return exact FASTQ record counts required for progress reporting."
    }
    $validationState = [ordered]@{
        schema = "fastq-validation-v3"
        identity = $identity
        read_mode = $checkedMode
        r1_content_sha256 = [string]$checked.r1_content_sha256
        r2_content_sha256 = $(if ($checkedMode -eq "paired") { [string]$checked.r2_content_sha256 } else { $null })
        r1_record_count = [long]$checked.r1_record_count
        r2_record_count = $(if ($checkedMode -eq "paired") { [long]$checked.r2_record_count } else { [long]0 })
        r1_decompressed_bytes = [long]$checked.r1_decompressed_bytes
        r2_decompressed_bytes = $(if ($checkedMode -eq "paired") { [long]$checked.r2_decompressed_bytes } else { [long]0 })
    } | ConvertTo-Json
    Set-AtomicText $state $validationState
    return $checked
}

function Test-VepCacheReady {
    $cacheDir = Join-Path $VepData "homo_sapiens\116_GRCh37"
    $structureReady = (
        (Test-Path -LiteralPath (Join-Path $cacheDir "info.txt") -PathType Leaf) -and
        (Test-Path -LiteralPath (Join-Path $cacheDir "1") -PathType Container)
    )
    if (-not $structureReady) { return $false }
    $marker = Join-Path $VepData ".dante_vep_116_grch37.complete"
    if (Test-Path -LiteralPath $marker -PathType Leaf) { return $true }
    # v7.3 removed the source archive only after tar returned successfully.
    # This safely adopts such a legacy installation.  If the archive is still
    # present, the directory may be a partial extraction and is not READY.
    $archive = Join-Path $VepData $VepCacheArchiveName
    return (-not (Test-Path -LiteralPath $archive -PathType Leaf))
}

function Test-FastaReady {
    $path = Join-Path $VepFastaData $FastaName
    return ((Test-Path -LiteralPath $path -PathType Leaf) -and ((Get-Item -LiteralPath $path).Length -gt 2147483648))
}

function Test-ClinVarReady {
    $vcf = Join-Path $ClinVarData "clinvar.vcf.gz"
    $index = "$vcf.tbi"
    return (
        (Test-Path -LiteralPath $vcf -PathType Leaf) -and
        ((Get-Item -LiteralPath $vcf).Length -gt 104857600) -and
        (Test-Path -LiteralPath $index -PathType Leaf)
    )
}

function Get-AnnotSVCompletionMarker {
    return (Join-Path $AnnotSVData $AnnotSVCompleteMarkerName)
}

function Test-AnnotSVHumanStructureReady {
    $genesDirectory = Join-Path $AnnotSVData "AnnotSV_annotations\Annotations_Human\Genes\GRCh37"
    if (-not (Test-Path -LiteralPath $genesDirectory -PathType Container)) { return $false }
    # AnnotSV 3.5.10 accepts either the downloaded UCSC source or the already
    # formatted RefSeq BED. The current official 3.5 human archive can contain
    # only the formatted files, so requiring a wildcard *.txt.gz rejects a
    # valid, fully extracted installation.
    $acceptedRefSeqFiles = @(
        "genes.RefSeq.sorted.bed",
        "ncbiRefSeq.txt.gz"
    )
    foreach ($name in $acceptedRefSeqFiles) {
        $path = Join-Path $genesDirectory $name
        if ((Test-Path -LiteralPath $path -PathType Leaf) -and
            ((Get-Item -LiteralPath $path).Length -gt 0)) {
            return $true
        }
    }
    return $false
}

function Test-AnnotSVReady {
    return (
        (Test-AnnotSVHumanStructureReady) -and
        (Test-Path -LiteralPath (Get-AnnotSVCompletionMarker) -PathType Leaf)
    )
}

function Get-AnnotSVMenuStatus {
    if (Test-AnnotSVReady) { return "READY" }
    $annotationDirectory = Join-Path $AnnotSVData "AnnotSV_annotations"
    if (Test-Path -LiteralPath $annotationDirectory -PathType Container) {
        return "INCOMPLETE; REPAIR PENDING"
    }
    return "MISSING"
}

function Test-AnnotSVOutputReady([string]$Path, [string]$ReportPath = "") {
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) { return $false }
    if ((Get-Item -LiteralPath $Path).Length -le 0) { return $false }
    try {
        $lines = @(Get-Content -LiteralPath $Path -Encoding UTF8 -TotalCount 10)
        if ($lines.Count -lt 1 -or -not $lines[0]) { return $false }
        $columns = @($lines[0] -split "`t")
        foreach ($required in @("AnnotSV_ID", "SV_chrom", "SV_start", "SV_end", "SV_type")) {
            if ($columns -notcontains $required) { return $false }
        }
        $dataLines = @($lines | Select-Object -Skip 1 | Where-Object { $_ -and $_.Trim() })
        if ($dataLines.Count -gt 0) { return $true }

        # A header-only table is complete only when preparation selected zero
        # input rows. Otherwise it is the characteristic interrupted AnnotSV
        # file that used to be mistaken for a completed result.
        if (Test-Path -LiteralPath $ReportPath -PathType Leaf) {
            $report = Get-Content -LiteralPath $ReportPath -Raw -Encoding UTF8
            return [bool]($report -match '(?m)^Selected rows:\s*0\s*$')
        }
    } catch {
        return $false
    }
    return $false
}

function Test-AnnotSVOutputCurrent(
    [string]$OutputPath,
    [string]$ReportPath,
    [System.IO.FileInfo]$InputFile
) {
    if (-not (Test-AnnotSVOutputReady $OutputPath $ReportPath)) { return $false }
    $outputItem = Get-Item -LiteralPath $OutputPath
    if ($outputItem.LastWriteTimeUtc -lt $InputFile.LastWriteTimeUtc) { return $false }
    $marker = Get-AnnotSVCompletionMarker
    if ((Test-Path -LiteralPath $marker -PathType Leaf) -and
        ($outputItem.LastWriteTimeUtc -lt (Get-Item -LiteralPath $marker).LastWriteTimeUtc)) {
        return $false
    }
    return $true
}

function Get-AnnotSVLayerMenuState([ValidateSet("cnv", "sv")][string]$Kind) {
    $upper = $Kind.ToUpper()
    $output = Join-Path $Results "Dante_$($upper)_AnnotSV.tsv"
    $part = Join-Path $Results "Dante_$($upper)_AnnotSV.part.tsv"
    $report = Join-Path $Qc "Dante_$($upper)_AnnotSV_input_report.txt"
    if (Test-AnnotSVOutputReady $output $report) {
        $outputItem = Get-Item -LiteralPath $output
        $pattern = if ($Kind -eq "cnv") { '\.cnv\.vcf(\.gz)?$' } else { '\.sv\.vcf(\.gz)?$' }
        $inputs = @(Find-DanteCandidates $pattern)
        if (($inputs.Count -eq 1) -and ($outputItem.LastWriteTimeUtc -lt $inputs[0].LastWriteTimeUtc)) {
            return "STALE"
        }
        $marker = Get-AnnotSVCompletionMarker
        if ((Test-Path -LiteralPath $marker -PathType Leaf) -and
            ($outputItem.LastWriteTimeUtc -lt (Get-Item -LiteralPath $marker).LastWriteTimeUtc)) {
            return "STALE"
        }
        return "READY"
    }
    if (Test-Path -LiteralPath $output -PathType Leaf) { return "INVALID" }
    if (Test-Path -LiteralPath $part -PathType Leaf) { return "INTERRUPTED" }
    return "NOT RUN"
}

function Get-AnnotSVOutputMenuStatus {
    $cnv = Get-AnnotSVLayerMenuState "cnv"
    $sv = Get-AnnotSVLayerMenuState "sv"
    if (($cnv -eq "READY") -and ($sv -eq "READY")) { return "READY" }
    $details = "CNV $cnv; SV $sv"
    if (($cnv -eq "INVALID") -or ($sv -eq "INVALID")) { return "INCOMPLETE: $details" }
    if (($cnv -eq "STALE") -and ($sv -eq "STALE")) { return "STALE: $details" }
    if (($cnv -eq "READY") -or ($sv -eq "READY") -or
        ($cnv -eq "STALE") -or ($sv -eq "STALE")) { return "PARTIAL: $details" }
    if (($cnv -eq "INTERRUPTED") -or ($sv -eq "INTERRUPTED")) { return "INTERRUPTED: $details" }
    return "NOT RUN"
}

function Test-GzipInVepData([string]$ArchiveName) {
    $archivePath = Join-Path $VepData $ArchiveName
    if (-not (Test-Path -LiteralPath $archivePath -PathType Leaf)) { return $false }
    try {
        Invoke-ObservedNative -Activity "Verify gzip integrity: $ArchiveName" -Program "docker" -Arguments @(
            "run", "--rm", "--mount", "type=bind,source=$VepData,target=/data",
            $VepImage, "bash", "-lc", "gzip -t /data/$ArchiveName"
        ) -WatchFile $archivePath
        return $true
    } catch {
        return $false
    }
}

function Test-GzipWithPython([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) { return $false }
    $python = Get-Python
    try {
        Invoke-ObservedNative -Activity "Verify gzip integrity: $(Split-Path -Leaf $Path)" -Program $python -Arguments @(
            "-c", "import gzip,sys; f=gzip.open(sys.argv[1],'rb'); [None for _ in iter(lambda:f.read(8*1024*1024),b'')]; f.close()", $Path
        ) -WatchFile $Path
        return $true
    } catch {
        return $false
    }
}

function Get-VepCacheMenuStatus {
    if (Test-VepCacheReady) { return "READY" }
    $archive = Join-Path $VepData $VepCacheArchiveName
    if (Test-Path -LiteralPath $archive -PathType Leaf) {
        $localLength = [long](Get-Item -LiteralPath $archive).Length
        $remoteLength = Get-VepCacheRemoteLength
        if ($remoteLength -gt 0) {
            if ($localLength -eq $remoteLength) {
                return "DOWNLOADED: $(Format-ByteSize $localLength); VERIFY/EXTRACT PENDING"
            }
            if ($localLength -lt $remoteLength) {
                $percent = 100.0 * $localLength / $remoteLength
                return ("PARTIAL: {0} / {1} ({2:N2}%)" -f (Format-ByteSize $localLength), (Format-ByteSize $remoteLength), $percent)
            }
            return "SIZE MISMATCH: local $(Format-ByteSize $localLength) > server $(Format-ByteSize $remoteLength)"
        }
        return "ARCHIVE PRESENT: $(Format-ByteSize $localLength); SIZE CHECK UNAVAILABLE"
    }
    return "MISSING"
}

function Get-RsidExportMenuStatus {
    $export = Join-Path $Results "Dante_23andMe_style.txt.gz"
    $report = Join-Path $Qc "Dante_23andMe_style_report.txt"
    $state = "$export.state.sha256"
    $vepOutput = Join-Path $Work "Dante_SNP_rsid_lookup.vep.vcf.gz"
    $vepPart = Join-Path $Work "Dante_SNP_rsid_lookup.part.vep.vcf.gz"
    if ((Test-Path -LiteralPath $export -PathType Leaf) -and
        (Test-Path -LiteralPath $report -PathType Leaf) -and
        (Test-Path -LiteralPath $state -PathType Leaf)) { return "READY" }
    if ((Test-Path -LiteralPath $export -PathType Leaf) -and
        (Test-Path -LiteralPath $report -PathType Leaf)) { return "LEGACY OUTPUT PRESENT" }
    if (Test-Path -LiteralPath $vepOutput -PathType Leaf) { return "VEP FILE PRESENT; VERIFY/EXPORT PENDING" }
    if (Test-Path -LiteralPath $vepPart -PathType Leaf) {
        return "INTERRUPTED VEP: $(Format-ByteSize ((Get-Item -LiteralPath $vepPart).Length)); RESTART REQUIRED"
    }
    return "NOT RUN"
}

function Test-NonEmptyFile([string]$Path) {
    return (
        (Test-Path -LiteralPath $Path -PathType Leaf) -and
        ((Get-Item -LiteralPath $Path).Length -gt 0)
    )
}

function Test-FileNotOlderThan([string]$Path, [string[]]$Dependencies) {
    if (-not (Test-NonEmptyFile $Path)) { return $false }
    $outputTime = (Get-Item -LiteralPath $Path).LastWriteTimeUtc
    foreach ($dependency in $Dependencies) {
        if (-not (Test-Path -LiteralPath $dependency -PathType Leaf)) { return $false }
        if ($outputTime -lt (Get-Item -LiteralPath $dependency).LastWriteTimeUtc) { return $false }
    }
    return $true
}

function Get-ArtifactMenuStatus([bool]$Ready, [string[]]$Paths) {
    if ($Ready) { return "READY" }
    foreach ($path in $Paths) {
        if (Test-Path -LiteralPath $path) { return "INCOMPLETE" }
    }
    return "NOT RUN"
}

function Test-NormalizedLayersReady {
    if (-not (Test-FastaReady)) { return $false }
    $fasta = Join-Path $VepFastaData $FastaName
    $referenceTargets = Join-Path $Qc $ReferenceContigTargetsName
    $referenceTargetsState = "$referenceTargets.state.sha256"
    $contigAudit = Join-Path $Qc $ReferenceContigAuditName
    foreach ($required in @($referenceTargets, $referenceTargetsState, $contigAudit)) {
        if (-not (Test-NonEmptyFile $required)) { return $false }
    }
    if (-not (Test-FileNotOlderThan -Path $referenceTargets -Dependencies @($fasta))) {
        return $false
    }
    $definitions = @(
        @("Dante_SNP_normalized.vcf.gz", '\.filtered\.snp\.vcf(\.gz)?$'),
        @("Dante_INDEL_normalized.vcf.gz", '\.filtered\.indel\.vcf(\.gz)?$'),
        @("Dante_RAW_normalized.vcf.gz", '\.raw\.vcf(\.gz)?$')
    )
    $normalized = @()
    $inputDependencies = @()
    foreach ($definition in $definitions) {
        $output = Join-Path $Work $definition[0]
        $normalized += $output
        $inputs = @(Find-DanteCandidates $definition[1])
        if ($inputs.Count -ne 1) { return $false }
        $inputDependencies += $inputs[0].FullName
        foreach ($required in @($output, "$output.tbi", "$output.state.sha256")) {
            if (-not (Test-NonEmptyFile $required)) { return $false }
        }
        if (-not (Test-FileNotOlderThan -Path $output -Dependencies @($inputs[0].FullName, $fasta, $referenceTargets))) {
            return $false
        }
    }

    $auditDependencies = @($referenceTargets) + $inputDependencies
    if (-not (Test-FileNotOlderThan -Path $contigAudit -Dependencies $auditDependencies)) {
        return $false
    }

    $residual = Join-Path $Work "Dante_RESIDUAL_normalized.vcf.gz"
    $residualReport = Join-Path $Qc "Dante_RESIDUAL_preparation_report.txt"
    foreach ($required in @($residual, "$residual.state.sha256", $residualReport)) {
        if (-not (Test-NonEmptyFile $required)) { return $false }
    }
    return (Test-FileNotOlderThan -Path $residual -Dependencies $normalized)
}

function Get-NormalizedLayersMenuStatus {
    $paths = @(
        (Join-Path $Work "Dante_SNP_normalized.vcf.gz"),
        (Join-Path $Work "Dante_SNP_normalized.part.vcf.gz"),
        (Join-Path $Work "Dante_INDEL_normalized.vcf.gz"),
        (Join-Path $Work "Dante_INDEL_normalized.part.vcf.gz"),
        (Join-Path $Work "Dante_RAW_normalized.vcf.gz"),
        (Join-Path $Work "Dante_RAW_normalized.part.vcf.gz"),
        (Join-Path $Work "Dante_RESIDUAL_normalized.vcf.gz"),
        (Join-Path $Work "Dante_RESIDUAL_normalized.vcf.gz.part"),
        (Join-Path $Qc $ReferenceContigTargetsName),
        (Join-Path $Qc $ReferenceContigAuditName)
    )
    return (Get-ArtifactMenuStatus -Ready (Test-NormalizedLayersReady) -Paths $paths)
}

function Test-FullVepOutputsReady {
    if (-not (Test-NormalizedLayersReady)) { return $false }
    foreach ($layer in @("SNP", "INDEL", "RESIDUAL")) {
        $input = Join-Path $Work "Dante_$($layer)_normalized.vcf.gz"
        $vep = Join-Path $FullVepResults "Dante_$($layer)_full.vep.vcf.gz"
        $summary = Join-Path $FullVepResults "Dante_$($layer)_VEP_summary.tsv.gz"
        $rsid = Join-Path $FullVepResults "Dante_$($layer)_rsid.vcf.gz"
        $report = Join-Path $Qc "Dante_$($layer)_VEP_report.txt"
        foreach ($required in @($vep, "$vep.state.sha256", $summary, $rsid, $report)) {
            if (-not (Test-NonEmptyFile $required)) { return $false }
        }
        if (-not (Test-FileNotOlderThan -Path $vep -Dependencies @($input))) { return $false }
        foreach ($postprocessed in @($summary, $rsid, $report)) {
            if (-not (Test-FileNotOlderThan -Path $postprocessed -Dependencies @($vep))) { return $false }
        }
    }
    return $true
}

function Get-FullVepMenuStatus {
    $paths = @(
        (Join-Path $FullVepResults "Dante_SNP_full.vep.vcf.gz"),
        (Join-Path $FullVepResults "Dante_SNP_full.part.vep.vcf.gz"),
        (Join-Path $FullVepResults "Dante_INDEL_full.vep.vcf.gz"),
        (Join-Path $FullVepResults "Dante_INDEL_full.part.vep.vcf.gz"),
        (Join-Path $FullVepResults "Dante_RESIDUAL_full.vep.vcf.gz"),
        (Join-Path $FullVepResults "Dante_RESIDUAL_full.part.vep.vcf.gz")
    )
    return (Get-ArtifactMenuStatus -Ready (Test-FullVepOutputsReady) -Paths $paths)
}

function Test-UnifiedCatalogReady {
    if (-not (Test-FullVepOutputsReady) -or -not (Test-ClinVarReady)) { return $false }
    $catalog = Join-Path $Results "Dante_small_variants_annotated.tsv.gz"
    $report = Join-Path $Qc "Dante_small_variants_catalog_report.txt"
    $manifestPath = Join-Path $Results "Dante_pipeline_manifest.json"
    foreach ($required in @($catalog, $report, $manifestPath)) {
        if (-not (Test-NonEmptyFile $required)) { return $false }
    }
    $dependencies = @(
        (Join-Path $FullVepResults "Dante_SNP_VEP_summary.tsv.gz"),
        (Join-Path $FullVepResults "Dante_INDEL_VEP_summary.tsv.gz"),
        (Join-Path $FullVepResults "Dante_RESIDUAL_VEP_summary.tsv.gz"),
        (Join-Path $ClinVarData "clinvar.vcf.gz")
    )
    try {
        $manifest = Get-Content -LiteralPath $manifestPath -Raw -Encoding UTF8 | ConvertFrom-Json
        if (($manifest.catalog.status -ne "complete") -or
            ($manifest.catalog.file -ne (Split-Path -Leaf $catalog)) -or
            ($manifest.catalog.sha256 -notmatch '^[0-9a-f]{64}$')) {
            return $false
        }
        $summaryIdentities = @($manifest.inputs.vep_summaries)
        if ($summaryIdentities.Count -ne 3) { return $false }
        foreach ($dependency in $dependencies[0..2]) {
            $name = Split-Path -Leaf $dependency
            $matches = @($summaryIdentities | Where-Object { $_.name -eq $name })
            if (($matches.Count -ne 1) -or
                ([long]$matches[0].size -ne [long](Get-Item -LiteralPath $dependency).Length)) {
                return $false
            }
        }
        $clinvar = $dependencies[3]
        if (($manifest.inputs.clinvar.name -ne (Split-Path -Leaf $clinvar)) -or
            ([long]$manifest.inputs.clinvar.size -ne [long](Get-Item -LiteralPath $clinvar).Length)) {
            return $false
        }
        return $true
    } catch {
        return $false
    }
}

function Get-UnifiedCatalogMenuStatus {
    $paths = @(
        (Join-Path $Results "Dante_small_variants_annotated.tsv.gz"),
        (Join-Path $Results "Dante_small_variants_annotated.tsv.gz.part"),
        (Join-Path $Results "Dante_pipeline_manifest.json")
    )
    return (Get-ArtifactMenuStatus -Ready (Test-UnifiedCatalogReady) -Paths $paths)
}

function Test-ShortlistReady {
    if (-not (Test-UnifiedCatalogReady)) { return $false }
    $shortlist = Join-Path $Results "Dante_priority_candidates.tsv"
    $genomewide = Join-Path $Results "Dante_genomewide_uncertain_coding.tsv"
    $report = Join-Path $Qc "Dante_priority_candidates_report.txt"
    $statePath = "$shortlist.state.json"
    $geneReference = Join-Path $Root "reference_data\CHESS_2.2_genes.tsv"
    foreach ($required in @($shortlist, $genomewide, $report, $statePath, $geneReference)) {
        if (-not (Test-NonEmptyFile $required)) { return $false }
    }
    try {
        $state = Get-Content -LiteralPath $statePath -Raw -Encoding UTF8 | ConvertFrom-Json
        $manifest = Get-Content -LiteralPath (Join-Path $Results "Dante_pipeline_manifest.json") -Raw -Encoding UTF8 | ConvertFrom-Json
        $geneHash = (Get-FileHash -LiteralPath $geneReference -Algorithm SHA256).Hash.ToLowerInvariant()
        return (
            ($state.status -eq "complete") -and
            ($state.rules_schema -eq "deterministic-shortlist-v3-chess-genomewide") -and
            ($state.catalog_sha256 -match '^[0-9a-f]{64}$') -and
            ($state.catalog_sha256 -eq $manifest.catalog.sha256) -and
            ([string]$state.gene_reference_sha256 -eq $geneHash) -and
            ([string]$state.genomewide_output -eq "Dante_genomewide_uncertain_coding.tsv")
        )
    } catch {
        return $false
    }
}

function Get-ShortlistMenuStatus {
    $paths = @(
        (Join-Path $Results "Dante_priority_candidates.tsv"),
        (Join-Path $Results "Dante_priority_candidates.tsv.part"),
        (Join-Path $Results "Dante_priority_candidates.tsv.state.json"),
        (Join-Path $Results "Dante_genomewide_uncertain_coding.tsv"),
        (Join-Path $Results "Dante_genomewide_uncertain_coding.tsv.part")
    )
    return (Get-ArtifactMenuStatus -Ready (Test-ShortlistReady) -Paths $paths)
}

function Get-FileIdentity([string]$Path) {
    $item = Get-Item -LiteralPath $Path
    return [ordered]@{
        name = $item.Name
        length = [long]$item.Length
        modified_utc_ticks = [long]$item.LastWriteTimeUtc.Ticks
    }
}

function Test-StateFileIdentity($Identity, [string]$Path) {
    if (($null -eq $Identity) -or -not (Test-Path -LiteralPath $Path -PathType Leaf)) { return $false }
    $item = Get-Item -LiteralPath $Path
    return (
        ([string]$Identity.name -eq $item.Name) -and
        ([long]$Identity.length -eq [long]$item.Length) -and
        ([long]$Identity.modified_utc_ticks -eq [long]$item.LastWriteTimeUtc.Ticks)
    )
}


function Test-StateFileContentIdentity($Identity, [string]$Path) {
    if (($null -eq $Identity) -or -not (Test-Path -LiteralPath $Path -PathType Leaf)) { return $false }
    $item = Get-Item -LiteralPath $Path
    if (([string]$Identity.name -ne $item.Name) -or
        ([long]$Identity.length -ne [long]$item.Length) -or
        ([string]$Identity.sha256 -notmatch '^[0-9a-f]{64}$')) {
        return $false
    }
    return ((Get-Sha256 $Path) -eq ([string]$Identity.sha256).ToLowerInvariant())
}

function Get-ReferencePanelFiles {
    if (-not (Test-Path -LiteralPath $ReferencePanelSources -PathType Container)) { return @() }
    return @(Get-ChildItem -LiteralPath $ReferencePanelSources -File | Where-Object {
        $_.Name -match '(?i)\.(txt|tsv|csv)(\.gz)?$'
    } | Sort-Object Name)
}

function Test-ReferenceTargetsReady {
    if (-not (Test-ClinVarReady)) { return $false }
    $targets = Join-Path $Qc $ReferenceTargetsName
    $bed = Join-Path $Qc $ReferenceTargetsBedName
    $report = Join-Path $Qc $ReferenceTargetsReportName
    $coordinateConflicts = Join-Path $Qc $ReferenceCoordinateConflictsName
    $genotypeConflicts = Join-Path $Qc $ReferenceSourceGenotypeConflictsName
    $statePath = "$targets.state.json"
    foreach ($required in @($targets, $bed, $report, $coordinateConflicts, $genotypeConflicts, $statePath)) {
        if (-not (Test-NonEmptyFile $required)) { return $false }
    }
    try {
        $state = Get-Content -LiteralPath $statePath -Raw -Encoding UTF8 | ConvertFrom-Json
        if (($state.status -ne "complete") -or ($state.schema -ne "reference-targets-v2")) { return $false }
        if (-not (Test-StateFileIdentity $state.inputs.clinvar (Join-Path $ClinVarData "clinvar.vcf.gz"))) { return $false }
        $currentPanels = @(Get-ReferencePanelFiles)
        $savedPanels = @($state.inputs.panels)
        if ($currentPanels.Count -ne $savedPanels.Count) { return $false }
        foreach ($panel in $currentPanels) {
            $matches = @($savedPanels | Where-Object { $_.name -eq $panel.Name })
            if (($matches.Count -ne 1) -or -not (Test-StateFileIdentity $matches[0] $panel.FullName)) { return $false }
        }
        return (
            (Test-StateFileIdentity $state.outputs.targets $targets) -and
            (Test-StateFileIdentity $state.outputs.bed $bed) -and
            (Test-StateFileIdentity $state.outputs.coordinate_conflicts $coordinateConflicts) -and
            (Test-StateFileIdentity $state.outputs.genotype_conflicts $genotypeConflicts) -and
            ([string]$state.outputs.targets.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.outputs.bed.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.outputs.coordinate_conflicts.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.outputs.genotype_conflicts.sha256 -match '^[0-9a-f]{64}$')
        )
    } catch { return $false }
}

function Get-ReferenceTargetsMenuStatus {
    if (-not (Test-ClinVarReady)) { return "BLOCKED: CLINVAR MISSING" }
    $paths = @(
        (Join-Path $Qc $ReferenceTargetsName),
        (Join-Path $Qc "$ReferenceTargetsName.part"),
        (Join-Path $Qc $ReferenceTargetsBedName),
        (Join-Path $Qc "$ReferenceTargetsBedName.part"),
        (Join-Path $Qc $ReferenceCoordinateConflictsName),
        (Join-Path $Qc $ReferenceSourceGenotypeConflictsName)
    )
    return (Get-ArtifactMenuStatus -Ready (Test-ReferenceTargetsReady) -Paths $paths)
}

function Test-ReferenceBamReady {
    if (-not (Test-ReferenceTargetsReady) -or -not (Test-FastaReady)) { return $false }
    $bam = Join-Path $Work $ReferenceBamName
    $statePath = "$bam.state.json"
    foreach ($required in @($bam, "$bam.bai", $statePath)) {
        if (-not (Test-NonEmptyFile $required)) { return $false }
    }
    try {
        $state = Get-Content -LiteralPath $statePath -Raw -Encoding UTF8 | ConvertFrom-Json
        $pair = Get-FastqPairInfo
        if (($state.status -ne "complete") -or ($state.schema -ne "reference-alignment-v3") -or -not (Test-FastqInputReady $pair)) { return $false }
        $readMode = [string]$pair.read_mode
        if ([string]$state.read_mode -ne $readMode) { return $false }
        $r2Ready = $true
        if ($readMode -eq "paired") {
            $r2Ready = (
                (Test-StateFileIdentity $state.inputs.r2 ([string]$pair.r2)) -and
                ([string]$state.inputs.r2.content_sha256 -match '^[0-9a-f]{64}$')
            )
        } elseif ($readMode -eq "single") {
            $r2Ready = ($null -eq $state.inputs.r2)
        } else {
            return $false
        }
        return (
            (Test-StateFileIdentity $state.inputs.r1 ([string]$pair.r1)) -and
            ([string]$state.inputs.r1.content_sha256 -match '^[0-9a-f]{64}$') -and
            $r2Ready -and
            (Test-StateFileIdentity $state.inputs.targets (Join-Path $Qc $ReferenceTargetsBedName)) -and
            (Test-StateFileIdentity $state.inputs.fasta (Join-Path $VepFastaData $FastaName)) -and
            ([string]$state.inputs.fasta.sha256 -match '^[0-9a-f]{64}$') -and
            (Test-ApprovedDockerImage ([string]$state.images.bwa) $BwaImage $BwaFallbackImages) -and
            (Test-ApprovedDockerImage ([string]$state.images.samtools) $SamtoolsImage $SamtoolsFallbackImages) -and
            (Test-StateFileIdentity $state.outputs.bam $bam) -and
            ([string]$state.outputs.bam.sha256 -match '^[0-9a-f]{64}$')
        )
    } catch { return $false }
}

function Get-ReferenceBamMenuStatus {
    $pair = Get-FastqPairInfo
    $pairStatus = Get-FastqPairMenuStatus $pair
    if (-not (Test-FastqInputReady $pair)) { return "BLOCKED: FASTQ $pairStatus" }
    if (-not (Test-FastaReady)) { return "BLOCKED: FASTA MISSING" }
    if (-not (Test-ReferenceTargetsReady)) { return "BLOCKED: TARGETS NOT READY" }
    $paths = @(
        (Join-Path $Work $ReferenceBamName),
        (Join-Path $Work ($ReferenceBamName -replace '\.bam$', '.part.bam'))
    )
    if (Test-ReferenceBamReady) { return "READY" }
    $checkpointDirectory = [Environment]::GetEnvironmentVariable("DANTE_REFERENCE_CHECKPOINT_DIR", "Process")
    if (-not $checkpointDirectory) { $checkpointDirectory = Join-Path $Work "reference_alignment_checkpoints" }
    $checkpointManifest = Join-Path $checkpointDirectory "checkpoint_manifest.json"
    if (Test-NonEmptyFile $checkpointManifest) {
        try {
            $checkpointState = Get-Content -LiteralPath $checkpointManifest -Raw -Encoding UTF8 | ConvertFrom-Json
            if (@($checkpointState.completed_chunks).Count -gt 0) { return "PARTIAL: RESUMABLE CHECKPOINTS" }
        } catch {
            return "PARTIAL: CHECKPOINT STATE"
        }
    }
    $bam = Join-Path $Work $ReferenceBamName
    $statePath = "$bam.state.json"
    if ((Test-NonEmptyFile $bam) -and (Test-NonEmptyFile $statePath)) {
        try {
            $state = Get-Content -LiteralPath $statePath -Raw -Encoding UTF8 | ConvertFrom-Json
            if ([string]$state.read_mode -ne [string]$pair.read_mode) { return "STALE: FASTQ MODE CHANGED" }
        } catch {
            # The generic artifact status below reports malformed or incomplete state.
        }
    }
    return (Get-ArtifactMenuStatus -Ready $false -Paths $paths)
}

function Test-ReferenceCallsReady {
    if (-not (Test-ReferenceBamReady)) { return $false }
    $calls = Join-Path $Work $ReferenceCallsName
    $metadata = Join-Path $Qc $ReferenceCallsMetadataName
    $uncallable = Join-Path $Qc $ReferenceCallsUncallableName
    $statePath = "$calls.state.json"
    foreach ($required in @($calls, "$calls.tbi", $metadata, $uncallable, $statePath)) {
        if (-not (Test-NonEmptyFile $required)) { return $false }
    }
    try {
        $state = Get-Content -LiteralPath $statePath -Raw -Encoding UTF8 | ConvertFrom-Json
        $currentSex = Get-DanteSampleSexMode
        if (($currentSex -ne "AUTO") -and ([string]$state.sample.sex_resolved -ne $currentSex)) { return $false }
        return (
            ($state.status -eq "complete") -and ($state.schema -eq "reference-calls-v5-known-alleles-ploidy-refq") -and
            ($state.calling_strategy -eq "resumable-v4-checkpoints-known-alleles-ploidy-plus-PL-reference-quality-v5") -and
            (Test-StateFileIdentity $state.inputs.bam (Join-Path $Work $ReferenceBamName)) -and
            (Test-StateFileIdentity $state.inputs.targets_bed (Join-Path $Qc $ReferenceTargetsBedName)) -and
            (Test-StateFileIdentity $state.inputs.targets_panel (Join-Path $Qc $ReferenceTargetsName)) -and
            (Test-StateFileIdentity $state.inputs.fasta (Join-Path $VepFastaData $FastaName)) -and
            (Test-StateFileContentIdentity $state.inputs.ploidy_definition (Join-Path $Root "reference_data\GRCh37_bcftools_ploidy.tsv")) -and
            ([string]$state.inputs.bam.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.inputs.targets_bed.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.inputs.targets_panel.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.inputs.fasta.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.inputs.ploidy_definition.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.sample.sex_resolved -match '^[MF]$') -and
            ([double]$state.call_qc.hom_ref_reference_quality_fraction -ge 0.95) -and
            ([double]$state.call_qc.alt_carried_gq_fraction -ge 0.95) -and
            ([int64]$state.call_qc.ploidy_mismatch_records -eq 0) -and
            (Test-ApprovedDockerImage ([string]$state.image) $BcftoolsImage $BcftoolsFallbackImages) -and
            (Test-ApprovedDockerImage ([string]$state.samtools_image) $SamtoolsImage $SamtoolsFallbackImages) -and
            (Test-StateFileIdentity $state.outputs.calls $calls) -and
            ([string]$state.outputs.calls.sha256 -match '^[0-9a-f]{64}$') -and
            (Test-StateFileIdentity $state.outputs.uncallable_reference_loci $uncallable) -and
            ([string]$state.outputs.uncallable_reference_loci.sha256 -match '^[0-9a-f]{64}$')
        )
    } catch { return $false }
}

function Get-ReferenceCallsMenuStatus {
    if (-not (Test-ReferenceBamReady)) { return "BLOCKED: ALIGNMENT NOT READY" }
    if (Test-ReferenceCallsReady) { return "READY" }
    $callsStatePath = "$(Join-Path $Work $ReferenceCallsName).state.json"
    if (Test-NonEmptyFile $callsStatePath) {
        try {
            $callsState = Get-Content -LiteralPath $callsStatePath -Raw -Encoding UTF8 | ConvertFrom-Json
            if ([string]$callsState.schema -in @("reference-calls-v1", "reference-calls-v2-known-alleles-ploidy", "reference-calls-v3-known-alleles-ploidy-nsafe", "reference-calls-v4-known-alleles-ploidy-coordinate-safe")) {
                return "STALE: PRE-REFQ CALL MODEL"
            }
        } catch {
            # Fall through to generic incomplete state handling.
        }
    }
    $checkpointDirectory = [Environment]::GetEnvironmentVariable("DANTE_REFERENCE_CALL_CHECKPOINT_DIR", "Process")
    if (-not $checkpointDirectory) { $checkpointDirectory = Join-Path $Work "reference_call_checkpoints" }
    $checkpointManifest = Join-Path $checkpointDirectory "checkpoint_manifest.json"
    if (Test-NonEmptyFile $checkpointManifest) {
        try {
            $checkpointState = Get-Content -LiteralPath $checkpointManifest -Raw -Encoding UTF8 | ConvertFrom-Json
            if ([string]$checkpointState.schema -ne "reference-call-checkpoints-v4-known-alleles-ploidy-coordinate-safe") {
                return "STALE: PRE-COORDINATE-SAFE CHECKPOINTS"
            }
            if (@($checkpointState.completed_chunks).Count -gt 0) { return "PARTIAL: RESUMABLE CHECKPOINTS" }
            return "PARTIAL: CHECKPOINT STATE"
        } catch {
            return "PARTIAL: CHECKPOINT STATE"
        }
    }
    $paths = @(
        (Join-Path $Work $ReferenceCallsName),
        (Join-Path $Work ($ReferenceCallsName -replace '\.vcf\.gz$', '.part.vcf.gz'))
    )
    return (Get-ArtifactMenuStatus -Ready $false -Paths $paths)
}

function Test-ReferenceSourceMergeReady {
    $export = Join-Path $Results $ReferenceExportName
    $report = Join-Path $Qc "Dante_reference_export_report.txt"
    $discrepancies = Join-Path $Qc "Dante_reference_export_discrepancies.tsv"
    $provenance = Join-Path $Qc $ReferenceProvenanceName
    $statePath = "$export.state.json"
    $base = Join-Path $Results "Dante_23andMe_style.txt.gz"
    if (-not (Test-ReferenceTargetsReady)) { return $false }
    foreach ($required in @($export, $report, $discrepancies, $provenance, $statePath, $base)) {
        if (-not (Test-NonEmptyFile $required)) { return $false }
    }
    try {
        $state = Get-Content -LiteralPath $statePath -Raw -Encoding UTF8 | ConvertFrom-Json
        if ($state.status -ne "complete") { return $false }
        if ([string]$state.mode -eq "fastq_enriched") {
            return (Test-ReferenceExportReady)
        }
        if (($state.schema -ne "reference-export-v3") -or ([string]$state.mode -ne "reference_sources_only")) { return $false }
        return (
            (Test-StateFileIdentity $state.inputs.targets (Join-Path $Qc $ReferenceTargetsName)) -and
            (Test-StateFileIdentity $state.inputs.base $base) -and
            ([string]$state.inputs.targets.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.inputs.base.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.fastq_source -eq "NONE") -and
            (Test-StateFileIdentity $state.outputs.export $export) -and
            (Test-StateFileIdentity $state.outputs.provenance $provenance) -and
            (Test-StateFileIdentity $state.outputs.discrepancies $discrepancies) -and
            ([string]$state.outputs.export.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.outputs.provenance.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.outputs.discrepancies.sha256 -match '^[0-9a-f]{64}$')
        )
    } catch { return $false }
}

function Get-ReferenceSourceMergeMenuStatus {
    if (-not (Test-ReferenceTargetsReady)) { return "BLOCKED: TARGETS NOT READY" }
    $base = Join-Path $Results "Dante_23andMe_style.txt.gz"
    if (-not (Test-NonEmptyFile $base)) { return "BLOCKED: BASE EXPORT MISSING" }
    $paths = @(
        (Join-Path $Results $ReferenceExportName),
        (Join-Path $Results "$ReferenceExportName.part")
    )
    return (Get-ArtifactMenuStatus -Ready (Test-ReferenceSourceMergeReady) -Paths $paths)
}

function Test-ReferenceExportReady {
    if (-not (Test-ReferenceCallsReady)) { return $false }
    $export = Join-Path $Results $ReferenceExportName
    $report = Join-Path $Qc "Dante_reference_export_report.txt"
    $discrepancies = Join-Path $Qc "Dante_reference_export_discrepancies.tsv"
    $provenance = Join-Path $Qc $ReferenceProvenanceName
    $evidence = Join-Path $Results $ReferenceFastqEvidenceName
    $statePath = "$export.state.json"
    $base = Join-Path $Results "Dante_23andMe_style.txt.gz"
    foreach ($required in @($export, $report, $discrepancies, $provenance, $evidence, $statePath, $base)) {
        if (-not (Test-NonEmptyFile $required)) { return $false }
    }
    try {
        $state = Get-Content -LiteralPath $statePath -Raw -Encoding UTF8 | ConvertFrom-Json
        $pair = Get-FastqPairInfo
        $expectedFastqSource = if ([string]$pair.read_mode -eq "single") { "R1" } else { "R1+R2" }
        return (
            ($state.status -eq "complete") -and ($state.schema -eq "reference-export-v4-tiered") -and
            ([string]$state.mode -eq "fastq_enriched") -and
            (Test-StateFileIdentity $state.inputs.calls (Join-Path $Work $ReferenceCallsName)) -and
            (Test-StateFileIdentity $state.inputs.targets (Join-Path $Qc $ReferenceTargetsName)) -and
            (Test-StateFileIdentity $state.inputs.base $base) -and
            ([string]$state.inputs.calls.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.inputs.targets.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.inputs.base.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.thresholds.export_min_tier -eq $ReferenceExportMinTier) -and
            ([int]$state.thresholds.tier_a.min_dp -eq $ReferenceTierAMinDp) -and
            ([int]$state.thresholds.tier_a.min_quality -eq $ReferenceTierAMinQuality) -and
            ([int]$state.thresholds.tier_a.min_allele_depth -eq $ReferenceTierAMinAlleleDepth) -and
            ([double]$state.thresholds.tier_a.hom_fraction -eq [double]$ReferenceTierAHomFraction) -and
            ([double]$state.thresholds.tier_a.het_min_fraction -eq [double]$ReferenceTierAHetMinFraction) -and
            ([double]$state.thresholds.tier_a.het_max_fraction -eq [double]$ReferenceTierAHetMaxFraction) -and
            ([int]$state.thresholds.tier_b.min_dp -eq $ReferenceTierBMinDp) -and
            ([int]$state.thresholds.tier_b.min_quality -eq $ReferenceTierBMinQuality) -and
            ([int]$state.thresholds.tier_b.hom_min_support -eq $ReferenceTierBHomMinSupport) -and
            ([double]$state.thresholds.tier_b.hom_fraction -eq [double]$ReferenceTierBHomFraction) -and
            ([int]$state.thresholds.tier_b.het_min_allele_depth -eq $ReferenceTierBHetMinAlleleDepth) -and
            ([double]$state.thresholds.tier_b.het_min_fraction -eq [double]$ReferenceTierBHetMinFraction) -and
            ([double]$state.thresholds.tier_b.het_max_fraction -eq [double]$ReferenceTierBHetMaxFraction) -and
            ([int]$state.thresholds.tier_c.min_dp -eq $ReferenceTierCMinDp) -and
            ([int]$state.thresholds.tier_c.min_quality -eq $ReferenceTierCMinQuality) -and
            ([string]$state.fastq_source -eq $expectedFastqSource) -and
            (Test-StateFileIdentity $state.outputs.export $export) -and
            (Test-StateFileIdentity $state.outputs.provenance $provenance) -and
            (Test-StateFileIdentity $state.outputs.discrepancies $discrepancies) -and
            (Test-StateFileIdentity $state.outputs.fastq_evidence $evidence) -and
            ([string]$state.outputs.export.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.outputs.provenance.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.outputs.discrepancies.sha256 -match '^[0-9a-f]{64}$') -and
            ([string]$state.outputs.fastq_evidence.sha256 -match '^[0-9a-f]{64}$')
        )
    } catch { return $false }
}

function Get-ReferenceExportMenuStatus {
    if (-not (Test-ReferenceCallsReady)) { return "BLOCKED: CALLS NOT READY" }
    $paths = @(
        (Join-Path $Results $ReferenceExportName),
        (Join-Path $Results "$ReferenceExportName.part")
    )
    return (Get-ArtifactMenuStatus -Ready (Test-ReferenceExportReady) -Paths $paths)
}

function Get-CompleteStageMenuStatus([string[]]$Statuses) {
    foreach ($status in $Statuses) {
        if ($status -ne "READY") { return "" }
    }
    return "READY"
}

function Show-AutomaticCheck {
    Write-Step "Automatic startup check"
    Write-CheckStatus "Active patient" $true "$Patient | $PatientRoot"
    Write-CheckStatus "Patient input" (Test-Path -LiteralPath $PatientInput -PathType Container) $PatientInput
    $dockerCli = [bool](Get-Command docker -ErrorAction SilentlyContinue)
    Write-CheckStatus "Docker CLI" $dockerCli $(if ($dockerCli) { "available" } else { "not found" })

    $dockerEngine = $false
    if ($dockerCli) {
        $saved = $ErrorActionPreference
        try {
            $ErrorActionPreference = "Continue"
            & docker info *> $null
            $dockerEngine = ($LASTEXITCODE -eq 0)
        } finally { $ErrorActionPreference = $saved }
    }
    Write-CheckStatus "Docker engine" $dockerEngine $(if ($dockerEngine) { "running" } else { "not running" })

    $pythonCommand = Get-Command py -ErrorAction SilentlyContinue
    if (-not $pythonCommand) { $pythonCommand = Get-Command python -ErrorAction SilentlyContinue }
    Write-CheckStatus "Python 3" ([bool]$pythonCommand) $(if ($pythonCommand) { $pythonCommand.Source } else { "not found" })
    $curlCommand = Get-Command curl.exe -ErrorAction SilentlyContinue
    Write-CheckStatus "curl.exe" ([bool]$curlCommand) $(if ($curlCommand) { "available" } else { "not found" })

    $inputDefinitions = @(
        @("SNP VCF", '\.filtered\.snp\.vcf(\.gz)?$'),
        @("CNV VCF", '\.cnv\.vcf(\.gz)?$'),
        @("SV VCF", '\.sv\.vcf(\.gz)?$'),
        @("RAW VCF", '\.raw\.vcf(\.gz)?$'),
        @("INDEL VCF", '\.filtered\.indel\.vcf(\.gz)?$')
    )
    foreach ($definition in $inputDefinitions) {
        # Keep the result array-shaped even when exactly one VCF is present.
        # Without @(...), StrictMode treats the single FileInfo as a scalar and
        # accessing .Count terminates the automatic startup check.
        $found = @(Find-DanteCandidates $definition[1])
        $ok = ($found.Count -eq 1)
        $detail = if ($found.Count -eq 1) { $found[0].Name } elseif ($found.Count -eq 0) { "missing" } else { "$($found.Count) matching files" }
        Write-CheckStatus $definition[0] $ok $detail
    }
    $fastq = Get-FastqPairInfo
    $fastqDetail = if (@($fastq.details).Count -gt 0) { (@($fastq.details) -join "; ") } else { [string]$fastq.status }
    Write-CheckStatus "FASTQ input" (Test-FastqInputReady $fastq) $fastqDetail
    Write-CheckStatus ".env file" (Test-Path -LiteralPath $EnvFile -PathType Leaf) $(if (Test-Path -LiteralPath $EnvFile) { "loaded" } else { "missing" })
    Write-CheckStatus "OpenAI API key" (Test-ApiKeyConfigured) $(if (Test-ApiKeyConfigured) { "configured (value hidden)" } else { "placeholder or missing; only stage 5 is affected" })
}

function Install-VepCacheData {
    Write-Step "Install/resume Ensembl VEP 116 GRCh37 cache"
    Assert-Docker
    New-Item -ItemType Directory -Force -Path $VepData | Out-Null
    $cacheDir = Join-Path $VepData "homo_sapiens\116_GRCh37"
    $cacheArchive = Join-Path $VepData $VepCacheArchiveName
    $readyMarker = Join-Path $VepData ".dante_vep_116_grch37.complete"
    if ((Test-VepCacheReady) -and -not $ForceDownload) {
        if (-not (Test-Path -LiteralPath $readyMarker -PathType Leaf)) {
            Set-AtomicText $readyMarker "complete|VEP=116|assembly=GRCh37"
        }
        Write-Host "VEP cache already exists: $cacheDir"
        return
    }
    [void](Resolve-DockerImage $VepImage @() "VEP 116")
    if ($ForceDownload) {
        if (Test-Path -LiteralPath $cacheArchive) { Remove-Item -LiteralPath $cacheArchive -Force }
        Remove-VepCacheTree $cacheDir
        if (Test-Path -LiteralPath $readyMarker) { Remove-Item -LiteralPath $readyMarker -Force }
    }
    $remoteLength = Get-VepCacheRemoteLength
    $localLength = if (Test-Path -LiteralPath $cacheArchive) { (Get-Item -LiteralPath $cacheArchive).Length } else { [long]0 }
    $archiveIntegrityValidated = $false
    $downloadRequired = $true
    if ($localLength -gt 0 -and $remoteLength -gt 0) {
        if ($localLength -eq $remoteLength) {
            $downloadRequired = $false
            Write-Host "The archive has the complete server size: $(Format-ByteSize $localLength)."
        } elseif ($localLength -gt $remoteLength) {
            throw "The local VEP archive is larger than the server file ($localLength > $remoteLength bytes). It was preserved. Re-run item 1 with -ForceDownload only if you want to discard it and download again."
        }
    } elseif ($localLength -gt 0 -and $remoteLength -le 0) {
        Write-Host "The server size is temporarily unavailable. Checking the existing archive before any download..."
        $archiveIntegrityValidated = Test-GzipInVepData $VepCacheArchiveName
        $downloadRequired = -not $archiveIntegrityValidated
    }
    if ($downloadRequired) {
        Write-Host "Downloading the official indexed cache. Network resets will resume automatically."
        Invoke-ResumableDownload -Url $VepCacheUrl -Destination $cacheArchive -RemoteLength $remoteLength
        $localLength = if (Test-Path -LiteralPath $cacheArchive) { [long](Get-Item -LiteralPath $cacheArchive).Length } else { [long]0 }
        if ($remoteLength -gt 0 -and $localLength -ne $remoteLength) {
            throw "The download command ended, but the archive size is $localLength bytes instead of the expected $remoteLength bytes. The downloaded bytes were preserved. Select item 1 to resume."
        }
    } else {
        Write-Host "A complete cache archive was found; reusing it."
    }
    if (-not $archiveIntegrityValidated) {
        Write-Host "Verifying gzip integrity once. This reads the complete 23-GiB archive and can take several minutes..."
        $archiveIntegrityValidated = Test-GzipInVepData $VepCacheArchiveName
    }
    if (-not $archiveIntegrityValidated) {
        throw "The VEP cache archive is incomplete or corrupt. Its partial bytes are preserved for another resume attempt."
    }
    if (Test-Path -LiteralPath $readyMarker) { Remove-Item -LiteralPath $readyMarker -Force }
    if (Test-Path -LiteralPath $cacheDir) {
        Write-Host "Removing a previous incomplete extraction before restarting tar..." -ForegroundColor Yellow
        Remove-VepCacheTree $cacheDir
    }
    Write-Host "Archive integrity is valid. Extracting the GRCh37 cache..."
    Invoke-ObservedNative -Activity "Extract VEP 116 GRCh37 cache" -Program "docker" -Arguments @(
        "run", "--rm", "--mount", "type=bind,source=$VepData,target=/data",
        $VepImage, "bash", "-lc", "tar -xzf /data/$VepCacheArchiveName -C /data"
    ) -WatchDirectory $cacheDir
    if (-not ((Test-Path -LiteralPath (Join-Path $cacheDir "info.txt") -PathType Leaf) -and
        (Test-Path -LiteralPath (Join-Path $cacheDir "1") -PathType Container))) {
        throw "The expected VEP 116 GRCh37 cache structure was not created. The source archive was preserved."
    }
    Set-AtomicText $readyMarker "complete|VEP=116|assembly=GRCh37"
    Remove-Item -LiteralPath $cacheArchive -Force
    Write-Host "The VEP 116 GRCh37 cache is ready."
}

function Ensure-VepCacheReady {
    if (Test-VepCacheReady) { return }
    $cacheArchive = Join-Path $VepData $VepCacheArchiveName
    if (Test-Path -LiteralPath $cacheArchive -PathType Leaf) {
        Write-Host "The VEP archive exists but is not installed yet. Verifying and extracting it before continuing..." -ForegroundColor Yellow
        Install-VepCacheData
        if (Test-VepCacheReady) { return }
    }
    throw "VEP GRCh37 cache is not installed. Select item 1 first."
}

function Test-AnnotSVFunctional {
    $probeInput = Join-Path $AnnotSVData ".dante_annotsv_smoke_input.bed"
    $probeDirectory = Join-Path $AnnotSVData ".dante_annotsv_smoke"
    $probeOutput = Join-Path $probeDirectory "smoke.annotated.tsv"
    if (Test-Path -LiteralPath $probeDirectory) {
        Remove-Item -LiteralPath $probeDirectory -Recurse -Force
    }
    New-Item -ItemType Directory -Force -Path $probeDirectory | Out-Null
    [System.IO.File]::WriteAllText(
        $probeInput,
        "1`t100000`t200000`tDEL`n",
        [System.Text.UTF8Encoding]::new($false)
    )
    try {
        Invoke-ObservedNative -Activity "Validate AnnotSV GRCh37 annotations" -Program "docker" -Arguments @(
            "run", "--rm", "--mount", "type=bind,source=$AnnotSVData,target=/data",
            $AnnotSVImage, "AnnotSV",
            "-SVinputFile", "/data/.dante_annotsv_smoke_input.bed",
            "-annotationsDir", "/data/AnnotSV_annotations",
            "-genomeBuild", "GRCh37", "-outputDir", "/data/.dante_annotsv_smoke",
            "-outputFile", "smoke.annotated.tsv", "-overwrite", "1", "-svtBEDcol", "4"
        ) -WatchDirectory $probeDirectory
        return (
            (Test-Path -LiteralPath $probeOutput -PathType Leaf) -and
            ((Get-Item -LiteralPath $probeOutput).Length -gt 0)
        )
    } catch {
        Write-Host ("AnnotSV validation failed: {0}" -f $_.Exception.Message) -ForegroundColor Yellow
        return $false
    } finally {
        if (Test-Path -LiteralPath $probeInput) { Remove-Item -LiteralPath $probeInput -Force }
        if (Test-Path -LiteralPath $probeDirectory) {
            Remove-Item -LiteralPath $probeDirectory -Recurse -Force
        }
    }
}

function Test-AnnotSVHumanArchive([string]$ArchivePath) {
    if (-not (Test-Path -LiteralPath $ArchivePath -PathType Leaf)) { return $false }
    try {
        Invoke-ObservedNative -Activity "Verify AnnotSV human annotation archive" -Program "docker" -Arguments @(
            "run", "--rm", "--mount", "type=bind,source=$AnnotSVData,target=/data",
            $AnnotSVImage, "bash", "-lc",
            "tar -tzf /data/AnnotSV_annotations/$AnnotSVHumanArchiveName >/dev/null"
        ) -WatchFile $ArchivePath
        return $true
    } catch {
        return $false
    }
}

function Install-AnnotSVData {
    Write-Step "Install AnnotSV 3.5.10 public annotations"
    Assert-Docker
    New-Item -ItemType Directory -Force -Path $AnnotSVData | Out-Null
    $readyMarker = Get-AnnotSVCompletionMarker
    if ((Test-AnnotSVReady) -and -not $ForceDownload) {
        Write-Host "AnnotSV annotations already passed the GRCh37 functional validation."
        return
    }
    if (Test-Path -LiteralPath $readyMarker) {
        Remove-Item -LiteralPath $readyMarker -Force
    }
    [void](Resolve-DockerImage $AnnotSVImage @() "AnnotSV 3.5.10")
    $annotationDirectory = Join-Path $AnnotSVData "AnnotSV_annotations"
    $humanArchive = Join-Path $annotationDirectory $AnnotSVHumanArchiveName
    if ((Test-AnnotSVHumanStructureReady) -and -not $ForceDownload) {
        Write-Host "An existing AnnotSV GRCh37 data tree was found. Validating it without re-extraction..." -ForegroundColor Yellow
    } else {
        if (Test-Path -LiteralPath $annotationDirectory -PathType Container) {
            Write-Host "An incomplete AnnotSV installation was found. Repairing it in place..." -ForegroundColor Yellow
        }
        New-Item -ItemType Directory -Force -Path $annotationDirectory | Out-Null
        if ($ForceDownload -and (Test-Path -LiteralPath $humanArchive)) {
            Remove-Item -LiteralPath $humanArchive -Force
        }
        $remoteLength = Get-RemoteContentLength $AnnotSVHumanArchiveUrl
        $archiveReady = $false
        if ((Test-Path -LiteralPath $humanArchive -PathType Leaf) -and $remoteLength -gt 0 -and
            ((Get-Item -LiteralPath $humanArchive).Length -gt $remoteLength)) {
            Write-Host "The local AnnotSV archive is larger than the official file and will be replaced." -ForegroundColor Yellow
            Remove-Item -LiteralPath $humanArchive -Force
        }
        $downloadCompleteBySize = (
            (Test-Path -LiteralPath $humanArchive -PathType Leaf) -and
            ($remoteLength -gt 0) -and
            ((Get-Item -LiteralPath $humanArchive).Length -eq $remoteLength)
        )
        if ($downloadCompleteBySize) {
            Write-Host "The AnnotSV human archive has the complete official size; validating it before reuse."
            $archiveReady = Test-AnnotSVHumanArchive $humanArchive
            if (-not $archiveReady) {
                Write-Host "The size-matched archive is corrupt and will be downloaded again." -ForegroundColor Yellow
                Remove-Item -LiteralPath $humanArchive -Force
            }
        } elseif ((Test-Path -LiteralPath $humanArchive -PathType Leaf) -and $remoteLength -le 0) {
            Write-Host "The official size is temporarily unavailable; checking the existing archive before any network retry."
            $archiveReady = Test-AnnotSVHumanArchive $humanArchive
        }
        if (-not $archiveReady) {
            Invoke-ResumableDownload $AnnotSVHumanArchiveUrl $humanArchive $remoteLength
        }
        if (-not $archiveReady -and -not (Test-AnnotSVHumanArchive $humanArchive)) {
            Remove-Item -LiteralPath $humanArchive -Force
            throw "The AnnotSV human archive failed tar/gzip validation and was removed. Rerun item 3 or 4 to download a clean copy."
        }
        Invoke-ObservedNative -Activity "Extract AnnotSV human annotations" -Program "docker" -Arguments @(
            "run", "--rm", "--mount", "type=bind,source=$AnnotSVData,target=/data",
            $AnnotSVImage, "bash", "-lc",
            "tar -xzf /data/AnnotSV_annotations/$AnnotSVHumanArchiveName -C /data/AnnotSV_annotations"
        ) -WatchDirectory $annotationDirectory
    }
    if (-not (Test-AnnotSVHumanStructureReady)) {
        throw "AnnotSV installation is incomplete: neither a non-empty genes.RefSeq.sorted.bed nor ncbiRefSeq.txt.gz was found for GRCh37. The partial directory was preserved; rerun item 3 or 4 to repair it."
    }
    Write-Host "Running a one-interval AnnotSV smoke test before publishing READY..."
    if (-not (Test-AnnotSVFunctional)) {
        throw "AnnotSV GRCh37 functional validation failed. No completion marker was written; rerun item 3 or 4 after checking the preceding AnnotSV error."
    }
    Set-AtomicText $readyMarker "complete|AnnotSV=3.5.10|human_annotations=validated|assembly=GRCh37|smoke=passed"
    if (Test-Path -LiteralPath $humanArchive) {
        Remove-Item -LiteralPath $humanArchive -Force
    }
    Write-Host "AnnotSV public annotations are complete and functionally validated."
}

function Ensure-AnnotSVReady {
    if (Test-AnnotSVReady) { return }
    Write-Host "AnnotSV annotations are missing or incomplete. Installing/repairing them before continuing..." -ForegroundColor Yellow
    Install-AnnotSVData
    if (-not (Test-AnnotSVReady)) {
        throw "AnnotSV GRCh37 annotations are still incomplete after repair."
    }
}

function Run-RsidExport {
    $cacheDir = Join-Path $VepData "homo_sapiens\116_GRCh37"
    Ensure-VepCacheReady
    Assert-Docker
    New-Item -ItemType Directory -Force -Path $Results, $Work, $Qc | Out-Null
    $inputs = Get-DanteInputs @("snp")
    $python = Get-Python
    $vepOutputName = "Dante_SNP_rsid_lookup.vep.vcf.gz"
    $vepOutput = Join-Path $Work $vepOutputName
    $vepPartName = "Dante_SNP_rsid_lookup.part.vep.vcf.gz"
    $vepPart = Join-Path $Work $vepPartName
    $vepState = "$vepOutput.state.sha256"
    $export = Join-Path $Results "Dante_23andMe_style.txt.gz"
    $exportPart = "$export.part"
    $report = Join-Path $Qc "Dante_23andMe_style_report.txt"
    $exportState = "$export.state.sha256"
    $cacheInfo = Join-Path $cacheDir "info.txt"

    Write-Host "Calculating the source identity for safe restart/reuse..."
    # v7.4.1 is a parser-only hotfix. Keep the v7.4 fingerprint so completed
    # v7.4-compatible intermediates remain reusable after the update.
    $fingerprint = "rsid-export-v7.4|VEP=116|GRCh37|$((Get-Sha256 $cacheInfo))|$((Get-Sha256 $inputs.snp.FullName))"

    if ((Test-Path -LiteralPath $export -PathType Leaf) -and
        (Test-Path -LiteralPath $report -PathType Leaf) -and
        (Test-Path -LiteralPath $exportState -PathType Leaf) -and
        ((Get-Content -LiteralPath $exportState -Raw).Trim() -eq $fingerprint) -and
        -not $ForceDownload) {
        Write-Host "The 23andMe-style export has matching state. Verifying gzip before reuse..."
        if (Test-GzipFile $export) {
            Write-Host "The 23andMe-style export is current for the same SNP VCF and VEP cache; reusing it."
            return
        }
        Write-Host "The state-matched export is not a complete gzip and will be rebuilt." -ForegroundColor Yellow
    }

    if ((Test-Path -LiteralPath $export -PathType Leaf) -and
        (Test-Path -LiteralPath $report -PathType Leaf) -and
        -not (Test-Path -LiteralPath $exportState -PathType Leaf) -and
        ((Get-Item -LiteralPath $export).LastWriteTimeUtc -ge $inputs.snp.LastWriteTimeUtc) -and
        -not $ForceDownload) {
        Write-Host "A completed v7.3 export was found. Verifying it once before adoption..." -ForegroundColor Yellow
        if (Test-GzipFile $export) {
            Set-AtomicText $exportState $fingerprint
            Write-Host "The completed v7.3 export was adopted; item 2 does not need to be recalculated."
            return
        }
        Write-Host "The legacy export is not a complete gzip and will be rebuilt." -ForegroundColor Yellow
    }

    if (Test-Path -LiteralPath $exportState) {
        Remove-Item -LiteralPath $exportState -Force
    }

    $vepReady = $false
    if ((Test-Path -LiteralPath $vepOutput -PathType Leaf) -and
        (Test-Path -LiteralPath $vepState -PathType Leaf) -and
        ((Get-Content -LiteralPath $vepState -Raw).Trim() -eq $fingerprint) -and
        -not $ForceDownload) {
        Write-Host "Checking the completed rsID VEP intermediate before reuse..."
        $vepReady = Test-GzipFile $vepOutput
    }
    if (-not $vepReady -and (Test-Path -LiteralPath $vepOutput -PathType Leaf) -and
        -not (Test-Path -LiteralPath $vepState -PathType Leaf) -and
        ((Get-Item -LiteralPath $vepOutput).LastWriteTimeUtc -ge $inputs.snp.LastWriteTimeUtc) -and
        -not $ForceDownload) {
        Write-Host "A v7.3 rsID VEP intermediate was found. Verifying it once..." -ForegroundColor Yellow
        if (Test-GzipFile $vepOutput) {
            Set-AtomicText $vepState $fingerprint
            $vepReady = $true
            Write-Host "The complete v7.3 VEP intermediate was adopted."
        } else {
            Write-Host "The v7.3 VEP intermediate is incomplete and cannot be resumed; it will be recalculated." -ForegroundColor Yellow
        }
    }

    if (-not $vepReady) {
        if (Test-Path -LiteralPath $vepPart) {
            Write-Host "Removing an interrupted VEP .part file; VEP itself cannot resume mid-file." -ForegroundColor Yellow
            Remove-Item -LiteralPath $vepPart -Force
        }
        if (Test-Path -LiteralPath $vepOutput) { Remove-Item -LiteralPath $vepOutput -Force }
        if (Test-Path -LiteralPath $vepState) { Remove-Item -LiteralPath $vepState -Force }
        Write-Step "Allele-specific rsID matching for Dante SNPs"
        Invoke-ObservedNative -Activity "Minimal VEP rsID lookup" -Program "docker" -Arguments @(
            "run", "--rm",
            "--mount", "type=bind,source=$Root,target=/work",
            "--mount", "type=bind,source=$VepData,target=/data",
            $VepImage, "vep", "--cache", "--offline", "--cache_version", $VepCacheVersion,
            "--species", "homo_sapiens", "--assembly", "GRCh37", "--dir_cache", "/data",
            "--format", "vcf", "--vcf", "--check_existing", "--exclude_null_alleles",
            "--allele_number", "--fields", "Allele,Consequence,Existing_variation,ALLELE_NUM",
            "--fork", "$Threads", "--buffer_size", "10000", "--no_stats",
            "--compress_output", "bgzip", "--force_overwrite",
            "--input_file", "/work/$($inputs.snp.Name)",
            "--output_file", "/work/Dante_work/$vepPartName",
            "--warning_file", "/work/Dante_qc/Dante_SNP_rsid_lookup.warnings.txt"
        ) -WatchFile $vepPart
        if (-not (Test-GzipFile $vepPart)) {
            throw "The rsID VEP command ended, but its .part output is not a complete gzip. It was preserved for diagnosis and will be replaced on the next run."
        }
        Move-Item -LiteralPath $vepPart -Destination $vepOutput -Force
        Set-AtomicText $vepState $fingerprint
        Write-Host "Minimal VEP intermediate is complete and restart-safe."
    } else {
        Write-Host "Minimal VEP intermediate is current; continuing with the four-column export."
    }

    Write-Step "Create strict four-column 23andMe-style export"
    if (Test-Path -LiteralPath $exportPart) { Remove-Item -LiteralPath $exportPart -Force }
    Invoke-ObservedNative -Activity "Build four-column 23andMe-style export" -Program $python -Arguments @(
        (Join-Path $Root "vep_to_23andme.py"), $vepOutput,
        "--output", $export,
        "--report", $report
    ) -WatchFile $exportPart
    if (-not (Test-Path -LiteralPath $export)) { throw "The 23andMe-style output was not created." }
    if (-not (Test-GzipFile $export)) { throw "The 23andMe-style output failed its gzip integrity check." }
    Set-AtomicText $exportState $fingerprint
    if (-not $KeepIntermediate) {
        Remove-Item -LiteralPath $vepOutput -Force
        if (Test-Path -LiteralPath $vepState) { Remove-Item -LiteralPath $vepState -Force }
    }
    Write-Host "Created: $export"
}

function Invoke-AnnotSVOne([string]$Python, [System.IO.FileInfo]$InputFile, [string]$Kind) {
    $upper = $Kind.ToUpper()
    $preparedName = "Dante_$($upper)_PASS_canonical.vcf"
    $prepared = Join-Path $Work $preparedName
    $report = Join-Path $Qc "Dante_$($upper)_AnnotSV_input_report.txt"
    $annotatedName = "Dante_$($upper)_AnnotSV.tsv"
    $annotated = Join-Path $Results $annotatedName
    if ((Test-AnnotSVOutputCurrent $annotated $report $InputFile) -and -not $ForceDownload) {
        Write-Host "[READY] $annotatedName passed the AnnotSV schema/data check and is newer than its source VCF and annotation marker; reusing it." -ForegroundColor Green
        return
    }
    Invoke-ObservedNative -Activity "Prepare $upper input for AnnotSV" -Program $Python -Arguments @(
        (Join-Path $Root "prepare_annotsv_input.py"), $InputFile.FullName, $prepared,
        "--kind", $Kind,
        "--report", $report
    ) -WatchFile "$prepared.part"
    $annotatedPartName = "Dante_$($upper)_AnnotSV.part.tsv"
    $annotatedPart = Join-Path $Results $annotatedPartName
    $annotsvContainerName = "dante-annotsv-$Kind-$PID"
    if (Test-Path -LiteralPath $annotatedPart) { Remove-Item -LiteralPath $annotatedPart -Force }
    Write-Step "AnnotSV annotation: $upper"
    Invoke-ObservedNative -Activity "AnnotSV annotation: $upper" -Program "docker" -Arguments @(
        "run", "--rm", "--name", $annotsvContainerName,
        "--mount", "type=bind,source=$Root,target=/work",
        $AnnotSVImage, "AnnotSV",
        "-SVinputFile", "/work/Dante_work/$preparedName",
        "-annotationsDir", "/work/annotation_data/annotsv/AnnotSV_annotations",
        "-genomeBuild", "GRCh37", "-outputDir", "/work/Dante_results",
        "-outputFile", $annotatedPartName, "-overwrite", "1", "-includeCI", "1",
        "-SVinputInfo", "1", "-annotationMode", "both", "-tx", "RefSeq", "-vcf", "0"
    ) -WatchFile $annotatedPart -ProgressProfile "annotsv" -ContainerName $annotsvContainerName
    if (-not (Test-Path -LiteralPath $annotatedPart -PathType Leaf)) { throw "AnnotSV did not create $annotatedPartName" }
    if (-not (Test-AnnotSVOutputReady $annotatedPart $report)) {
        throw "AnnotSV ended, but $annotatedPartName failed the schema/data completeness check. The partial file was preserved and was not published as READY."
    }
    Move-Item -LiteralPath $annotatedPart -Destination $annotated -Force
    if (-not $KeepIntermediate) { Remove-Item -LiteralPath $prepared -Force }
    Write-Host "Created: $annotated"
}

function Run-StructuralAnnotation {
    Ensure-AnnotSVReady
    Assert-Docker
    New-Item -ItemType Directory -Force -Path $Results, $Work, $Qc | Out-Null
    $inputs = Get-DanteInputs @("cnv", "sv")
    $python = Get-Python
    Invoke-AnnotSVOne $python $inputs.cnv "cnv"
    Invoke-AnnotSVOne $python $inputs.sv "sv"
}

function Install-FastaData {
    Write-Step "Install GRCh37 primary-assembly reference FASTA"
    Assert-Docker
    New-Item -ItemType Directory -Force -Path $VepFastaData | Out-Null
    $archive = Join-Path $VepFastaData $FastaArchiveName
    $fasta = Join-Path $VepFastaData $FastaName
    if ((Test-FastaReady) -and -not $ForceDownload) {
        Write-Host "Reference FASTA already exists: $fasta"
        return
    }
    [void](Resolve-DockerImage $VepImage @() "VEP 116")
    if ($ForceDownload) {
        if (Test-Path -LiteralPath $archive) { Remove-Item -LiteralPath $archive -Force }
        if (Test-Path -LiteralPath $fasta) { Remove-Item -LiteralPath $fasta -Force }
    }
    $remoteLength = Get-RemoteContentLength $FastaUrl
    $localLength = if (Test-Path -LiteralPath $archive) { (Get-Item -LiteralPath $archive).Length } else { [long]0 }
    $archiveComplete = if ($remoteLength -gt 0) {
        ($localLength -eq $remoteLength -and (Test-GzipInVepData "fasta/$FastaArchiveName"))
    } else {
        (Test-GzipInVepData "fasta/$FastaArchiveName")
    }
    if (-not $archiveComplete) {
        if ($remoteLength -gt 0 -and $localLength -gt $remoteLength) {
            Remove-Item -LiteralPath $archive -Force
        }
        Invoke-ResumableDownload -Url $FastaUrl -Destination $archive -RemoteLength $remoteLength
    }
    if (-not (Test-GzipInVepData "fasta/$FastaArchiveName")) {
        throw "The FASTA archive is incomplete or corrupt. Its partial bytes are preserved."
    }
    $fastaPart = "$fasta.part"
    if (Test-Path -LiteralPath $fastaPart) { Remove-Item -LiteralPath $fastaPart -Force }
    Invoke-ObservedNative -Activity "Extract GRCh37 reference FASTA" -Program "docker" -Arguments @(
        "run", "--rm", "--mount", "type=bind,source=$VepData,target=/data",
        $VepImage, "bash", "-lc",
        "gzip -dc /data/fasta/$FastaArchiveName > /data/fasta/$FastaName.part && mv /data/fasta/$FastaName.part /data/fasta/$FastaName"
    ) -WatchFile $fastaPart
    if (-not (Test-FastaReady)) { throw "The FASTA was extracted, but its expected size/structure was not found." }
    Remove-Item -LiteralPath $archive -Force
    Write-Host "GRCh37 FASTA is ready: $fasta"
}

function Install-ClinVarData {
    Write-Step "Install/update current ClinVar GRCh37"
    New-Item -ItemType Directory -Force -Path $ClinVarData | Out-Null
    $curl = Get-Curl
    $vcf = Join-Path $ClinVarData "clinvar.vcf.gz"
    $index = "$vcf.tbi"
    $md5 = "$vcf.md5"
    $md5Part = "$md5.part"
    Invoke-NativeWithRetry $curl @("--fail", "--location", "--output", $md5Part, $ClinVarMd5Url)
    Move-Item -LiteralPath $md5Part -Destination $md5 -Force
    $md5Text = Get-Content -LiteralPath $md5 -Raw
    $match = [regex]::Match($md5Text, '(?i)\b[a-f0-9]{32}\b')
    if (-not $match.Success) { throw "Could not parse the official ClinVar MD5 file." }
    $expectedMd5 = $match.Value.ToLowerInvariant()
    $current = $false
    if ((Test-Path -LiteralPath $vcf -PathType Leaf) -and -not $ForceDownload) {
        $actualMd5 = Get-FileHashTracked -Path $vcf -Algorithm "MD5" -Label "ClinVar VCF"
        $current = ($actualMd5 -eq $expectedMd5)
    }
    if (-not $current) {
        if (Test-Path -LiteralPath $vcf) {
            if ($ForceDownload -or (Test-GzipWithPython $vcf)) {
                # A complete old weekly ClinVar file cannot be range-resumed
                # into the current release. An incomplete gzip is preserved.
                Remove-Item -LiteralPath $vcf -Force
            }
        }
        Invoke-ResumableDownload $ClinVarUrl $vcf
        $actualMd5 = Get-FileHashTracked -Path $vcf -Algorithm "MD5" -Label "ClinVar VCF"
        if ($actualMd5 -ne $expectedMd5) {
            Remove-Item -LiteralPath $vcf -Force
            throw "ClinVar download failed its official MD5 check. The invalid hybrid file was removed; rerun item 7."
        }
    } else {
        Write-Host "The installed ClinVar VCF already matches the current official MD5."
    }
    $indexPart = "$index.part"
    Invoke-NativeWithRetry $curl @("--fail", "--location", "--output", $indexPart, $ClinVarIndexUrl)
    Move-Item -LiteralPath $indexPart -Destination $index -Force
    if (-not (Test-ClinVarReady)) { throw "ClinVar VCF or its tabix index is missing after installation." }
    Write-Host "ClinVar GRCh37 is ready: $vcf"
}

function Invoke-BcftoolsNormalize(
    [System.IO.FileInfo]$InputFile,
    [string]$OutputName,
    [string]$Label,
    [string]$FastaHash,
    [string]$ReferenceTargetsHash,
    [string]$BcftoolsRuntimeImage
) {
    $output = Join-Path $Work $OutputName
    $state = "$output.state.sha256"
    $fingerprint = "$BcftoolsToolContract|norm-m-any|keep-sum-AD|reference-target-filter-v1|$FastaHash|$ReferenceTargetsHash|$(Get-Sha256 $InputFile.FullName)"
    if ((Test-Path -LiteralPath $output -PathType Leaf) -and (Test-Path -LiteralPath $state -PathType Leaf) -and
        ((Get-Content -LiteralPath $state -Raw).Trim() -eq $fingerprint) -and -not $ForceDownload) {
        Write-Host "$Label normalized layer is current; reusing $OutputName"
        return
    }
    $partName = ($OutputName -replace '\.vcf\.gz$', '.part.vcf.gz')
    $partPath = Join-Path $Work $partName
    if (Test-Path -LiteralPath $partPath) { Remove-Item -LiteralPath $partPath -Force }
    if (Test-Path -LiteralPath "$partPath.tbi") { Remove-Item -LiteralPath "$partPath.tbi" -Force }
    Write-Step "Normalize and split multiallelic $Label variants"
    Invoke-ObservedNative -Activity "bcftools normalize $Label" -Program "docker" -Arguments @(
        "run", "--rm",
        "--mount", "type=bind,source=$Root,target=/work",
        "--mount", "type=bind,source=$VepData,target=/data",
        "--entrypoint", "bcftools", $BcftoolsRuntimeImage, "norm", "--check-ref", "e",
        "--fasta-ref", "/data/fasta/$FastaName", "--multiallelics", "-any",
        "--keep-sum", "AD", "--targets-file", "/work/Dante_qc/$ReferenceContigTargetsName",
        "--output-type", "z",
        "--output", "/work/Dante_work/$partName", "/work/$($InputFile.Name)"
    ) -WatchFile $partPath
    Invoke-ObservedNative -Activity "Index normalized $Label layer" -Program "docker" -Arguments @(
        "run", "--rm", "--mount", "type=bind,source=$Root,target=/work",
        "--entrypoint", "bcftools", $BcftoolsRuntimeImage, "index", "--tbi", "--force", "/work/Dante_work/$partName"
    ) -WatchFile "$partPath.tbi"
    Move-Item -LiteralPath $partPath -Destination $output -Force
    Move-Item -LiteralPath "$partPath.tbi" -Destination "$output.tbi" -Force
    Set-AtomicText $state $fingerprint
}

function Assert-BcftoolsCompatibility([string]$BcftoolsRuntimeImage) {
    Write-Host "Checking approved bcftools 1.24 image..."

    $versionProbe = Invoke-NativeCapture "docker" @(
        "run", "--rm", "--entrypoint", "bcftools", $BcftoolsRuntimeImage, "--version"
    )
    if (($versionProbe.ExitCode -ne 0) -or ($versionProbe.Output -notmatch '(?m)^bcftools 1\.24(?:\s|$)')) {
        throw ("The selected approved bcftools image does not report exact bcftools 1.24. " +
            "Image: $BcftoolsRuntimeImage. Exit: $($versionProbe.ExitCode). Output: " +
            (Get-CommandProbeTail $versionProbe.Output))
    }
    Write-Host "Detected bcftools 1.24" -ForegroundColor Green

    # Stage 2 uses a bash streaming/helper path. Probe the actual executable
    # directly instead of parsing shell help text through another shell.
    $bashProbe = Invoke-NativeCapture "docker" @(
        "run", "--rm", "--entrypoint", "bash", $BcftoolsRuntimeImage, "--version"
    )
    if (($bashProbe.ExitCode -ne 0) -or ($bashProbe.Output -notmatch 'GNU bash')) {
        throw ("The selected approved bcftools image lacks bash required by the current runtime. " +
            "Image: $BcftoolsRuntimeImage. Output: " + (Get-CommandProbeTail $bashProbe.Output))
    }
    Write-Host "bcftools 1.24 runtime check passed." -ForegroundColor Green
}

function Ensure-ReferenceContigTargets([string]$FastaHash, [string]$BcftoolsRuntimeImage) {
    $targets = Join-Path $Qc $ReferenceContigTargetsName
    $state = "$targets.state.sha256"
    $fingerprint = "reference-contig-targets-v1|$FastaHash"
    if ((Test-NonEmptyFile $targets) -and (Test-NonEmptyFile $state) -and
        ((Get-Content -LiteralPath $state -Raw).Trim() -eq $fingerprint) -and -not $ForceDownload) {
        Write-Host "Reference-contig target list is current; reusing $ReferenceContigTargetsName"
        return $targets
    }

    $part = "$targets.part"
    if (Test-Path -LiteralPath $part) { Remove-Item -LiteralPath $part -Force }
    Write-Step "Build normalization target list from installed GRCh37 FASTA"
    Invoke-ObservedNative -Activity "Read GRCh37 FASTA contigs" -Program "docker" -Arguments @(
        "run", "--rm",
        "--mount", "type=bind,source=$Root,target=/work",
        "--mount", "type=bind,source=$VepData,target=/data",
        "--entrypoint", "sh", $BcftoolsRuntimeImage, "-c",
        "if [ -s /data/fasta/$FastaName.fai ] && [ /data/fasta/$FastaName.fai -nt /data/fasta/$FastaName ]; then sed 's/^\([^\t]*\)\t\([0-9]*\).*/\1\t1\t\2/' /data/fasta/$FastaName.fai; else sed -n 's/^>\([^[:space:]]*\).*/\1\t1\t2147483647/p' /data/fasta/$FastaName; fi > /work/Dante_qc/$ReferenceContigTargetsName.part"
    ) -WatchFile $part
    if (-not (Test-NonEmptyFile $part)) {
        throw "No reference contigs were extracted from the installed GRCh37 FASTA."
    }
    $targetCount = @(Get-Content -LiteralPath $part -Encoding UTF8).Count
    if ($targetCount -lt 24) {
        throw "Only $targetCount reference contigs were extracted from the GRCh37 FASTA; refusing to normalize."
    }
    Move-Item -LiteralPath $part -Destination $targets -Force
    Set-AtomicText $state $fingerprint
    Write-Host "Reference-contig target list is ready: $targetCount FASTA sequences" -ForegroundColor Green
    return $targets
}

function Prepare-NormalizedLayers {
    if (-not (Test-FastaReady)) { throw "GRCh37 FASTA is not installed. Select item 6 first." }
    Assert-Docker
    $python = Get-Python
    $inputs = Get-DanteInputs @("raw", "snp", "indel")
    New-Item -ItemType Directory -Force -Path $Work, $Qc | Out-Null
    $bcftoolsRuntimeImage = Get-BcftoolsRuntimeImage
    Assert-BcftoolsCompatibility $bcftoolsRuntimeImage
    $fasta = Join-Path $VepFastaData $FastaName
    Write-Host "Hashing the reference FASTA once for reproducible step state..."
    $fastaHash = Get-Sha256 $fasta
    $referenceTargets = Ensure-ReferenceContigTargets $fastaHash $bcftoolsRuntimeImage
    $referenceTargetsHash = Get-Sha256 $referenceTargets
    Invoke-Native $python @(
        (Join-Path $Root "reference_contig_audit.py"),
        "--targets", $referenceTargets,
        "--report", (Join-Path $Qc $ReferenceContigAuditName),
        "--vcf", "SNP=$($inputs.snp.FullName)",
        "--vcf", "INDEL=$($inputs.indel.FullName)",
        "--vcf", "RAW=$($inputs.raw.FullName)"
    )
    Invoke-BcftoolsNormalize $inputs.snp "Dante_SNP_normalized.vcf.gz" "SNP" $fastaHash $referenceTargetsHash $bcftoolsRuntimeImage
    Invoke-BcftoolsNormalize $inputs.indel "Dante_INDEL_normalized.vcf.gz" "INDEL" $fastaHash $referenceTargetsHash $bcftoolsRuntimeImage
    Invoke-BcftoolsNormalize $inputs.raw "Dante_RAW_normalized.vcf.gz" "RAW" $fastaHash $referenceTargetsHash $bcftoolsRuntimeImage

    $snp = Join-Path $Work "Dante_SNP_normalized.vcf.gz"
    $indel = Join-Path $Work "Dante_INDEL_normalized.vcf.gz"
    $raw = Join-Path $Work "Dante_RAW_normalized.vcf.gz"
    $residual = Join-Path $Work "Dante_RESIDUAL_normalized.vcf.gz"
    $residualState = "$residual.state.sha256"
    $residualFingerprint = "residual-v7-allele-normalized|$(Get-Sha256 $snp)|$(Get-Sha256 $indel)|$(Get-Sha256 $raw)"
    if ((Test-Path -LiteralPath $residual) -and (Test-Path -LiteralPath $residualState) -and
        ((Get-Content -LiteralPath $residualState -Raw).Trim() -eq $residualFingerprint) -and -not $ForceDownload) {
        Write-Host "Residual layer is current; reusing it."
    } else {
        Write-Step "Extract normalized raw residual/complex variants"
        Invoke-ObservedNative -Activity "Extract normalized residual variants" -Program $python -Arguments @(
            (Join-Path $Root "extract_raw_residual.py"),
            "--raw", $raw, "--snp", $snp, "--indel", $indel,
            "--output", $residual,
            "--report", (Join-Path $Qc "Dante_RESIDUAL_preparation_report.txt")
        ) -WatchFile "$residual.part"
        Set-AtomicText $residualState $residualFingerprint
    }
    Write-Host "Prepared normalized SNP, INDEL and residual layers in $Work"
}

function Test-GzipFile([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) { return $false }
    $directory = Split-Path -Parent $Path
    $name = Split-Path -Leaf $Path
    try {
        Invoke-ObservedNative -Activity "Verify gzip integrity: $name" -Program "docker" -Arguments @(
            "run", "--rm", "--mount", "type=bind,source=$directory,target=/check",
            $VepImage, "bash", "-lc", "gzip -t /check/$name"
        ) -WatchFile $Path
        return $true
    } catch {
        return $false
    }
}

function Invoke-FullVepLayer([string]$Layer, [string]$InputName) {
    $input = Join-Path $Work $InputName
    if (-not (Test-Path -LiteralPath $input -PathType Leaf)) { throw "Missing prepared layer: $input" }
    $upper = $Layer.ToUpper()
    $outputName = "Dante_$($upper)_full.vep.vcf.gz"
    $output = Join-Path $FullVepResults $outputName
    $state = "$output.state.sha256"
    $cacheInfo = Join-Path $VepData "homo_sapiens\116_GRCh37\info.txt"
    $referenceIdentity = "$(Get-Sha256 $cacheInfo)|$((Get-Item -LiteralPath (Join-Path $VepFastaData $FastaName)).Length)"
    $fingerprint = "VEP-116-everything-GRCh37|$referenceIdentity|$(Get-Sha256 $input)"
    if ((Test-GzipFile $output) -and (Test-Path -LiteralPath $state) -and
        ((Get-Content -LiteralPath $state -Raw).Trim() -eq $fingerprint) -and -not $ForceDownload) {
        Write-Host "$upper full VEP output is current; reusing it."
    } else {
        $partName = "Dante_$($upper)_full.part.vep.vcf.gz"
        $part = Join-Path $FullVepResults $partName
        if (Test-Path -LiteralPath $part) { Remove-Item -LiteralPath $part -Force }
        Write-Step "Full VEP annotation: $upper"
        Invoke-ObservedNative -Activity "Full VEP annotation: $upper" -Program "docker" -Arguments @(
            "run", "--rm",
            "--mount", "type=bind,source=$Root,target=/work",
            "--mount", "type=bind,source=$VepData,target=/data",
            $VepImage, "vep", "--cache", "--offline", "--cache_version", $VepCacheVersion,
            "--species", "homo_sapiens", "--assembly", "GRCh37", "--dir_cache", "/data",
            "--fasta", "/data/fasta/$FastaName", "--format", "vcf", "--vcf",
            "--everything", "--allele_number", "--show_ref_allele", "--uploaded_allele",
            "--transcript_version", "--check_ref", "--exclude_null_alleles",
            "--fork", "$Threads", "--buffer_size", "5000", "--no_stats",
            "--compress_output", "bgzip", "--force_overwrite",
            "--input_file", "/work/Dante_work/$InputName",
            "--output_file", "/work/Dante_results/full_vep/$partName",
            "--warning_file", "/work/Dante_qc/Dante_$($upper)_full_vep.warnings.txt"
        ) -WatchFile $part
        if (-not (Test-GzipFile $part)) { throw "$upper VEP output failed its gzip integrity check." }
        Move-Item -LiteralPath $part -Destination $output -Force
        Set-AtomicText $state $fingerprint
    }

    $python = Get-Python
    $summary = Join-Path $FullVepResults "Dante_$($upper)_VEP_summary.tsv.gz"
    $rsidVcf = Join-Path $FullVepResults "Dante_$($upper)_rsid.vcf.gz"
    Invoke-ObservedNative -Activity "Postprocess full VEP layer: $upper" -Program $python -Arguments @(
        (Join-Path $Root "vep_postprocess.py"), $output,
        "--rsid-vcf", $rsidVcf, "--summary", $summary,
        "--report", (Join-Path $Qc "Dante_$($upper)_VEP_report.txt"),
        "--layer", $Layer
    ) -WatchFile "$summary.part"
}

function Run-FullVepAnnotation {
    if (-not (Test-VepCacheReady)) { throw "VEP 116 GRCh37 cache is missing. Select item 1 first." }
    if (-not (Test-FastaReady)) { throw "GRCh37 FASTA is missing. Select item 6 first." }
    Assert-Docker
    New-Item -ItemType Directory -Force -Path $Results, $FullVepResults, $Work, $Qc | Out-Null
    Invoke-FullVepLayer "snp" "Dante_SNP_normalized.vcf.gz"
    Invoke-FullVepLayer "indel" "Dante_INDEL_normalized.vcf.gz"
    Invoke-FullVepLayer "residual" "Dante_RESIDUAL_normalized.vcf.gz"
    Write-Host "Full VEP layers are ready in $FullVepResults"
}

function Build-UnifiedCatalog {
    if (-not (Test-ClinVarReady)) { throw "Current ClinVar GRCh37 is missing. Select item 7 first." }
    $python = Get-Python
    $summaries = @(
        (Join-Path $FullVepResults "Dante_SNP_VEP_summary.tsv.gz"),
        (Join-Path $FullVepResults "Dante_INDEL_VEP_summary.tsv.gz"),
        (Join-Path $FullVepResults "Dante_RESIDUAL_VEP_summary.tsv.gz")
    )
    foreach ($summary in $summaries) {
        if (-not (Test-Path -LiteralPath $summary -PathType Leaf)) { throw "Missing VEP summary: $summary. Select item 9 first." }
    }
    New-Item -ItemType Directory -Force -Path $Results, $Qc | Out-Null
    Write-Step "Build unified annotated small-variant catalog"
    $catalogOutput = Join-Path $Results "Dante_small_variants_annotated.tsv.gz"
    Invoke-ObservedNative -Activity "Build unified annotated small-variant catalog" -Program $python -Arguments @(
        (Join-Path $Root "build_unified_catalog.py"),
        "--summary", $summaries[0], "--summary", $summaries[1], "--summary", $summaries[2],
        "--clinvar", (Join-Path $ClinVarData "clinvar.vcf.gz"),
        "--output", $catalogOutput,
        "--report", (Join-Path $Qc "Dante_small_variants_catalog_report.txt"),
        "--manifest", (Join-Path $Results "Dante_pipeline_manifest.json"),
        "--vep-version", "116", "--assembly", "GRCh37.p13"
    ) -WatchFile "$catalogOutput.part"
}

function Build-DeterministicShortlist {
    $catalog = Join-Path $Results "Dante_small_variants_annotated.tsv.gz"
    if (-not (Test-Path -LiteralPath $catalog -PathType Leaf)) { throw "Unified catalog is missing. Select item 10 first." }
    $geneReference = Join-Path $Root "reference_data\CHESS_2.2_genes.tsv"
    if (-not (Test-NonEmptyFile $geneReference)) { throw "Embedded CHESS 2.2 gene reference is missing: $geneReference" }
    $python = Get-Python
    Write-Step "Build deterministic shortlist + genome-wide VUS/unclassified coding sweep (CHESS 2.2 enriched)"
    $shortlistOutput = Join-Path $Results "Dante_priority_candidates.tsv"
    $genomewideOutput = Join-Path $Results "Dante_genomewide_uncertain_coding.tsv"
    Invoke-ObservedNative -Activity "Build deterministic + genome-wide uncertain-coding shortlist" -Program $python -Arguments @(
        (Join-Path $Root "build_shortlist.py"), $catalog,
        "--output", $shortlistOutput,
        "--genomewide-output", $genomewideOutput,
        "--report", (Join-Path $Qc "Dante_priority_candidates_report.txt"),
        "--gene-reference", $geneReference,
        "--unclassified-max-af", "0.01"
    ) -WatchFile "$shortlistOutput.part"
}

function Remove-Stage15OrphanContainers {
    $probe = Invoke-NativeCapture "docker" @(
        "ps", "-aq", "--filter", "label=biodata.dante.stage=15", "--filter", "label=biodata.dante.runtime=reference-alignment"
    )
    if ($probe.ExitCode -ne 0) { return }
    $ids = @($probe.Output -split "`r?`n" | Where-Object { $_ -match '^[0-9a-f]{12,64}$' })
    if ($ids.Count -eq 0) { return }
    Write-Host "Cleaning $($ids.Count) orphan Stage 15 Docker container(s)..." -ForegroundColor Yellow
    foreach ($id in $ids) {
        $cleanup = Invoke-NativeCapture "docker" @("rm", "-f", $id)
        if ($cleanup.ExitCode -ne 0) {
            Write-Host "Warning: could not remove orphan Stage 15 container $id" -ForegroundColor Yellow
        }
    }
}

function Assert-ReferenceAlignmentCompatibility {
    Write-Host "Checking approved BWA and samtools images before processing FASTQ..."
    Write-Stage15Progress 0.0 "TOOL PREFLIGHT | resolve BWA 0.7.19 and samtools 1.22.1 containers"
    $bwaRuntimeImage = Get-BwaRuntimeImage
    $samtoolsRuntimeImage = Get-SamtoolsRuntimeImage

    # IMPORTANT: do not use `sh -c` for these probes. Windows PowerShell native
    # argument quoting can corrupt shell command-substitution strings, which
    # caused the v7.5.6 regression. Probe each executable
    # directly through Docker ENTRYPOINT instead.
    $bwaProbe = Invoke-NativeCapture "docker" @(
        "run", "--rm", "--entrypoint", "bwa", $bwaRuntimeImage
    )
    $bwaVersionMatch = [regex]::Match($bwaProbe.Output, 'Version:\s*(0\.7\.19(?:[^\r\n\s]*)?)')
    if (-not $bwaVersionMatch.Success) {
        throw ("BWA container preflight failed. Expected BWA 0.7.19. " +
            "No input or finalized output was changed. Image: $bwaRuntimeImage. " +
            "Exit: $($bwaProbe.ExitCode). Output: " + (Get-CommandProbeTail $bwaProbe.Output))
    }
    Write-Host "Detected BWA $($bwaVersionMatch.Groups[1].Value)" -ForegroundColor Green

    $samVersionProbe = Invoke-NativeCapture "docker" @(
        "run", "--rm", "--entrypoint", "samtools", $samtoolsRuntimeImage, "--version"
    )
    if (($samVersionProbe.ExitCode -ne 0) -or ($samVersionProbe.Output -notmatch '(?m)^samtools 1\.22\.1(?:\s|$)')) {
        throw ("samtools container preflight failed. Expected samtools 1.22.1. " +
            "No input or finalized output was changed. Image: $samtoolsRuntimeImage. " +
            "Exit: $($samVersionProbe.ExitCode). Output: " + (Get-CommandProbeTail $samVersionProbe.Output))
    }
    Write-Host "Detected samtools 1.22.1" -ForegroundColor Green

    $bashProbe = Invoke-NativeCapture "docker" @(
        "run", "--rm", "--entrypoint", "bash", $samtoolsRuntimeImage, "--version"
    )
    if (($bashProbe.ExitCode -ne 0) -or ($bashProbe.Output -notmatch 'GNU bash')) {
        throw ("samtools image lacks bash required by the current pipefail streaming pipeline. " +
            "Image: $samtoolsRuntimeImage. Output: " + (Get-CommandProbeTail $bashProbe.Output))
    }

    Write-Host "BWA 0.7.19 / samtools 1.22.1 capability check passed." -ForegroundColor Green
    return [pscustomobject]@{ bwa = $bwaRuntimeImage; samtools = $samtoolsRuntimeImage }
}

function Assert-ReferenceCallingCompatibility {
    Write-Host "Checking approved bcftools 1.24 image before reference-aware calling..."
    $bcftoolsRuntimeImage = Get-BcftoolsRuntimeImage

    $versionProbe = Invoke-NativeCapture "docker" @(
        "run", "--rm", "--entrypoint", "bcftools", $bcftoolsRuntimeImage, "--version"
    )
    if (($versionProbe.ExitCode -ne 0) -or ($versionProbe.Output -notmatch '(?m)^bcftools 1\.24(?:\s|$)')) {
        throw ("Reference-calling preflight failed: expected bcftools 1.24. " +
            "Image: $bcftoolsRuntimeImage. Exit: $($versionProbe.ExitCode). Output: " +
            (Get-CommandProbeTail $versionProbe.Output))
    }
    Write-Host "Detected bcftools 1.24" -ForegroundColor Green

    $bashProbe = Invoke-NativeCapture "docker" @(
        "run", "--rm", "--entrypoint", "bash", $bcftoolsRuntimeImage, "--version"
    )
    if (($bashProbe.ExitCode -ne 0) -or ($bashProbe.Output -notmatch 'GNU bash')) {
        throw ("Reference-calling preflight failed: bcftools image lacks bash required by the streaming caller. " +
            "Image: $bcftoolsRuntimeImage. Output: " + (Get-CommandProbeTail $bashProbe.Output))
    }

    # bcftools is pinned to exact 1.24; Stage 16 therefore relies on that
    # versioned command contract instead of fragile parsing of localized/help
    # output for options that are part of the pinned release.
    $bgzipProbe = Invoke-NativeCapture "docker" @(
        "run", "--rm", "--entrypoint", "bgzip", $bcftoolsRuntimeImage, "--version"
    )
    if ($bgzipProbe.ExitCode -ne 0) {
        throw ("Reference-calling preflight failed: bcftools image lacks bgzip required for -C alleles targets. " +
            "Image: $bcftoolsRuntimeImage. Output: " + (Get-CommandProbeTail $bgzipProbe.Output))
    }
    $tabixProbe = Invoke-NativeCapture "docker" @(
        "run", "--rm", "--entrypoint", "tabix", $bcftoolsRuntimeImage, "--version"
    )
    if ($tabixProbe.ExitCode -ne 0) {
        throw ("Reference-calling preflight failed: bcftools image lacks tabix required for -C alleles targets. " +
            "Image: $bcftoolsRuntimeImage. Output: " + (Get-CommandProbeTail $tabixProbe.Output))
    }

    Write-Host "bcftools 1.24 + bgzip/tabix known-alleles calling runtime check passed." -ForegroundColor Green
    return $bcftoolsRuntimeImage
}

function Build-ReferenceTargetPanel {
    Write-Step "Build GRCh37 target panel for reference-aware genotyping"
    if (-not (Test-ClinVarReady)) {
        throw "Current ClinVar GRCh37 is missing. Select item 7 first."
    }
    New-Item -ItemType Directory -Force -Path $Qc, $ReferencePanelSources | Out-Null
    if ((Test-ReferenceTargetsReady) -and -not $ForceDownload) {
        Write-Host "The target panel is current for the installed ClinVar and optional panel sources; reusing it."
        return
    }
    $python = Get-Python
    $clinvar = Join-Path $ClinVarData "clinvar.vcf.gz"
    $targets = Join-Path $Qc $ReferenceTargetsName
    $bed = Join-Path $Qc $ReferenceTargetsBedName
    $report = Join-Path $Qc $ReferenceTargetsReportName
    $coordinateConflicts = Join-Path $Qc $ReferenceCoordinateConflictsName
    $genotypeConflicts = Join-Path $Qc $ReferenceSourceGenotypeConflictsName
    foreach ($part in @("$targets.part", "$bed.part", "$report.part", "$coordinateConflicts.part", "$genotypeConflicts.part")) {
        if (Test-Path -LiteralPath $part) { Remove-Item -LiteralPath $part -Force }
    }
    Invoke-ObservedNative -Activity "Build reference genotype target panel" -Program $python -Arguments @(
        (Join-Path $Root "build_reference_targets.py"),
        "--clinvar", $clinvar,
        "--panel-dir", $ReferencePanelSources,
        "--targets", $targets,
        "--bed", $bed,
        "--report", $report,
        "--coordinate-conflicts", $coordinateConflicts,
        "--genotype-conflicts", $genotypeConflicts
    ) -WatchFile $targets
    foreach ($required in @($targets, $bed, $report, $coordinateConflicts, $genotypeConflicts)) {
        if (-not (Test-NonEmptyFile $required)) { throw "Target-panel builder did not create a complete file: $required" }
    }
    $header = Get-Content -LiteralPath $targets -TotalCount 1 -Encoding UTF8
    if ($header -ne "rsid`tchromosome`tposition`tref`talts`tsources`treference_genotype`treference_genotype_status`treference_genotype_sources") {
        throw "The target panel failed its schema check."
    }
    $panels = @()
    foreach ($panel in @(Get-ReferencePanelFiles)) {
        $identity = Get-FileIdentity $panel.FullName
        $panels += [ordered]@{
            name = $identity.name
            length = $identity.length
            modified_utc_ticks = $identity.modified_utc_ticks
            sha256 = Get-PersistentFileSha256 $panel.FullName
        }
    }
    $clinvarIdentity = Get-FileIdentity $clinvar
    $state = [ordered]@{
        status = "complete"
        schema = "reference-targets-v2"
        assembly = "GRCh37"
        inputs = [ordered]@{
            clinvar = [ordered]@{
                name = $clinvarIdentity.name
                length = $clinvarIdentity.length
                modified_utc_ticks = $clinvarIdentity.modified_utc_ticks
                sha256 = Get-PersistentFileSha256 $clinvar
            }
            panels = $panels
        }
        outputs = [ordered]@{
            targets = [ordered]@{
                name = (Get-Item -LiteralPath $targets).Name
                length = [long](Get-Item -LiteralPath $targets).Length
                modified_utc_ticks = [long](Get-Item -LiteralPath $targets).LastWriteTimeUtc.Ticks
                sha256 = Get-Sha256 $targets
            }
            bed = [ordered]@{
                name = (Get-Item -LiteralPath $bed).Name
                length = [long](Get-Item -LiteralPath $bed).Length
                modified_utc_ticks = [long](Get-Item -LiteralPath $bed).LastWriteTimeUtc.Ticks
                sha256 = Get-Sha256 $bed
            }
            coordinate_conflicts = [ordered]@{
                name = (Get-Item -LiteralPath $coordinateConflicts).Name
                length = [long](Get-Item -LiteralPath $coordinateConflicts).Length
                modified_utc_ticks = [long](Get-Item -LiteralPath $coordinateConflicts).LastWriteTimeUtc.Ticks
                sha256 = Get-Sha256 $coordinateConflicts
            }
            genotype_conflicts = [ordered]@{
                name = (Get-Item -LiteralPath $genotypeConflicts).Name
                length = [long](Get-Item -LiteralPath $genotypeConflicts).Length
                modified_utc_ticks = [long](Get-Item -LiteralPath $genotypeConflicts).LastWriteTimeUtc.Ticks
                sha256 = Get-Sha256 $genotypeConflicts
            }
        }
    }
    Set-AtomicText "$targets.state.json" ($state | ConvertTo-Json -Depth 8)
    if (-not (Test-ReferenceTargetsReady)) { throw "Target panel was created but failed the READY-state validator." }
    Write-Host "Created target panel: $targets"
}

function Run-ReferenceAlignment {
    Write-Step "Align FASTQ reads and retain target-overlapping evidence"
    if (-not (Test-ReferenceTargetsReady)) { throw "Reference target panel is not READY. Select item 13 first." }
    if (-not (Test-FastaReady)) { throw "GRCh37 reference FASTA is missing. Select item 6 first." }
    if ((Test-ReferenceBamReady) -and -not $ForceDownload) {
        Write-Host "The target-only alignment is current for the same FASTQ input, target panel, reference and tool images; reusing it."
        return
    }
    Assert-Docker
    # Resolve all external tool containers before the expensive full FASTQ pass.
    # If a registry is unavailable, fail fast without spending minutes reading 100+ GiB.
    $alignmentTools = Assert-ReferenceAlignmentCompatibility
    $pair = Confirm-FastqInputComplete
    New-Item -ItemType Directory -Force -Path $Work, $BwaIndexData | Out-Null
    $bam = Join-Path $Work $ReferenceBamName
    $alignmentArguments = @(
        (Join-Path $Root "reference_alignment.py"),
        "--root", $Root,
        "--data-dir", $VepData,
        "--index-dir", $BwaIndexData,
        "--fasta-name", $FastaName,
        "--fasta-sha256", (Get-PersistentFileSha256 (Join-Path $VepFastaData $FastaName)),
        "--r1", ([string]$pair.r1),
        "--r1-records", "$([long]$pair.r1_record_count)",
        "--r1-sha256", ([string]$pair.r1_content_sha256)
    )
    if ([string]$pair.read_mode -eq "paired") {
        $alignmentArguments += @(
            "--r2", ([string]$pair.r2),
            "--r2-records", "$([long]$pair.r2_record_count)",
            "--r2-sha256", ([string]$pair.r2_content_sha256)
        )
    }
    $checkpointDirectory = [Environment]::GetEnvironmentVariable("DANTE_REFERENCE_CHECKPOINT_DIR", "Process")
    if (-not $checkpointDirectory) { $checkpointDirectory = Join-Path $Work "reference_alignment_checkpoints" }
    $alignmentArguments += @(
        "--targets-bed", (Join-Path $Qc $ReferenceTargetsBedName),
        "--output", $bam,
        "--bwa-image", ([string]$alignmentTools.bwa),
        "--bwa-contract", $BwaToolContract,
        "--samtools-image", ([string]$alignmentTools.samtools),
        "--samtools-contract", $SamtoolsToolContract,
        "--checkpoint-dir", $checkpointDirectory,
        "--chunk-records", "10000000",
        "--disk-reserve-gb", "15",
        "--threads", "$Threads",
        "--bwa-threads", "$($script:BwaThreads)",
        "--sort-threads", "$($script:SamtoolsSortThreads)",
        "--sort-memory", "$($script:SamtoolsSortMemory)",
        "--index-threads", "$($script:SamtoolsIndexThreads)"
    )
    Remove-Stage15OrphanContainers
    try {
        Invoke-ObservedNative -Activity "Reference-aware $($pair.read_mode)-end FASTQ alignment" -Program (Get-Python) -Arguments $alignmentArguments -WatchFile (Join-Path $Work ($ReferenceBamName -replace '\.bam$', '.part.bam'))
    } finally {
        # reference_alignment.py labels its long-running Docker containers. If Ctrl+C,
        # Python or Docker CLI dies between cleanup points, remove only those Stage 15 runtimes.
        Remove-Stage15OrphanContainers
    }
    foreach ($required in @($bam, "$bam.bai")) {
        if (-not (Test-NonEmptyFile $required)) { throw "Reference alignment did not create a complete file: $required" }
    }
    $r1Identity = Get-FileIdentity ([string]$pair.r1)
    $r2State = $null
    if ([string]$pair.read_mode -eq "paired") {
        $r2Identity = Get-FileIdentity ([string]$pair.r2)
        $r2State = [ordered]@{ name = $r2Identity.name; length = $r2Identity.length; modified_utc_ticks = $r2Identity.modified_utc_ticks; content_sha256 = [string]$pair.r2_content_sha256; record_count = [long]$pair.r2_record_count }
    }
    $state = [ordered]@{
        status = "complete"
        schema = "reference-alignment-v3"
        assembly = "GRCh37"
        read_mode = [string]$pair.read_mode
        alignment_strategy = "resumable-chunks-global-markdup-v1"
        inputs = [ordered]@{
            r1 = [ordered]@{ name = $r1Identity.name; length = $r1Identity.length; modified_utc_ticks = $r1Identity.modified_utc_ticks; content_sha256 = [string]$pair.r1_content_sha256; record_count = [long]$pair.r1_record_count }
            r2 = $r2State
            targets = Get-FileIdentity (Join-Path $Qc $ReferenceTargetsBedName)
            fasta = [ordered]@{
                name = (Get-Item -LiteralPath (Join-Path $VepFastaData $FastaName)).Name
                length = [long](Get-Item -LiteralPath (Join-Path $VepFastaData $FastaName)).Length
                modified_utc_ticks = [long](Get-Item -LiteralPath (Join-Path $VepFastaData $FastaName)).LastWriteTimeUtc.Ticks
                sha256 = Get-PersistentFileSha256 (Join-Path $VepFastaData $FastaName)
            }
        }
        images = [ordered]@{ bwa = [string]$alignmentTools.bwa; samtools = [string]$alignmentTools.samtools; bwa_contract = $BwaToolContract; samtools_contract = $SamtoolsToolContract }
        outputs = [ordered]@{ bam = [ordered]@{ name = (Get-Item -LiteralPath $bam).Name; length = [long](Get-Item -LiteralPath $bam).Length; modified_utc_ticks = [long](Get-Item -LiteralPath $bam).LastWriteTimeUtc.Ticks; sha256 = Get-Sha256 $bam } }
    }
    Set-AtomicText "$bam.state.json" ($state | ConvertTo-Json -Depth 8)
    if (-not (Test-ReferenceBamReady)) { throw "Alignment files were created but failed the READY-state validator." }
    # Only now, after BAM + BAI + provenance/state are all validated, release resumable
    # alignment checkpoints.  If Python or PowerShell crashes earlier, they remain usable.
    if (Test-Path -LiteralPath $checkpointDirectory -PathType Container) {
        $checkpointFiles = @(Get-ChildItem -LiteralPath $checkpointDirectory -Filter "chunk_*.bam" -File -ErrorAction SilentlyContinue)
        $releasedBytes = [long]0
        foreach ($checkpointFile in $checkpointFiles) {
            $releasedBytes += [long]$checkpointFile.Length
            Remove-Item -LiteralPath $checkpointFile.FullName -Force -ErrorAction SilentlyContinue
        }
        if ($checkpointFiles.Count -gt 0) {
            Write-Host ("Released {0} of Stage 15 checkpoint storage after final state validation." -f (Format-ByteSize $releasedBytes)) -ForegroundColor DarkGray
        }
    }
    Write-Host "Created target-only BAM: $bam"
}

function Run-ReferenceGenotypeCalling {
    Write-Step "Call strict reference and alternate genotypes at target rsIDs"
    if (-not (Test-ReferenceBamReady)) { throw "Reference alignment is not READY. Select item 15 first." }
    if ((Test-ReferenceCallsReady) -and -not $ForceDownload) {
        Write-Host "Known-alleles, sex-aware reference calls are current for the same BAM, panel, reference and tool contracts; reusing them."
        return
    }

    # Fail fast on a mixed/partial software update.  This prevents an older
    # reference_calling.py from being launched by a newer PowerShell menu.
    $referenceCallingScript = Join-Path $Root "reference_calling.py"
    if (-not (Test-Path -LiteralPath $referenceCallingScript -PathType Leaf)) {
        throw "reference_calling.py is missing from the pipeline root: $referenceCallingScript"
    }
    $contractProbe = Invoke-NativeCapture (Get-Python) @($referenceCallingScript, "--print-contract")
    $expectedCallsSchema = "reference-calls-v5-known-alleles-ploidy-refq"
    $expectedCheckpointSchema = "reference-call-checkpoints-v4-known-alleles-ploidy-coordinate-safe"
    $contractLines = @(
        $contractProbe.Output -split "`r?`n" |
        Where-Object { $_ -match '\|reference-calls-v5-known-alleles-ploidy-refq\|reference-call-checkpoints-v4-known-alleles-ploidy-coordinate-safe$' }
    )
    $actualCallingContract = if ($contractLines.Count -gt 0) { [string]$contractLines[-1] } else { "" }
    $contractParts = @($actualCallingContract -split '\|')
    $contractReady = (
        ($contractProbe.ExitCode -eq 0) -and
        ($contractParts.Count -eq 3) -and
        ($contractParts[1] -eq $expectedCallsSchema) -and
        ($contractParts[2] -eq $expectedCheckpointSchema)
    )
    if (-not $contractReady) {
        throw (
            "Mixed/incomplete pipeline installation detected. Stage 16 requires analytical contracts " +
            "$expectedCallsSchema / $expectedCheckpointSchema. The package release number is intentionally not part of this compatibility test. " +
            "Extract the current archive directly into the pipeline folder and replace existing program files. Probe output: " +
            ($contractProbe.Output -replace "`r?`n", " | ")
        )
    }
    Write-Host ("Stage 16 analytical contract: {0} / {1} [OK] (caller release {2})" -f $expectedCallsSchema, $expectedCheckpointSchema, $contractParts[0]) -ForegroundColor Green

    Assert-Docker
    $bcftoolsRuntimeImage = Assert-ReferenceCallingCompatibility
    $samtoolsRuntimeImage = Get-SamtoolsRuntimeImage
    $calls = Join-Path $Work $ReferenceCallsName
    $metadata = Join-Path $Qc $ReferenceCallsMetadataName
    $uncallable = Join-Path $Qc $ReferenceCallsUncallableName
    $ploidyFile = Join-Path $Root "reference_data\GRCh37_bcftools_ploidy.tsv"
    if (-not (Test-NonEmptyFile $ploidyFile)) { throw "GRCh37 ploidy definition is missing from the distribution: $ploidyFile" }
    $sampleSex = Get-DanteSampleSexMode
    $bamState = Get-Content -LiteralPath ("$(Join-Path $Work $ReferenceBamName).state.json") -Raw -Encoding UTF8 | ConvertFrom-Json
    $targetsState = Get-Content -LiteralPath ("$(Join-Path $Qc $ReferenceTargetsName).state.json") -Raw -Encoding UTF8 | ConvertFrom-Json

    $callCheckpointDirectory = [Environment]::GetEnvironmentVariable("DANTE_REFERENCE_CALL_CHECKPOINT_DIR", "Process")
    if (-not $callCheckpointDirectory) { $callCheckpointDirectory = Join-Path $Work "reference_call_checkpoints" }
    New-Item -ItemType Directory -Force -Path $callCheckpointDirectory | Out-Null

    $callChunkSites = 250000
    $callChunkSitesEnv = [Environment]::GetEnvironmentVariable("DANTE_REFERENCE_CALL_CHUNK_SITES", "Process")
    if ($callChunkSitesEnv -and ($callChunkSitesEnv -match '^\d+$')) {
        $parsedCallChunkSites = [int64]$callChunkSitesEnv
        if ($parsedCallChunkSites -gt 0 -and $parsedCallChunkSites -le 5000000) {
            $callChunkSites = [int]$parsedCallChunkSites
        }
    }

    Invoke-ObservedNative -Activity "Known-alleles sex-aware target genotype calling" -Program (Get-Python) -Arguments @(
        (Join-Path $Root "reference_calling.py"),
        "--root", $Root,
        "--data-dir", $VepData,
        "--fasta-name", $FastaName,
        "--fasta-sha256", ([string]$bamState.inputs.fasta.sha256),
        "--bam", (Join-Path $Work $ReferenceBamName),
        "--bam-sha256", ([string]$bamState.outputs.bam.sha256),
        "--targets-bed", (Join-Path $Qc $ReferenceTargetsBedName),
        "--targets-sha256", ([string]$targetsState.outputs.bed.sha256),
        "--targets-tsv", (Join-Path $Qc $ReferenceTargetsName),
        "--targets-tsv-sha256", ([string]$targetsState.outputs.targets.sha256),
        "--ploidy-file", $ploidyFile,
        "--ploidy-sha256", (Get-Sha256 $ploidyFile),
        "--sample-sex", $sampleSex,
        "--samtools-image", $samtoolsRuntimeImage,
        "--samtools-contract", $SamtoolsToolContract,
        "--metadata-output", $metadata,
        "--uncallable-output", $uncallable,
        "--output", $calls,
        "--checkpoint-dir", $callCheckpointDirectory,
        "--bcftools-image", $bcftoolsRuntimeImage,
        "--bcftools-contract", $BcftoolsToolContract,
        "--threads", "$($script:BcftoolsThreads)",
        "--min-mapq", "30",
        "--min-baseq", "20",
        "--chunk-sites", "$callChunkSites",
        "--disk-reserve-gb", "5"
    ) -WatchDirectory $callCheckpointDirectory
    foreach ($required in @($calls, "$calls.tbi", $metadata, $uncallable)) {
        if (-not (Test-NonEmptyFile $required)) { throw "Reference calling did not create a complete file: $required" }
    }
    $callMetadata = Get-Content -LiteralPath $metadata -Raw -Encoding UTF8 | ConvertFrom-Json
    if ([string]$callMetadata.schema -ne "reference-calls-v5-known-alleles-ploidy-refq") {
        throw "Reference-call metadata schema is not the expected v5 known-alleles/ploidy/REFQ contract."
    }
    $bamIdentity = Get-FileIdentity (Join-Path $Work $ReferenceBamName)
    $bedIdentity = Get-FileIdentity (Join-Path $Qc $ReferenceTargetsBedName)
    $targetsIdentity = Get-FileIdentity (Join-Path $Qc $ReferenceTargetsName)
    $fastaIdentity = Get-FileIdentity (Join-Path $VepFastaData $FastaName)
    $ploidyIdentity = Get-FileIdentity $ploidyFile
    $state = [ordered]@{
        status = "complete"
        schema = "reference-calls-v5-known-alleles-ploidy-refq"
        assembly = "GRCh37"
        calling_strategy = "resumable-v4-checkpoints-known-alleles-ploidy-plus-PL-reference-quality-v5"
        sample = $callMetadata.sample
        known_alleles = $callMetadata.known_alleles
        call_qc = $callMetadata.call_qc
        inputs = [ordered]@{
            bam = [ordered]@{ name = $bamIdentity.name; length = $bamIdentity.length; modified_utc_ticks = $bamIdentity.modified_utc_ticks; sha256 = [string]$bamState.outputs.bam.sha256 }
            targets_bed = [ordered]@{ name = $bedIdentity.name; length = $bedIdentity.length; modified_utc_ticks = $bedIdentity.modified_utc_ticks; sha256 = [string]$targetsState.outputs.bed.sha256 }
            targets_panel = [ordered]@{ name = $targetsIdentity.name; length = $targetsIdentity.length; modified_utc_ticks = $targetsIdentity.modified_utc_ticks; sha256 = [string]$targetsState.outputs.targets.sha256 }
            fasta = [ordered]@{ name = $fastaIdentity.name; length = $fastaIdentity.length; modified_utc_ticks = $fastaIdentity.modified_utc_ticks; sha256 = [string]$bamState.inputs.fasta.sha256 }
            ploidy_definition = [ordered]@{ name = $ploidyIdentity.name; length = $ploidyIdentity.length; modified_utc_ticks = $ploidyIdentity.modified_utc_ticks; sha256 = Get-Sha256 $ploidyFile }
        }
        image = $bcftoolsRuntimeImage
        samtools_image = $samtoolsRuntimeImage
        samtools_contract = $SamtoolsToolContract
        tool_contract = $BcftoolsToolContract
        filters = [ordered]@{ minimum_mapping_quality = 30; minimum_base_quality = 20 }
        outputs = [ordered]@{
            calls = [ordered]@{ name = (Get-Item -LiteralPath $calls).Name; length = [long](Get-Item -LiteralPath $calls).Length; modified_utc_ticks = [long](Get-Item -LiteralPath $calls).LastWriteTimeUtc.Ticks; sha256 = Get-Sha256 $calls }
            uncallable_reference_loci = [ordered]@{ name = (Get-Item -LiteralPath $uncallable).Name; length = [long](Get-Item -LiteralPath $uncallable).Length; modified_utc_ticks = [long](Get-Item -LiteralPath $uncallable).LastWriteTimeUtc.Ticks; sha256 = Get-Sha256 $uncallable }
        }
    }
    Set-AtomicText "$calls.state.json" ($state | ConvertTo-Json -Depth 10)
    if (-not (Test-ReferenceCallsReady)) { throw "Reference calls were created but failed the v5 READY-state validator." }

    if (Test-Path -LiteralPath $callCheckpointDirectory -PathType Container) {
        $checkpointFiles = @(
            Get-ChildItem -LiteralPath $callCheckpointDirectory -File -ErrorAction SilentlyContinue |
                Where-Object { $_.Name -match '^(call_chunk_.*\.vcf\.gz(?:\.tbi)?|current_targets_.*\.bed|current_alleles_.*\.tsv(?:\.gz(?:\.tbi)?)?|sample_sex\.tsv|checkpoint_manifest\.json(?:\.part)?)$' }
        )
        $releasedBytes = [long]0
        foreach ($checkpointFile in $checkpointFiles) {
            $releasedBytes += [long]$checkpointFile.Length
            Remove-Item -LiteralPath $checkpointFile.FullName -Force -ErrorAction SilentlyContinue
        }
        if ($checkpointFiles.Count -gt 0) {
            Write-Host ("Released {0} of Stage 16 checkpoint storage after v4 state validation." -f (Format-ByteSize $releasedBytes)) -ForegroundColor DarkGray
        }
    }
    Write-Host "Created known-alleles target-call VCF: $calls"
    Write-Host ("Resolved sample sex/ploidy: {0} ({1})" -f $callMetadata.sample.sex_resolved, $callMetadata.sample.sex_source) -ForegroundColor DarkGray
    Write-Host ("Homozygous-reference calls with usable GQ/REFQ: {0:P2}" -f [double]$callMetadata.call_qc.hom_ref_reference_quality_fraction) -ForegroundColor DarkGray
    Write-Host ("ALT-carried calls with native GQ: {0:P2}" -f [double]$callMetadata.call_qc.alt_carried_gq_fraction) -ForegroundColor DarkGray
}

function Build-ReferenceSourceMerge {
    Write-Step "Merge known array genotypes before FASTQ processing"
    if (-not (Test-ReferenceTargetsReady)) { throw "Reference target panel is not READY. Select item 13 first." }
    $base = Join-Path $Results "Dante_23andMe_style.txt.gz"
    if (-not (Test-NonEmptyFile $base)) { throw "The variant-only base export is missing. Select item 2 first." }
    if ((Test-ReferenceSourceMergeReady) -and -not $ForceDownload) {
        Write-Host "Known non-conflicting reference-panel genotypes are already merged into the current summary export; reusing it."
        return
    }
    New-Item -ItemType Directory -Force -Path $Results, $Qc | Out-Null
    $export = Join-Path $Results $ReferenceExportName
    $report = Join-Path $Qc "Dante_reference_export_report.txt"
    $discrepancies = Join-Path $Qc "Dante_reference_export_discrepancies.tsv"
    $provenance = Join-Path $Qc $ReferenceProvenanceName
    Invoke-ObservedNative -Activity "Merge known reference-panel genotypes" -Program (Get-Python) -Arguments @(
        (Join-Path $Root "reference_genotype_export.py"),
        "--targets", (Join-Path $Qc $ReferenceTargetsName),
        "--base", $base,
        "--output", $export,
        "--report", $report,
        "--discrepancies", $discrepancies,
        "--provenance", $provenance,
        "--fastq-source", "NONE"
    ) -WatchFile "$export.part"
    foreach ($required in @($export, $report, $discrepancies, $provenance)) {
        if (-not (Test-NonEmptyFile $required)) { throw "Reference-source merge did not create a complete file: $required" }
    }
    if (-not (Test-GzipWithPython $export)) { throw "The preliminary reference-complete export failed its full gzip integrity check." }
    $targetsState = Get-Content -LiteralPath ("$(Join-Path $Qc $ReferenceTargetsName).state.json") -Raw -Encoding UTF8 | ConvertFrom-Json
    $targetsIdentity = Get-FileIdentity (Join-Path $Qc $ReferenceTargetsName)
    $baseIdentity = Get-FileIdentity $base
    $state = [ordered]@{
        status = "complete"
        schema = "reference-export-v3"
        mode = "reference_sources_only"
        assembly = "GRCh37"
        inputs = [ordered]@{
            targets = [ordered]@{ name = $targetsIdentity.name; length = $targetsIdentity.length; modified_utc_ticks = $targetsIdentity.modified_utc_ticks; sha256 = [string]$targetsState.outputs.targets.sha256 }
            base = [ordered]@{ name = $baseIdentity.name; length = $baseIdentity.length; modified_utc_ticks = $baseIdentity.modified_utc_ticks; sha256 = Get-PersistentFileSha256 $base }
        }
        fastq_source = "NONE"
        outputs = [ordered]@{
            export = [ordered]@{ name = (Get-Item -LiteralPath $export).Name; length = [long](Get-Item -LiteralPath $export).Length; modified_utc_ticks = [long](Get-Item -LiteralPath $export).LastWriteTimeUtc.Ticks; sha256 = Get-Sha256 $export }
            provenance = [ordered]@{ name = (Get-Item -LiteralPath $provenance).Name; length = [long](Get-Item -LiteralPath $provenance).Length; modified_utc_ticks = [long](Get-Item -LiteralPath $provenance).LastWriteTimeUtc.Ticks; sha256 = Get-Sha256 $provenance }
            discrepancies = [ordered]@{ name = (Get-Item -LiteralPath $discrepancies).Name; length = [long](Get-Item -LiteralPath $discrepancies).Length; modified_utc_ticks = [long](Get-Item -LiteralPath $discrepancies).LastWriteTimeUtc.Ticks; sha256 = Get-Sha256 $discrepancies }
        }
    }
    Set-AtomicText "$export.state.json" ($state | ConvertTo-Json -Depth 8)
    if (-not (Test-ReferenceSourceMergeReady)) { throw "Reference-source merge was created but failed the READY-state validator." }
    Write-Host "Created preliminary summary: $export"
    Write-Host "It can be upgraded later with R1/R1+R2 evidence without losing source provenance."
}

function Build-ReferenceCompleteExport {
    Write-Step "Build tiered reference-complete 23andMe/Promethease export"
    if (-not (Test-ReferenceCallsReady)) { throw "Reference calls are not READY. Select item 16 first." }
    $base = Join-Path $Results "Dante_23andMe_style.txt.gz"
    if (-not (Test-NonEmptyFile $base)) { throw "The variant-only base export is missing. Select item 2 first." }
    if ((Test-ReferenceExportReady) -and -not $ForceDownload) {
        Write-Host "The tiered reference-complete export is current for the same base export, calls, panel and A/B/C confidence contract; reusing it."
        return
    }
    New-Item -ItemType Directory -Force -Path $Results, $Qc | Out-Null
    $export = Join-Path $Results $ReferenceExportName
    $report = Join-Path $Qc "Dante_reference_export_report.txt"
    $discrepancies = Join-Path $Qc "Dante_reference_export_discrepancies.tsv"
    $provenance = Join-Path $Qc $ReferenceProvenanceName
    $evidence = Join-Path $Results $ReferenceFastqEvidenceName
    $pair = Get-FastqPairInfo
    $fastqSource = if ([string]$pair.read_mode -eq "single") { "R1" } else { "R1+R2" }
    Invoke-ObservedNative -Activity "Build tiered reference-complete genotype export" -Program (Get-Python) -Arguments @(
        (Join-Path $Root "reference_genotype_export.py"),
        "--targets", (Join-Path $Qc $ReferenceTargetsName),
        "--calls", (Join-Path $Work $ReferenceCallsName),
        "--base", $base,
        "--output", $export,
        "--report", $report,
        "--discrepancies", $discrepancies,
        "--provenance", $provenance,
        "--evidence", $evidence,
        "--fastq-source", $fastqSource,
        "--export-min-tier", $ReferenceExportMinTier,
        "--tier-a-min-dp", "$ReferenceTierAMinDp",
        "--tier-a-min-quality", "$ReferenceTierAMinQuality",
        "--tier-a-min-allele-depth", "$ReferenceTierAMinAlleleDepth",
        "--tier-a-hom-fraction", "$ReferenceTierAHomFraction",
        "--tier-a-het-min-fraction", "$ReferenceTierAHetMinFraction",
        "--tier-a-het-max-fraction", "$ReferenceTierAHetMaxFraction",
        "--tier-b-min-dp", "$ReferenceTierBMinDp",
        "--tier-b-min-quality", "$ReferenceTierBMinQuality",
        "--tier-b-hom-min-support", "$ReferenceTierBHomMinSupport",
        "--tier-b-hom-fraction", "$ReferenceTierBHomFraction",
        "--tier-b-het-min-allele-depth", "$ReferenceTierBHetMinAlleleDepth",
        "--tier-b-het-min-fraction", "$ReferenceTierBHetMinFraction",
        "--tier-b-het-max-fraction", "$ReferenceTierBHetMaxFraction",
        "--tier-c-min-dp", "$ReferenceTierCMinDp",
        "--tier-c-min-quality", "$ReferenceTierCMinQuality"
    ) -WatchFile "$export.part"
    foreach ($required in @($export, $report, $discrepancies, $provenance, $evidence)) {
        if (-not (Test-NonEmptyFile $required)) { throw "Tiered reference export did not create a complete file: $required" }
    }
    if (-not (Test-GzipWithPython $export)) { throw "The reference-complete export failed its full gzip integrity check." }
    if (-not (Test-GzipWithPython $evidence)) { throw "The tiered FASTQ evidence file failed its full gzip integrity check." }
    $callsState = Get-Content -LiteralPath ("$(Join-Path $Work $ReferenceCallsName).state.json") -Raw -Encoding UTF8 | ConvertFrom-Json
    $targetsState = Get-Content -LiteralPath ("$(Join-Path $Qc $ReferenceTargetsName).state.json") -Raw -Encoding UTF8 | ConvertFrom-Json
    $callsIdentity = Get-FileIdentity (Join-Path $Work $ReferenceCallsName)
    $targetsIdentity = Get-FileIdentity (Join-Path $Qc $ReferenceTargetsName)
    $baseIdentity = Get-FileIdentity $base
    $state = [ordered]@{
        status = "complete"
        schema = "reference-export-v4-tiered"
        mode = "fastq_enriched"
        assembly = "GRCh37"
        inputs = [ordered]@{
            calls = [ordered]@{ name = $callsIdentity.name; length = $callsIdentity.length; modified_utc_ticks = $callsIdentity.modified_utc_ticks; sha256 = [string]$callsState.outputs.calls.sha256 }
            targets = [ordered]@{ name = $targetsIdentity.name; length = $targetsIdentity.length; modified_utc_ticks = $targetsIdentity.modified_utc_ticks; sha256 = [string]$targetsState.outputs.targets.sha256 }
            base = [ordered]@{ name = $baseIdentity.name; length = $baseIdentity.length; modified_utc_ticks = $baseIdentity.modified_utc_ticks; sha256 = Get-PersistentFileSha256 $base }
        }
        fastq_source = $fastqSource
        thresholds = [ordered]@{
            export_min_tier = $ReferenceExportMinTier
            tier_a = [ordered]@{ min_dp = $ReferenceTierAMinDp; min_quality = $ReferenceTierAMinQuality; min_allele_depth = $ReferenceTierAMinAlleleDepth; hom_fraction = $ReferenceTierAHomFraction; het_min_fraction = $ReferenceTierAHetMinFraction; het_max_fraction = $ReferenceTierAHetMaxFraction }
            tier_b = [ordered]@{ min_dp = $ReferenceTierBMinDp; min_quality = $ReferenceTierBMinQuality; hom_min_support = $ReferenceTierBHomMinSupport; hom_fraction = $ReferenceTierBHomFraction; het_min_allele_depth = $ReferenceTierBHetMinAlleleDepth; het_min_fraction = $ReferenceTierBHetMinFraction; het_max_fraction = $ReferenceTierBHetMaxFraction }
            tier_c = [ordered]@{ min_dp = $ReferenceTierCMinDp; min_quality = $ReferenceTierCMinQuality; research_only = $true }
        }
        outputs = [ordered]@{
            export = [ordered]@{ name = (Get-Item -LiteralPath $export).Name; length = [long](Get-Item -LiteralPath $export).Length; modified_utc_ticks = [long](Get-Item -LiteralPath $export).LastWriteTimeUtc.Ticks; sha256 = Get-Sha256 $export }
            provenance = [ordered]@{ name = (Get-Item -LiteralPath $provenance).Name; length = [long](Get-Item -LiteralPath $provenance).Length; modified_utc_ticks = [long](Get-Item -LiteralPath $provenance).LastWriteTimeUtc.Ticks; sha256 = Get-Sha256 $provenance }
            discrepancies = [ordered]@{ name = (Get-Item -LiteralPath $discrepancies).Name; length = [long](Get-Item -LiteralPath $discrepancies).Length; modified_utc_ticks = [long](Get-Item -LiteralPath $discrepancies).LastWriteTimeUtc.Ticks; sha256 = Get-Sha256 $discrepancies }
            fastq_evidence = [ordered]@{ name = (Get-Item -LiteralPath $evidence).Name; length = [long](Get-Item -LiteralPath $evidence).Length; modified_utc_ticks = [long](Get-Item -LiteralPath $evidence).LastWriteTimeUtc.Ticks; sha256 = Get-Sha256 $evidence }
        }
    }
    Set-AtomicText "$export.state.json" ($state | ConvertTo-Json -Depth 10)
    if (-not (Test-ReferenceExportReady)) { throw "Tiered reference-complete export was created but failed the READY-state validator." }
    Write-Host "Created: $export"
    Write-Host "Created tiered FASTQ evidence: $evidence" -ForegroundColor DarkGray
    Write-Host "Default export policy: Tier A + B; Tier C is retained for research review only." -ForegroundColor DarkGray
}

function Run-Stage3 {
    if (-not (Test-FastaReady)) { Install-FastaData }
    if (-not (Test-ClinVarReady)) { Install-ClinVarData }
    if (-not (Test-NonEmptyFile (Join-Path $Results "Dante_23andMe_style.txt.gz"))) { Run-RsidExport }
    Build-ReferenceTargetPanel
    Build-ReferenceSourceMerge
    Run-ReferenceAlignment
    Run-ReferenceGenotypeCalling
    Build-ReferenceCompleteExport
}

function Ensure-OpenAiSdk([string]$Python) {
    $available = $true
    try { Invoke-Native $Python @("-c", "import openai") } catch { $available = $false }
    if (-not $available) {
        Write-Host "The OpenAI Python SDK is missing; installing/updating it now..." -ForegroundColor Yellow
        try {
            Invoke-NativeWithRetry $Python @("-m", "pip", "install", "--upgrade", "openai")
        } catch {
            throw ("Stage 4 could not install the optional OpenAI Python SDK. All local genome outputs remain valid; " +
                "fix Python package-network access and retry item 19. " + $_.Exception.Message)
        }
    }
}


function Get-TelomereSampleDate {
    $value = [Environment]::GetEnvironmentVariable("DANTE_SAMPLE_DATE", "Process")
    return $(if ($value) { $value.Trim() } else { "" })
}

function Get-TelomereSampleYear {
    $value = [Environment]::GetEnvironmentVariable("DANTE_SAMPLE_YEAR", "Process")
    return $(if ($value) { $value.Trim() } else { "" })
}

function Get-TelomereSampleMaterial {
    $value = [Environment]::GetEnvironmentVariable("DANTE_SAMPLE_MATERIAL", "Process")
    return $(if ($value) { $value.Trim() } else { "" })
}

function Test-TelomereEstimateReady {
    $report = Join-Path $Results $TelomereReportName
    $metrics = Join-Path $Qc $TelomereMetricsName
    $statePath = Join-Path $Qc $TelomereStateName
    if (-not (Test-NonEmptyFile $report) -or -not (Test-NonEmptyFile $metrics) -or -not (Test-NonEmptyFile $statePath)) { return $false }
    $pair = Get-FastqPairInfo
    if (-not (Test-FastqInputReady $pair)) { return $false }
    try {
        $state = Get-Content -LiteralPath $statePath -Raw -Encoding UTF8 | ConvertFrom-Json
        if ([string]$state.schema -ne "dante-telomere-report-state-v1") { return $false }
        if ([string]$state.read_mode -ne [string]$pair.read_mode) { return $false }
        $r1 = Get-FileIdentity ([string]$pair.r1)
        if (($state.inputs.r1.name -ne $r1.name) -or ([long]$state.inputs.r1.length -ne [long]$r1.length) -or ([long]$state.inputs.r1.modified_utc_ticks -ne [long]$r1.modified_utc_ticks)) { return $false }
        if ([string]$pair.read_mode -eq "paired") {
            if (-not $state.inputs.r2) { return $false }
            $r2 = Get-FileIdentity ([string]$pair.r2)
            if (($state.inputs.r2.name -ne $r2.name) -or ([long]$state.inputs.r2.length -ne [long]$r2.length) -or ([long]$state.inputs.r2.modified_utc_ticks -ne [long]$r2.modified_utc_ticks)) { return $false }
        } elseif ($state.inputs.r2) { return $false }
        if ([string]$state.sample_date -ne (Get-TelomereSampleDate)) { return $false }
        if ([string]$state.sample_year -ne (Get-TelomereSampleYear)) { return $false }
        if ([string]$state.sample_material -ne (Get-TelomereSampleMaterial)) { return $false }
        return $true
    } catch {
        return $false
    }
}

function Get-TelomereMenuStatus {
    $pair = Get-FastqPairInfo
    if (-not (Test-FastqInputReady $pair)) { return "FASTQ MISSING" }
    if (Test-TelomereEstimateReady) { return $(if ([string]$pair.read_mode -eq "paired") { "READY: R1+R2" } else { "READY: R1" }) }
    $report = Join-Path $Results $TelomereReportName
    $metrics = Join-Path $Qc $TelomereMetricsName
    if ((Test-NonEmptyFile $report) -or (Test-NonEmptyFile $metrics)) { return "STALE" }
    return $(if ([string]$pair.read_mode -eq "paired") { "READY TO RUN: R1+R2" } else { "READY TO RUN: R1" })
}

function Run-TelomereEstimate {
    Write-Step "Estimate mean telomere length from raw WGS FASTQ"
    $pair = Confirm-FastqInputComplete -ProgressStage 4
    New-Item -ItemType Directory -Force -Path $Results, $Qc | Out-Null
    $report = Join-Path $Results $TelomereReportName
    $metrics = Join-Path $Qc $TelomereMetricsName
    $statePath = Join-Path $Qc $TelomereStateName
    $handler = Join-Path $Root "telomere_estimate.py"
    if (-not (Test-Path -LiteralPath $handler -PathType Leaf)) { throw "Telomere estimator is missing: $handler" }

    $arguments = @(
        $handler,
        "--r1", ([string]$pair.r1),
        "--r1-sha256", ([string]$pair.r1_content_sha256),
        "--r1-records", "$([long]$pair.r1_record_count)",
        "--metrics", $metrics,
        "--report", $report,
        "--metadata-root", $PatientInput
    )

    # Optional metadata must be appended only when non-empty. In PowerShell an
    # empty-string function return can collapse to $null inside an argument array,
    # leaving a bare --sample-date token and causing argparse to fail before the scan.
    $sampleDate = Get-TelomereSampleDate
    $sampleYear = Get-TelomereSampleYear
    $sampleMaterial = Get-TelomereSampleMaterial
    if ($sampleDate) { $arguments += @("--sample-date", [string]$sampleDate) }
    if ($sampleYear) { $arguments += @("--sample-year", [string]$sampleYear) }
    if ($sampleMaterial) { $arguments += @("--sample-material", [string]$sampleMaterial) }
    if ([string]$pair.read_mode -eq "paired") {
        $arguments += @(
            "--r2", ([string]$pair.r2),
            "--r2-sha256", ([string]$pair.r2_content_sha256),
            "--r2-records", "$([long]$pair.r2_record_count)"
        )
    }
    Invoke-ObservedNative -Activity "Stage 4 telomere FASTQ scan/report" -Program (Get-Python) -Arguments $arguments -WatchFile "$report.part"
    if (-not (Test-NonEmptyFile $report) -or -not (Test-NonEmptyFile $metrics)) { throw "Telomere stage did not create the required TXT report and metrics cache." }

    $r1Identity = Get-FileIdentity ([string]$pair.r1)
    $r2State = $null
    if ([string]$pair.read_mode -eq "paired") {
        $r2Identity = Get-FileIdentity ([string]$pair.r2)
        $r2State = [ordered]@{ name = $r2Identity.name; length = $r2Identity.length; modified_utc_ticks = $r2Identity.modified_utc_ticks; content_sha256 = [string]$pair.r2_content_sha256 }
    }
    $state = [ordered]@{
        schema = "dante-telomere-report-state-v1"
        status = "complete"
        read_mode = [string]$pair.read_mode
        method = "TelSeq-compatible FASTQ estimator v1"
        sample_date = [string]$sampleDate
        sample_year = [string]$sampleYear
        sample_material = [string]$sampleMaterial
        inputs = [ordered]@{
            r1 = [ordered]@{ name = $r1Identity.name; length = $r1Identity.length; modified_utc_ticks = $r1Identity.modified_utc_ticks; content_sha256 = [string]$pair.r1_content_sha256 }
            r2 = $r2State
        }
        outputs = [ordered]@{
            report = [ordered]@{ name = (Get-Item -LiteralPath $report).Name; length = [long](Get-Item -LiteralPath $report).Length; sha256 = Get-Sha256 $report }
            metrics = [ordered]@{ name = (Get-Item -LiteralPath $metrics).Name; length = [long](Get-Item -LiteralPath $metrics).Length; sha256 = Get-Sha256 $metrics }
        }
    }
    Set-AtomicText $statePath ($state | ConvertTo-Json -Depth 8)
    if (-not (Test-TelomereEstimateReady)) { throw "Telomere output was created but failed its READY-state validator." }
    Write-Host "Telomere TXT report: $report" -ForegroundColor Green
    if (-not $sampleDate -and -not $sampleYear) {
        Write-Host "Sample date remains UNKNOWN. Add DANTE_SAMPLE_DATE=YYYY-MM-DD or DANTE_SAMPLE_YEAR=YYYY to .env later; rerunning item 19 will reuse the telomere metrics without rereading FASTQ." -ForegroundColor Yellow
    }
}

function Run-AiResearch {
    if (-not (Test-ApiKeyConfigured)) {
        throw "OPENAI_API_KEY is still a placeholder or missing. Edit .env, then restart this menu."
    }
    $shortlist = Join-Path $Results "Dante_priority_candidates.tsv"
    if (-not (Test-Path -LiteralPath $shortlist -PathType Leaf)) { throw "Deterministic shortlist is missing. Select item 11 first." }
    $python = Get-Python
    Ensure-OpenAiSdk $python
    $maxVariants = [Environment]::GetEnvironmentVariable("OPENAI_MAX_VARIANTS", "Process")
    if (-not $maxVariants) { $maxVariants = "20" }
    $reasoning = [Environment]::GetEnvironmentVariable("OPENAI_REASONING", "Process")
    if (-not $reasoning) { $reasoning = "high" }
    Write-Step "AI evidence research for a small deterministic shortlist"
    Invoke-ObservedNative -Activity "AI evidence research" -Program $python -Arguments @(
        (Join-Path $Root "variant_research_agent.py"), $shortlist,
        "--output-dir", (Join-Path $Results "variant_research"),
        "--max-variants", $maxVariants, "--reasoning", $reasoning
    ) -WatchDirectory (Join-Path $Results "variant_research")
}

function Run-Stage1 {
    Install-VepCacheData
    Install-AnnotSVData
    Run-RsidExport
    Run-StructuralAnnotation
}

function Run-Stage2 {
    if (-not (Test-VepCacheReady)) { Install-VepCacheData }
    Install-FastaData
    Install-ClinVarData
    Prepare-NormalizedLayers
    Run-FullVepAnnotation
    Build-UnifiedCatalog
    Build-DeterministicShortlist
}


function Test-OpenReferenceReady {
    if (-not (Test-Path -LiteralPath $OpenReferenceManifest -PathType Leaf)) { return $false }
    try {
        $manifest = Get-Content -LiteralPath $OpenReferenceManifest -Raw -Encoding UTF8 | ConvertFrom-Json
        if (($manifest.status -ne "complete") -or ($manifest.schema -ne "open-reference-bundle-v1")) { return $false }
        $requiredSources = @("HGNC", "dbSNP", "GIAB", "AlphaMissense", "GenCC", "ClinGen", "Mondo", "HPO")
        foreach ($source in $requiredSources) {
            if ($null -eq $manifest.sources.$source) { return $false }
            foreach ($property in $manifest.sources.$source.artifacts.PSObject.Properties) {
                $artifact = $property.Value
                if (-not $artifact.path) { continue }
                $path = Join-Path $OpenReferenceData ([string]$artifact.path)
                if (-not (Test-Path -LiteralPath $path -PathType Leaf)) { return $false }
                if ([long](Get-Item -LiteralPath $path).Length -ne [long]$artifact.length) { return $false }
            }
        }
        $attr = Join-Path $OpenReferenceData ([string]$manifest.attributions.path)
        if (-not (Test-Path -LiteralPath $attr -PathType Leaf)) { return $false }
        return $true
    } catch {
        return $false
    }
}

function Get-OpenReferenceMenuStatus {
    if (Test-OpenReferenceReady) { return "READY" }
    if (Test-Path -LiteralPath $OpenReferenceManifest -PathType Leaf) { return "INCOMPLETE" }
    if (Test-Path -LiteralPath $OpenReferenceData -PathType Container) {
        $files = @(Get-ChildItem -LiteralPath $OpenReferenceData -File -Recurse -ErrorAction SilentlyContinue)
        if ($files.Count -gt 0) { return "PARTIAL / RESUMABLE" }
    }
    return "NOT RUN"
}

function Install-OpenReferenceKnowledge {
    Write-Step "Install/update open reference knowledge bundle"
    if (-not (Test-Path -LiteralPath $OpenReferenceManager -PathType Leaf)) {
        throw "Open-reference manager is missing: $OpenReferenceManager"
    }
    $python = Get-Python
    $bcftoolsRuntimeImage = Get-BcftoolsRuntimeImage
    $args = @(
        $OpenReferenceManager,
        "install",
        "--root", $Root,
        "--output", $OpenReferenceData,
        "--bcftools-image", $bcftoolsRuntimeImage
    )
    if ($ForceDownload) { $args += "--force" }
    Invoke-ObservedNative -Activity "Open reference knowledge bundle" -Program $python -Arguments $args -WatchDirectory $OpenReferenceData -HeartbeatSeconds 20
    if (-not (Test-OpenReferenceReady)) {
        throw "Open reference knowledge installation finished but the bundle failed its READY validator."
    }
    Write-Host "Installed open reference knowledge bundle: $OpenReferenceData" -ForegroundColor Green
    Write-Host "Attribution ledger: $(Join-Path $OpenReferenceData 'ATTRIBUTIONS.md')" -ForegroundColor DarkGray
    Write-Host "Full dbSNP VCFs are remote-indexed only; the 26/28 GiB files are intentionally not duplicated locally." -ForegroundColor DarkGray
}


function Test-ResearchReferenceReady {
    if (-not (Test-Path -LiteralPath $ResearchReferenceManifest -PathType Leaf)) { return $false }
    try {
        $manifest = Get-Content -LiteralPath $ResearchReferenceManifest -Raw -Encoding UTF8 | ConvertFrom-Json
        if (($manifest.status -ne "complete") -or
            ($manifest.schema -ne "research-reference-bundle-v1") -or
            ($manifest.bundle_scope -ne "RESEARCH_ONLY") -or
            ([bool]$manifest.commercial_allowed)) { return $false }
        foreach ($source in @("PanelApp", "CADD")) {
            if ($null -eq $manifest.sources.$source) { return $false }
            foreach ($property in $manifest.sources.$source.artifacts.PSObject.Properties) {
                $artifact = $property.Value
                if (-not $artifact.path) { continue }
                $path = Join-Path $ResearchReferenceData ([string]$artifact.path)
                if (-not (Test-Path -LiteralPath $path -PathType Leaf)) { return $false }
                if ([long](Get-Item -LiteralPath $path).Length -ne [long]$artifact.length) { return $false }
            }
        }
        return $true
    } catch { return $false }
}

function Get-ResearchReferenceMenuStatus {
    if (Test-ResearchReferenceReady) { return "READY: RESEARCH ONLY" }
    if (Test-Path -LiteralPath $ResearchReferenceManifest -PathType Leaf) { return "INCOMPLETE" }
    if (Test-Path -LiteralPath $ResearchReferenceData -PathType Container) {
        $files = @(Get-ChildItem -LiteralPath $ResearchReferenceData -File -Recurse -ErrorAction SilentlyContinue)
        if ($files.Count -gt 0) { return "PARTIAL / RESUMABLE" }
    }
    return "NOT RUN"
}

function Install-ResearchReferenceKnowledge {
    Write-Step "Install research only knowledge bundle (NO COMMERCIAL USE ALLOWED)"
    if (-not (Test-Path -LiteralPath $ResearchReferenceManager -PathType Leaf)) {
        throw "Research-reference manager is missing: $ResearchReferenceManager"
    }
    $python = Get-Python
    $args = @($ResearchReferenceManager, "install", "--root", $Root, "--output", $ResearchReferenceData)
    if ($ForceDownload) { $args += "--force" }
    Invoke-ObservedNative -Activity "Research-only knowledge bundle" -Program $python -Arguments $args -WatchDirectory $ResearchReferenceData -HeartbeatSeconds 20
    if (-not (Test-ResearchReferenceReady)) {
        throw "Research-only knowledge installation finished but failed its READY validator."
    }
    Write-Host "Installed research-only human knowledge bundle: $ResearchReferenceData" -ForegroundColor Green
    Write-Host "NO COMMERCIAL USE. PanelApp additionally prohibits diagnostic/medical-decision/healthcare-service use without a separate GEL agreement." -ForegroundColor Yellow
    Write-Host "UCSC Chain Files are intentionally not auto-installed under the current EULA; target build reconciliation will use dbSNP dual-build data instead." -ForegroundColor DarkGray
}


function Test-GenomeBrowserReady {
    if (-not (Test-Path -LiteralPath $GenomeBrowserApp -PathType Leaf)) { return $false }
    $bam = Join-Path $Work $ReferenceBamName
    $fasta = Join-Path $VepFastaData $FastaName
    return (
        (Test-NonEmptyFile $bam) -and
        (Test-NonEmptyFile "$bam.bai") -and
        (Test-NonEmptyFile $fasta) -and
        (Test-NonEmptyFile "$fasta.fai")
    )
}

function Get-GenomeBrowserMenuStatus {
    if (Test-GenomeBrowserReady) { return "READY: TARGET BAM" }
    if (-not (Test-FastaReady)) { return "BLOCKED: FASTA NOT READY" }
    if (-not (Test-ReferenceBamReady)) { return "BLOCKED: STAGE 15 BAM NOT READY" }
    return "INCOMPLETE"
}

function Ensure-GenomeBrowserRuntime {
    $python = Get-Python
    $probe = Invoke-NativeCapture $python @(
        "-c",
        "import importlib.metadata as m; print(m.version('flask'))"
    )
    if (($probe.ExitCode -eq 0) -and ($probe.Output -match '^\s*3\.1\.3\s*$')) {
        return
    }

    Write-Host "Installing local genome-browser dependency Flask 3.1.3..." -ForegroundColor DarkCyan
    $install = Invoke-NativeCapture $python @(
        "-m", "pip", "install",
        "--disable-pip-version-check",
        "--user",
        "Flask==3.1.3"
    )
    if ($install.ExitCode -ne 0) {
        throw (
            "Failed to install Flask 3.1.3 for the local genome browser. pip output: " +
            ($install.Output -replace "`r?`n", " | ")
        )
    }

    $verify = Invoke-NativeCapture $python @(
        "-c",
        "import importlib.metadata as m; print(m.version('flask'))"
    )
    if (($verify.ExitCode -ne 0) -or ($verify.Output -notmatch '^\s*3\.1\.3\s*$')) {
        throw "Flask installation completed but the pinned 3.1.3 runtime is not importable by the pipeline Python."
    }
}

function Launch-LocalGenomeBrowser {
    Write-Step "Launch local IGV genome browser"
    if (-not (Test-GenomeBrowserReady)) {
        throw "Genome browser requires the indexed GRCh37 FASTA and Stage 15 target BAM/BAI."
    }
    Ensure-GenomeBrowserRuntime
    $python = Get-Python
    Write-Host "Starting browser on localhost only. Close it with the red 'Stop server' button or Ctrl+C." -ForegroundColor DarkGray
    Write-Host "The current BAM is target-only; empty inter-target regions are not evidence of absent original WGS coverage." -ForegroundColor Yellow
    & $python $GenomeBrowserApp --root $Root --patient-root $PatientRoot --port 8765 --open-browser
    if ($LASTEXITCODE -ne 0) {
        throw "Local genome browser exited with code $LASTEXITCODE."
    }
}

function Select-Mode {
    Show-AutomaticCheck
    Write-Host ""
    Write-Host "Dante WGS local pipeline v8.5.2" -ForegroundColor Cyan
    Write-Host ("Patient: {0}" -f $Patient) -ForegroundColor Green
    Write-Host ("Workspace: {0}" -f $PatientRoot) -ForegroundColor DarkGray
    Write-Host ""
    $vepCacheStatus = Get-VepCacheMenuStatus
    $rsidStatus = Get-RsidExportMenuStatus
    $annotSvDataStatus = Get-AnnotSVMenuStatus
    $annotSvOutputStatus = Get-AnnotSVOutputMenuStatus
    $stage1Status = Get-CompleteStageMenuStatus -Statuses @(
        $vepCacheStatus, $rsidStatus, $annotSvDataStatus, $annotSvOutputStatus
    )

    $fastaStatus = if (Test-FastaReady) { "READY" } else { "MISSING" }
    $clinVarStatus = if (Test-ClinVarReady) { "READY" } else { "MISSING" }
    $normalizedStatus = Get-NormalizedLayersMenuStatus
    $fullVepStatus = Get-FullVepMenuStatus
    $catalogStatus = Get-UnifiedCatalogMenuStatus
    $shortlistStatus = Get-ShortlistMenuStatus
    $stage2Status = Get-CompleteStageMenuStatus -Statuses @(
        $fastaStatus, $clinVarStatus, $normalizedStatus,
        $fullVepStatus, $catalogStatus, $shortlistStatus
    )
    $referenceTargetsStatus = Get-ReferenceTargetsMenuStatus
    $referenceSourceMergeStatus = Get-ReferenceSourceMergeMenuStatus
    $referenceBamStatus = Get-ReferenceBamMenuStatus
    $referenceCallsStatus = Get-ReferenceCallsMenuStatus
    $referenceExportStatus = Get-ReferenceExportMenuStatus
    $stage3Status = Get-CompleteStageMenuStatus -Statuses @(
        $referenceTargetsStatus, $referenceSourceMergeStatus, $referenceBamStatus,
        $referenceCallsStatus, $referenceExportStatus
    )
    $telomereStatus = Get-TelomereMenuStatus
    $allLocalStatus = Get-CompleteStageMenuStatus -Statuses @($stage1Status, $stage2Status)
    $openReferenceStatus = Get-OpenReferenceMenuStatus
    $researchReferenceStatus = Get-ResearchReferenceMenuStatus
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
    Write-Host "STAGE 5 - REFERENCE KNOWLEDGE INSTALLATION (DOWNLOAD / SNAPSHOT ONLY)" -ForegroundColor DarkCyan
    Write-MenuItem "22" "Install/update open reference knowledge bundle" $openReferenceStatus
    Write-MenuItem "23" "Install research only knowledge bundle (NO COMMERCIAL USE ALLOWED)" $researchReferenceStatus
    Write-Host "     NOTE: installation does not itself create findings; a downstream enrichment stage must scan the full genomic outputs before AI/report publication." -ForegroundColor DarkGray
    Write-Host ""
    Write-Host "STAGE 6 - AI EVIDENCE RESEARCH" -ForegroundColor DarkCyan
    Write-MenuItem "20" "Research selected shortlist variants  [API KEY REQUIRED]"
    Write-Host ""
    Write-Host "STAGE 7 - LOCAL GENOME BROWSER" -ForegroundColor DarkCyan
    Write-MenuItem "24" "Launch local IGV.js browser: GRCh37 + BAM reads + calls + ClinVar" $genomeBrowserStatus
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
        "20" { return "ai" }
        "21" { return "local-all" }
        "22" { return "open-reference" }
        "23" { return "research-reference" }
        "24" { return "genome-browser" }
        "0" { return "exit" }
        default { throw "Unknown menu choice: $choice" }
    }
}

Import-DotEnv $EnvFile
Import-DotEnv $PatientEnvFile
$referenceExportTierEnv = [Environment]::GetEnvironmentVariable("DANTE_REFERENCE_EXPORT_MIN_TIER", "Process")
if ($referenceExportTierEnv) {
    $referenceExportTierEnv = $referenceExportTierEnv.Trim().ToUpperInvariant()
    if ($referenceExportTierEnv -notin @("A", "B")) {
        throw "DANTE_REFERENCE_EXPORT_MIN_TIER must be A or B."
    }
    $ReferenceExportMinTier = $referenceExportTierEnv
}
$exitCode = 0
try {
    if ($Mode -eq "menu") { $Mode = Select-Mode }
    switch ($Mode) {
        "exit" { }
        "setup-vep-cache" { Install-VepCacheData }
        "rsid" { Run-RsidExport }
        "setup-annotsv" { Install-AnnotSVData }
        "annotsv" { Run-StructuralAnnotation }
        "stage1" { Run-Stage1 }
        "setup-fasta" { Install-FastaData }
        "setup-clinvar" { Install-ClinVarData }
        "prepare-small" { Prepare-NormalizedLayers }
        "vep-full" { Run-FullVepAnnotation }
        "catalog" { Build-UnifiedCatalog }
        "shortlist" { Build-DeterministicShortlist }
        "stage2" { Run-Stage2 }
        "reference-targets" { Build-ReferenceTargetPanel }
        "reference-source-merge" { Build-ReferenceSourceMerge }
        "reference-align" { Run-ReferenceAlignment }
        "reference-call" { Run-ReferenceGenotypeCalling }
        "reference-export" { Build-ReferenceCompleteExport }
        "stage3" { Run-Stage3 }
        "telomere" { Run-TelomereEstimate }
        "ai" { Run-AiResearch }
        "local-all" { Run-Stage1; Run-Stage2 }
        "open-reference" { Install-OpenReferenceKnowledge }
        "research-reference" { Install-ResearchReferenceKnowledge }
        "genome-browser" { Launch-LocalGenomeBrowser }
    }
    if ($Mode -ne "exit") {
        Write-Step "Done"
        Write-Host "Patient: $Patient"
        Write-Host "Patient workspace: $PatientRoot"
        Write-Host "Results: $Results"
        Write-Host "Quality-control reports: $Qc"
        Write-Host "Original Dante VCF files were not modified."
        if ($Mode -notin @("ai", "research-reference", "genome-browser")) { Write-Host "No API requests were made." }
        if ($Mode -eq "research-reference") { Write-Host "Only public knowledge API requests were made to build the research snapshot; no OpenAI evidence-research API calls were made." }
        if ($Mode -eq "genome-browser") { Write-Host "Genome browser genomic data stayed on localhost. Network access is used only to install/cache pinned browser runtime dependencies when missing." }
    }
} catch {
    $exitCode = 1
    Write-Host ""
    Write-Host "ERROR: $($_.Exception.Message)" -ForegroundColor Red
    Write-Host $_.ScriptStackTrace -ForegroundColor DarkRed
}

if (-not $NoPause) { [void](Read-Host "Press Enter to close") }
exit $exitCode
