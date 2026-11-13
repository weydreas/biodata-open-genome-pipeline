param(
    [Parameter(Mandatory=$true)][string]$Root,
    [string]$PythonCommand = "py",
    [string]$Repository = "https://github.com/weydreas/biodata-open-genome-pipeline.git",
    [string]$Branch = "main"
)

$ErrorActionPreference = "Stop"
$Root = [System.IO.Path]::GetFullPath($Root)
$Cache = Join-Path $Root ".pipeline_repo"
$ManifestRelative = "release\pipeline_manifest.json"

function Write-Info([string]$Message) { Write-Host "[BOOT] $Message" -ForegroundColor DarkCyan }
function Write-Warn([string]$Message) { Write-Host "[BOOT WARNING] $Message" -ForegroundColor Yellow }

function Invoke-Git([string[]]$GitArgs, [switch]$AllowFailure) {
    & git @GitArgs | Out-Host
    $code = [int]$LASTEXITCODE
    if ($code -ne 0 -and -not $AllowFailure) { throw "git $($GitArgs -join ' ') failed with exit code $code" }
    return $code
}

function Invoke-GitCapture([string[]]$GitArgs, [switch]$AllowFailure) {
    $raw = @(& git @GitArgs 2>&1)
    $code = [int]$LASTEXITCODE
    if ($code -ne 0 -and -not $AllowFailure) { throw "git $($GitArgs -join ' ') failed with exit code $code" }
    return [pscustomobject]@{
        ExitCode = $code
        Output = (($raw | ForEach-Object { [string]$_ }) -join "`n").Trim()
    }
}

function Get-HashOrEmpty([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) { return "" }
    return (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
}

function Test-RootGitCheckout {
    # Legacy installs keep their private checkout in .pipeline_repo, not in $Root.
    # Skip probing the non-repository root: Windows PowerShell 5.1 treats Git stderr
    # as a terminating NativeCommandError while $ErrorActionPreference is Stop.
    if (-not (Test-Path -LiteralPath (Join-Path $Root ".git"))) { return $false }
    if (-not (Get-Command git -ErrorAction SilentlyContinue)) { return $false }
    $probe = Invoke-GitCapture @("-C", $Root, "rev-parse", "--is-inside-work-tree") -AllowFailure
    return (($probe.ExitCode -eq 0) -and ($probe.Output -eq "true"))
}

function Get-RootGitBranch {
    $probe = Invoke-GitCapture @("-C", $Root, "symbolic-ref", "--quiet", "--short", "HEAD") -AllowFailure
    if ($probe.ExitCode -ne 0) { return "" }
    return $probe.Output.Trim()
}

function Test-RootGitDirty {
    $probe = Invoke-GitCapture @("-C", $Root, "status", "--porcelain") -AllowFailure
    if ($probe.ExitCode -ne 0) { return $true }
    return (-not [string]::IsNullOrWhiteSpace($probe.Output))
}

function Get-RootGitSha([string]$Ref) {
    $probe = Invoke-GitCapture @("-C", $Root, "rev-parse", $Ref) -AllowFailure
    if ($probe.ExitCode -ne 0) { return "" }
    return $probe.Output.Trim()
}

function Sync-GitWorkstationCheckout {
    $currentBranch = Get-RootGitBranch
    if (-not $currentBranch) {
        $sha = Get-RootGitSha "HEAD"
        $shortSha = if ($sha.Length -ge 8) { $sha.Substring(0, 8) } else { $sha }
        Write-Info "git-native checkout is detached at $shortSha; using checked-out code unchanged"
        return
    }

    if ($currentBranch -ne $Branch) {
        Write-Info "development branch detected: $currentBranch"
        Write-Info "running branch code unchanged; origin/$Branch will not be deployed over active development"
        return
    }

    Write-Info "git-native workstation on $Branch; checking origin/$Branch"
    $fetchCode = Invoke-Git @("-C", $Root, "fetch", "--quiet", "origin", $Branch) -AllowFailure
    if ($fetchCode -ne 0) {
        Write-Warn "cannot fetch origin/$Branch; continuing with local $Branch"
        return
    }

    if (Test-RootGitDirty) {
        Write-Warn "working tree has local changes; automatic fast-forward is skipped to protect them"
        return
    }

    $headSha = Get-RootGitSha "HEAD"
    $remoteSha = Get-RootGitSha "origin/$Branch"
    if (-not $headSha -or -not $remoteSha) {
        Write-Warn "cannot resolve local/remote $Branch SHA; continuing with local checkout"
        return
    }
    if ($headSha -eq $remoteSha) {
        Write-Info "$Branch is current at $($headSha.Substring(0, 8))"
        return
    }

    $ancestor = Invoke-GitCapture @("-C", $Root, "merge-base", "--is-ancestor", "HEAD", "origin/$Branch") -AllowFailure
    if ($ancestor.ExitCode -eq 0) {
        $frozenDiff = Invoke-GitCapture @(
            "-C", $Root, "diff", "--name-only", "HEAD", "origin/$Branch", "--", "start.cmd", "bootstrap_git.cmd"
        ) -AllowFailure
        if ($frozenDiff.ExitCode -ne 0) {
            Write-Warn "cannot validate frozen launcher delta; automatic fast-forward is skipped"
            return
        }
        if (-not [string]::IsNullOrWhiteSpace($frozenDiff.Output)) {
            Write-Warn "origin/$Branch changes a frozen launcher ($($frozenDiff.Output -replace "`n", ', ')); live fast-forward is blocked"
            Write-Warn "close the launcher and update this workstation explicitly before continuing"
            return
        }

        Write-Info "fast-forwarding local $Branch to origin/$Branch"
        $mergeCode = Invoke-Git @("-C", $Root, "merge", "--ff-only", "origin/$Branch") -AllowFailure
        if ($mergeCode -ne 0) {
            Write-Warn "fast-forward failed; continuing with local $Branch without destructive reset"
            return
        }
        $newSha = Get-RootGitSha "HEAD"
        Write-Info "$Branch updated to $($newSha.Substring(0, 8))"
        return
    }

    $remoteAncestor = Invoke-GitCapture @("-C", $Root, "merge-base", "--is-ancestor", "origin/$Branch", "HEAD") -AllowFailure
    if ($remoteAncestor.ExitCode -eq 0) {
        Write-Warn "local $Branch is ahead of origin/$Branch; bootstrap will not rewrite unpublished commits"
    } else {
        Write-Warn "local $Branch diverged from origin/$Branch; resolve Git history manually before release work"
    }
}

function Deploy-TrackedCode {
    $manifestPath = Join-Path $Cache $ManifestRelative
    if (-not (Test-Path -LiteralPath $manifestPath -PathType Leaf)) {
        throw "Repository deployment manifest is missing: $manifestPath"
    }
    $manifest = Get-Content -LiteralPath $manifestPath -Raw -Encoding UTF8 | ConvertFrom-Json
    $changed = 0
    foreach ($rel in @($manifest.deployment_paths)) {
        # start.cmd and bootstrap_git.cmd are frozen launchers. Replacing a batch file
        # while cmd.exe is still executing it can make the interpreter continue from
        # the same byte offset in a different revision and execute arbitrary text.
        # Legacy non-Git installs therefore update only the mutable payload below src/
        # and other explicitly deployable files; launcher migration is separate.
        $normalizedRel = ([string]$rel).Replace('/', '\')
        if ($normalizedRel -ieq 'start.cmd' -or $normalizedRel -ieq 'bootstrap_git.cmd') {
            Write-Warn "refusing live self-update of frozen launcher: $rel"
            continue
        }
        $src = Join-Path $Cache ([string]$rel)
        $dst = Join-Path $Root ([string]$rel)
        if (-not (Test-Path -LiteralPath $src -PathType Leaf)) { throw "Tracked deployment file missing in repo cache: $rel" }
        $srcHash = Get-HashOrEmpty $src
        $dstHash = Get-HashOrEmpty $dst
        if ($srcHash -eq $dstHash) { continue }
        $parent = Split-Path -Parent $dst
        if ($parent) { New-Item -ItemType Directory -Force -Path $parent | Out-Null }
        $tmp = "$dst.update"
        Copy-Item -LiteralPath $src -Destination $tmp -Force
        Move-Item -LiteralPath $tmp -Destination $dst -Force
        Write-Info "updated $rel"
        $changed++
    }
    Write-Info "legacy code sync complete; changed files: $changed; repo release: $($manifest.display_version)"
}

function Sync-LegacyRuntimeCache {
    try {
        if (-not (Test-Path -LiteralPath (Join-Path $Cache ".git") -PathType Container)) {
            Write-Info "legacy runtime detected; initializing private pipeline code cache from GitHub"
            if (Test-Path -LiteralPath $Cache) { Remove-Item -LiteralPath $Cache -Recurse -Force }
            & git clone --depth 1 --branch $Branch $Repository $Cache | Out-Host
            if ($LASTEXITCODE -ne 0) { throw "initial git clone failed" }
        } else {
            Write-Info "legacy runtime detected; checking GitHub for pipeline updates"
            & git -C $Cache fetch --quiet origin $Branch | Out-Host
            if ($LASTEXITCODE -ne 0) { throw "git fetch failed" }
            & git -C $Cache reset --hard "origin/$Branch" | Out-Null
            if ($LASTEXITCODE -ne 0) { throw "git reset failed" }
        }
        Deploy-TrackedCode
    } catch {
        Write-Warn ("GitHub update unavailable; continuing with existing local code: " + $_.Exception.Message)
    }
}

function Ensure-RepositoryCurrent {
    if ($env:DANTE_SKIP_GIT_UPDATE -eq "1") {
        Write-Info "Git update skipped by compatibility override DANTE_SKIP_GIT_UPDATE=1"
        return
    }
    if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
        Write-Warn "git is not installed; continuing with local pipeline files"
        return
    }
    if (Test-RootGitCheckout) {
        Sync-GitWorkstationCheckout
        return
    }
    Sync-LegacyRuntimeCache
}

function Invoke-PythonManager([string]$Manager, [string[]]$ManagerArgs) {
    $parts = @($Manager) + $ManagerArgs
    if ($PythonCommand -eq "py") {
        & py @parts | Out-Host
    } else {
        & python @parts | Out-Host
    }
    $exitCode = [int]$LASTEXITCODE
    return $exitCode
}

function Show-ReferenceBundleStatus([string]$Name, [string]$ManagerRelative, [string]$OutputRelative, [string]$MenuItem) {
    $manager = Join-Path $Root $ManagerRelative
    $output = Join-Path $Root $OutputRelative
    if (-not (Test-Path -LiteralPath $manager -PathType Leaf)) {
        Write-Warn "$Name manager missing: $manager"
        return
    }
    $statusArgs = @("status", "--root", $Root, "--output", $output)
    $code = Invoke-PythonManager $manager $statusArgs
    if ($code -eq 0) {
        Write-Info "$Name reference bundle READY"
        return
    }
    Write-Warn "$Name reference bundle is not READY; startup will not install it automatically"
    Write-Warn "install/update it explicitly with menu item $MenuItem when that analytical stage is needed"
}

Ensure-RepositoryCurrent

if ($env:DANTE_SKIP_REFERENCE_CHECK -ne "1") {
    Show-ReferenceBundleStatus "open" "src\open_reference_manager.py" "reference_data\open_reference" "22"
    Show-ReferenceBundleStatus "research" "src\research_reference_manager.py" "reference_data\research_reference" "23"
} else {
    Write-Info "reference checks skipped by compatibility override DANTE_SKIP_REFERENCE_CHECK=1"
}
