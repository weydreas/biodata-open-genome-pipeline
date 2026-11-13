param(
    [string]$RepositoryRoot = ""
)

Set-StrictMode -Version 2.0
$ErrorActionPreference = "Stop"

if (-not $RepositoryRoot) {
    $RepositoryRoot = Split-Path -Parent $PSScriptRoot
}
$RepositoryRoot = [System.IO.Path]::GetFullPath($RepositoryRoot)
$BootstrapSource = Join-Path $RepositoryRoot "src\bootstrap_update.ps1"
if (-not (Test-Path -LiteralPath $BootstrapSource -PathType Leaf)) {
    throw "bootstrap_update.ps1 is missing: $BootstrapSource"
}

function Invoke-GitChecked([string]$WorkingDirectory, [string[]]$GitArgs) {
    & git -C $WorkingDirectory @GitArgs | Out-Host
    if ($LASTEXITCODE -ne 0) {
        throw "git -C $WorkingDirectory $($GitArgs -join ' ') failed with exit code $LASTEXITCODE"
    }
}

function Invoke-GitGlobalChecked([string[]]$GitArgs) {
    & git @GitArgs | Out-Host
    if ($LASTEXITCODE -ne 0) {
        throw "git $($GitArgs -join ' ') failed with exit code $LASTEXITCODE"
    }
}

function Commit-All([string]$WorkingDirectory, [string]$Message) {
    Invoke-GitChecked $WorkingDirectory @("add", "-A")
    Invoke-GitChecked $WorkingDirectory @("commit", "-m", $Message)
}

function Run-Bootstrap([string]$WorkingDirectory) {
    $bootstrap = Join-Path $WorkingDirectory "src\bootstrap_update.ps1"
    & powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File $bootstrap -Root $WorkingDirectory -PythonCommand "python" | Out-Host
    if ($LASTEXITCODE -ne 0) {
        throw "bootstrap_update.ps1 failed with exit code $LASTEXITCODE"
    }
}

function Read-Trimmed([string]$Path) {
    return (Get-Content -LiteralPath $Path -Raw -Encoding UTF8).Trim()
}

$TempRoot = Join-Path ([System.IO.Path]::GetTempPath()) ("dante-git-native-" + [Guid]::NewGuid().ToString("N"))
$Origin = Join-Path $TempRoot "origin.git"
$Seed = Join-Path $TempRoot "seed"
$Work = Join-Path $TempRoot "work"

try {
    New-Item -ItemType Directory -Force -Path $TempRoot | Out-Null
    Invoke-GitGlobalChecked @("init", "--bare", $Origin)
    Invoke-GitGlobalChecked @("init", $Seed)
    Invoke-GitChecked $Seed @("checkout", "-b", "main")
    Invoke-GitChecked $Seed @("config", "user.name", "Dante CI")
    Invoke-GitChecked $Seed @("config", "user.email", "dante-ci@example.invalid")

    New-Item -ItemType Directory -Force -Path (Join-Path $Seed "src") | Out-Null
    Copy-Item -LiteralPath $BootstrapSource -Destination (Join-Path $Seed "src\bootstrap_update.ps1") -Force
    [System.IO.File]::WriteAllText((Join-Path $Seed "start.cmd"), "@echo off`r`necho launcher-v1`r`n", (New-Object System.Text.UTF8Encoding($false)))
    [System.IO.File]::WriteAllText((Join-Path $Seed "payload.txt"), "v1`n", (New-Object System.Text.UTF8Encoding($false)))
    Commit-All $Seed "initial"
    Invoke-GitChecked $Seed @("remote", "add", "origin", $Origin)
    Invoke-GitChecked $Seed @("push", "-u", "origin", "main")
    Invoke-GitChecked $Origin @("symbolic-ref", "HEAD", "refs/heads/main")

    Invoke-GitGlobalChecked @("clone", $Origin, $Work)
    Invoke-GitChecked $Work @("config", "user.name", "Dante CI")
    Invoke-GitChecked $Work @("config", "user.email", "dante-ci@example.invalid")

    Write-Host "[SELFTEST] git-native main auto fast-forward"
    [System.IO.File]::WriteAllText((Join-Path $Seed "payload.txt"), "v2`n", (New-Object System.Text.UTF8Encoding($false)))
    Commit-All $Seed "main v2"
    Invoke-GitChecked $Seed @("push", "origin", "main")
    Run-Bootstrap $Work
    if ((Read-Trimmed (Join-Path $Work "payload.txt")) -ne "v2") {
        throw "git-native main did not fast-forward to origin/main"
    }
    if (Test-Path -LiteralPath (Join-Path $Work ".pipeline_repo")) {
        throw "git-native checkout unexpectedly created .pipeline_repo"
    }

    Write-Host "[SELFTEST] frozen launcher delta blocks live fast-forward"
    $beforeFrozenSha = (& git -C $Work rev-parse HEAD).Trim()
    [System.IO.File]::WriteAllText((Join-Path $Seed "start.cmd"), "@echo off`r`necho launcher-v2`r`n", (New-Object System.Text.UTF8Encoding($false)))
    [System.IO.File]::WriteAllText((Join-Path $Seed "payload.txt"), "v3`n", (New-Object System.Text.UTF8Encoding($false)))
    Commit-All $Seed "main v3 changes frozen launcher"
    Invoke-GitChecked $Seed @("push", "origin", "main")
    Run-Bootstrap $Work
    $afterFrozenSha = (& git -C $Work rev-parse HEAD).Trim()
    if ($afterFrozenSha -ne $beforeFrozenSha) {
        throw "bootstrap advanced main despite a frozen launcher delta"
    }
    if ((Read-Trimmed (Join-Path $Work "payload.txt")) -ne "v2") {
        throw "payload changed despite frozen launcher block"
    }
    if ((Read-Trimmed (Join-Path $Work "start.cmd")) -notmatch "launcher-v1") {
        throw "start.cmd changed during blocked live fast-forward"
    }

    Write-Host "[SELFTEST] main resumes automatic updates when frozen launcher final state is unchanged"
    [System.IO.File]::WriteAllText((Join-Path $Seed "start.cmd"), "@echo off`r`necho launcher-v1`r`n", (New-Object System.Text.UTF8Encoding($false)))
    [System.IO.File]::WriteAllText((Join-Path $Seed "payload.txt"), "v4`n", (New-Object System.Text.UTF8Encoding($false)))
    Commit-All $Seed "main v4 restores frozen launcher"
    Invoke-GitChecked $Seed @("push", "origin", "main")
    Run-Bootstrap $Work
    if ((Read-Trimmed (Join-Path $Work "payload.txt")) -ne "v4") {
        throw "main did not resume safe automatic fast-forward"
    }

    Write-Host "[SELFTEST] feature branch is protected from origin/main deployment"
    Invoke-GitChecked $Work @("switch", "-c", "feat/local-test")
    [System.IO.File]::WriteAllText((Join-Path $Work "payload.txt"), "feature-local`n", (New-Object System.Text.UTF8Encoding($false)))
    Commit-All $Work "feature local"

    [System.IO.File]::WriteAllText((Join-Path $Seed "payload.txt"), "v5`n", (New-Object System.Text.UTF8Encoding($false)))
    Commit-All $Seed "main v5"
    Invoke-GitChecked $Seed @("push", "origin", "main")
    Run-Bootstrap $Work

    $branch = (& git -C $Work branch --show-current).Trim()
    if ($branch -ne "feat/local-test") {
        throw "bootstrap changed the active development branch: $branch"
    }
    if ((Read-Trimmed (Join-Path $Work "payload.txt")) -ne "feature-local") {
        throw "origin/main was deployed over active feature-branch code"
    }
    if (Test-Path -LiteralPath (Join-Path $Work ".pipeline_repo")) {
        throw "feature-branch git checkout unexpectedly created .pipeline_repo"
    }

    Write-Host "GIT-NATIVE BOOTSTRAP SELFTEST PASS | main auto-FF + frozen-launcher guard + feature isolation + no cache" -ForegroundColor Green
    exit 0
} finally {
    Remove-Item -LiteralPath $TempRoot -Recurse -Force -ErrorAction SilentlyContinue
}
