$ErrorActionPreference = "Stop"
$CodeRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$InstallRoot = Split-Path -Parent $CodeRoot
$Tests = Join-Path $InstallRoot "tests"

Write-Host "=== Python unittest discovery ==="
python -m unittest discover -s $Tests -p "test_*.py" -v
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "`n=== PowerShell parser check ==="
powershell -ExecutionPolicy Bypass -File (Join-Path $Tests "Test-PowerShellSyntax.ps1")
exit $LASTEXITCODE
