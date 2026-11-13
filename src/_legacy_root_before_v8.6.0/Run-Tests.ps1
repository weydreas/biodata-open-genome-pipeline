$ErrorActionPreference = "Stop"
Write-Host "=== Python unittest discovery ==="
python -m unittest discover -s tests -p "test_*.py" -v
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "`n=== PowerShell parser check ==="
powershell -ExecutionPolicy Bypass -File .\tests\Test-PowerShellSyntax.ps1
exit $LASTEXITCODE
