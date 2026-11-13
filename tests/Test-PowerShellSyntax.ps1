param(
    [string]$ScriptPath = (Join-Path (Join-Path (Split-Path -Parent $PSScriptRoot) "src") "dante_wgs_pipeline.ps1")
)

$tokens = $null
$parseErrors = $null
[void][System.Management.Automation.Language.Parser]::ParseFile(
    $ScriptPath,
    [ref]$tokens,
    [ref]$parseErrors
)
if ($parseErrors.Count -gt 0) {
    foreach ($parseError in $parseErrors) {
        Write-Error ("{0}:{1}:{2}: {3}" -f
            $ScriptPath,
            $parseError.Extent.StartLineNumber,
            $parseError.Extent.StartColumnNumber,
            $parseError.Message)
    }
    exit 1
}
Write-Host "PowerShell parser: OK"
exit 0
