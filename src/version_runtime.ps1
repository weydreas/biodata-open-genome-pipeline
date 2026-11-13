# Canonical UI release identity for all runtime profiles.
# The release manifest is the single source of truth; menus must never hardcode
# a release number independently from release/pipeline_manifest.json.

$script:DanteDisplayVersionCache = $null

function Get-DanteDisplayVersion {
    if (-not [string]::IsNullOrWhiteSpace([string]$script:DanteDisplayVersionCache)) {
        return [string]$script:DanteDisplayVersionCache
    }

    $pipelineRoot = Split-Path -Parent $CodeRoot
    $manifestPath = Join-Path $pipelineRoot "release\pipeline_manifest.json"
    if (-not (Test-Path -LiteralPath $manifestPath -PathType Leaf)) {
        return "UNKNOWN"
    }

    try {
        $manifest = Get-Content -LiteralPath $manifestPath -Raw -Encoding UTF8 | ConvertFrom-Json
        $version = [string]$manifest.display_version
        if ([string]::IsNullOrWhiteSpace($version)) { return "UNKNOWN" }
        $script:DanteDisplayVersionCache = $version.Trim()
        return [string]$script:DanteDisplayVersionCache
    } catch {
        return "UNKNOWN"
    }
}
