from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WRAPPER = ROOT / "src" / "dante_wgs_pipeline.ps1"


def main() -> int:
    text = WRAPPER.read_text(encoding="utf-8-sig")
    version_hook = '. (Join-Path $CodeRoot "version_runtime.ps1")'
    report_hook = '. (Join-Path $CodeRoot "report_runtime.ps1")'
    risk_hook = '. (Join-Path $CodeRoot "clinical_risk_runtime.ps1")'
    array_hook = '. (Join-Path $CodeRoot "array_profile_runtime.ps1")'
    array_qc_hook = '. (Join-Path $CodeRoot "array_normalization_qc_runtime.ps1")'
    txt_hook = '. (Join-Path $CodeRoot "genotype_txt_profile_runtime.ps1")'
    vcf_hook = '. (Join-Path $CodeRoot "vcf_export_runtime.ps1")'
    hook_sequence = "\n".join(
        [
            version_hook,
            report_hook,
            risk_hook,
            array_hook,
            array_qc_hook,
            txt_hook,
            vcf_hook,
            "",
            "Import-DotEnv $EnvFile",
        ]
    )

    if version_hook not in text:
        raise RuntimeError("production wrapper does not load version_runtime.ps1")
    wgs_version_replacement = (
        "$Text = $Text.Replace('Dante WGS local pipeline v9.0.0', "
        "'Dante WGS local pipeline v$(Get-DanteDisplayVersion)')"
    )
    if wgs_version_replacement not in text:
        raise RuntimeError("production wrapper does not replace the seeded WGS version label from the release manifest")
    if report_hook not in text:
        raise RuntimeError("production wrapper does not load report_runtime.ps1")
    if risk_hook not in text:
        raise RuntimeError("production wrapper does not load clinical_risk_runtime.ps1")
    if array_hook not in text:
        raise RuntimeError("production wrapper does not load array_profile_runtime.ps1")
    if array_qc_hook not in text:
        raise RuntimeError("production wrapper does not load array_normalization_qc_runtime.ps1")
    if txt_hook not in text:
        raise RuntimeError("production wrapper does not load genotype_txt_profile_runtime.ps1")
    if vcf_hook not in text:
        raise RuntimeError("production wrapper does not load vcf_export_runtime.ps1")
    if hook_sequence not in text:
        raise RuntimeError(
            "runtime hook block must be version -> report -> clinical risk -> array profile -> array normalization QC -> genotype TXT -> generic VCF -> environment init"
        )
    if '"commercial-report" { Run-GenomeReport "commercial" }' not in text:
        raise RuntimeError("production wrapper lost commercial report switch")
    if '"research-report" { Run-GenomeReport "research" }' not in text:
        raise RuntimeError("production wrapper lost research report switch")
    if "[switch]$MaterializeOnly" not in text:
        raise RuntimeError("production wrapper lacks isolated materialize-only contract")
    if "-MaterializeOnly requires -MaterializedRuntimePath" not in text:
        raise RuntimeError("materialize-only mode must require an explicit output path")
    if "Set-Content -LiteralPath $materialized -Value $Text" not in text:
        raise RuntimeError("materialize-only path does not emit the production runtime text")

    print(
        "PRODUCTION WRAPPER PASS | manifest-backed WGS version + version -> report -> clinical risk -> array profile -> array normalization QC -> genotype TXT -> generic VCF -> environment init; isolated materialization enabled"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
