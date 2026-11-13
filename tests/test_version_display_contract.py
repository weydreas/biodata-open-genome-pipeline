from __future__ import annotations

import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "release" / "pipeline_manifest.json"
VERSION_RUNTIME = ROOT / "src" / "version_runtime.ps1"
WRAPPER = ROOT / "src" / "dante_wgs_pipeline.ps1"
PROFILE_MENUS = (
    ROOT / "src" / "array_profile_runtime.ps1",
    ROOT / "src" / "genotype_txt_profile_runtime.ps1",
    ROOT / "src" / "vcf_export_runtime.ps1",
)


class VersionDisplayContractTests(unittest.TestCase):
    def test_manifest_is_single_ui_release_source(self) -> None:
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
        self.assertEqual(manifest["display_version"], "9.1.0")

        helper = VERSION_RUNTIME.read_text(encoding="utf-8-sig")
        self.assertIn("function Get-DanteDisplayVersion", helper)
        self.assertIn('release\\pipeline_manifest.json', helper)
        self.assertIn("$manifest.display_version", helper)
        self.assertNotIn('return "9.1.0"', helper)

    def test_every_profile_menu_uses_manifest_version(self) -> None:
        for path in PROFILE_MENUS:
            text = path.read_text(encoding="utf-8-sig")
            with self.subTest(path=path.name):
                self.assertIn("Get-DanteDisplayVersion", text)
                self.assertNotIn('Write-Host "Dante genomic pipeline v9.0.0"', text)
                self.assertNotIn('Write-Host "Dante genomic pipeline v9.1.0"', text)

    def test_seeded_wgs_menu_uses_dynamic_version_in_canonical_wrapper(self) -> None:
        text = WRAPPER.read_text(encoding="utf-8-sig")
        canonical_wrapper = "[switch]$MaterializeOnly" in text
        if not canonical_wrapper:
            self.skipTest(
                "historical ci_materialize_runtime.py intentionally replaces the canonical wrapper; canonical WGS UI is validated before isolated regression and by the Windows production materializer"
            )

        self.assertIn('. (Join-Path $CodeRoot "version_runtime.ps1")', text)
        self.assertIn(
            "$Text.Replace('Dante WGS local pipeline v9.0.0', 'Dante WGS local pipeline v$(Get-DanteDisplayVersion)')",
            text,
        )
        self.assertIn(
            "$Text.Replace('Dante genomic pipeline v9.0.0', 'Dante genomic pipeline v$(Get-DanteDisplayVersion)')",
            text,
        )
        self.assertNotIn(
            "$Text.Replace('Dante WGS local pipeline v9.0.0', 'Dante WGS local pipeline v9.1.0')",
            text,
        )


if __name__ == "__main__":
    unittest.main()
