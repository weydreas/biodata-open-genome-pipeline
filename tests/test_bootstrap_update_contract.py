import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
BOOTSTRAP = ROOT / "src" / "bootstrap_update.ps1"
MANIFEST = ROOT / "release" / "pipeline_manifest.json"
GITIGNORE = ROOT / ".gitignore"


class BootstrapUpdateContractTests(unittest.TestCase):
    def test_powershell_automatic_args_variable_is_not_used_as_parameter(self) -> None:
        text = BOOTSTRAP.read_text(encoding="utf-8-sig")
        self.assertNotIn("[string[]]$Args", text)
        self.assertIn("[string[]]$ManagerArgs", text)
        self.assertIn("@($Manager) + $ManagerArgs", text)

    def test_reference_startup_is_status_only_and_never_auto_installs(self) -> None:
        text = BOOTSTRAP.read_text(encoding="utf-8-sig")
        self.assertIn('function Show-ReferenceBundleStatus', text)
        self.assertIn('$statusArgs = @("status", "--root", $Root, "--output", $output)', text)
        self.assertIn("Invoke-PythonManager $manager $statusArgs", text)
        self.assertIn("startup will not install it automatically", text)
        self.assertNotIn('$installArgs = @("install"', text)
        self.assertNotIn("running resumable manager", text)

    def test_native_stdout_is_not_captured_as_exit_code(self) -> None:
        text = BOOTSTRAP.read_text(encoding="utf-8-sig")
        self.assertIn("& py @parts | Out-Host", text)
        self.assertIn("& python @parts | Out-Host", text)
        self.assertIn("$exitCode = [int]$LASTEXITCODE", text)
        self.assertIn("return $exitCode", text)

    def test_git_checkout_is_native_and_feature_branch_is_never_overwritten(self) -> None:
        text = BOOTSTRAP.read_text(encoding="utf-8-sig")
        self.assertIn('rev-parse", "--is-inside-work-tree"', text)
        self.assertIn('symbolic-ref", "--quiet", "--short", "HEAD"', text)
        self.assertIn("development branch detected", text)
        self.assertIn("running branch code unchanged", text)
        self.assertIn('"merge", "--ff-only", "origin/$Branch"', text)
        self.assertIn("working tree has local changes", text)
        self.assertIn("Sync-LegacyRuntimeCache", text)

    def test_frozen_launchers_are_protected_in_both_git_native_and_legacy_modes(self) -> None:
        text = BOOTSTRAP.read_text(encoding="utf-8-sig")
        self.assertIn('"diff", "--name-only", "HEAD", "origin/$Branch", "--", "start.cmd", "bootstrap_git.cmd"', text)
        self.assertIn("live fast-forward is blocked", text)
        self.assertIn("frozen launchers", text)
        self.assertIn("-ieq 'start.cmd'", text)
        self.assertIn("-ieq 'bootstrap_git.cmd'", text)

        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
        deployed = {str(path).replace("/", "\\").lower() for path in manifest["deployment_paths"]}
        self.assertNotIn("start.cmd", deployed)
        self.assertNotIn("bootstrap_git.cmd", deployed)
        self.assertIn("src\\bootstrap_update.ps1", deployed)
        self.assertIn("tests\\test_bootstrap_update_contract.py", deployed)

    def test_local_patient_selector_state_is_globally_ignored(self) -> None:
        ignored = {
            line.strip()
            for line in GITIGNORE.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        }
        self.assertIn(".active_patient.json", ignored)
        self.assertIn(".pipeline_repo/", ignored)


if __name__ == "__main__":
    unittest.main()
