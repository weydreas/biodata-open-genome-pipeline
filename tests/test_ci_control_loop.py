from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github" / "workflows" / "regression-gate.yml"
DEPENDABOT = ROOT / ".github" / "dependabot.yml"


class CIControlLoopProofTests(unittest.TestCase):
    def test_full_regression_gate_is_stable(self) -> None:
        text = WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("name: python-regression", text)
        self.assertIn("name: windows-powershell", text)
        self.assertIn("name: full-regression", text)
        self.assertIn("needs: [python_regression, windows_powershell]", text)
        self.assertIn('needs.python_regression.result', text)
        self.assertIn('needs.windows_powershell.result', text)

    def test_actions_are_immutable_and_checkout_credentials_do_not_persist(self) -> None:
        text = WORKFLOW.read_text(encoding="utf-8")
        uses = re.findall(r"^\s*uses:\s*([^\s#]+)", text, flags=re.MULTILINE)
        self.assertTrue(uses)
        for action in uses:
            self.assertRegex(action, r"^[^@]+@[0-9a-f]{40}$")
        self.assertNotRegex(text, r"uses:\s*[^\s]+@v\d+")
        self.assertIn("permissions:\n  contents: read", text)
        self.assertGreaterEqual(text.count("persist-credentials: false"), 2)

    def test_materialization_and_windows_execution_are_explicit_isolated_contracts(self) -> None:
        text = WORKFLOW.read_text(encoding="utf-8")
        self.assertNotIn("python tests/ci_materialize_runtime.py\n", text)
        self.assertIn("python tests/ci_run_isolated_regression.py", text)
        self.assertIn("Validate runtime/data Git index boundary", text)
        self.assertIn("Materialize production runtime in isolated staging", text)
        self.assertIn("Parse isolated materialized production runtime", text)
        self.assertIn("Execute Windows production-path selftest", text)
        self.assertIn("Test-WindowsPlatformSelfTest.ps1", text)

    def test_pr_concurrency_does_not_cancel_post_merge_main_proof(self) -> None:
        text = WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("cancel-in-progress: ${{ github.event_name == 'pull_request' }}", text)
        self.assertIn("push:\n    branches: [main]", text)
        self.assertIn("pull_request:\n    branches: [main]", text)

    def test_dependabot_tracks_github_actions_updates(self) -> None:
        text = DEPENDABOT.read_text(encoding="utf-8")
        self.assertIn('package-ecosystem: "github-actions"', text)
        self.assertIn('interval: "weekly"', text)


if __name__ == "__main__":
    unittest.main()
