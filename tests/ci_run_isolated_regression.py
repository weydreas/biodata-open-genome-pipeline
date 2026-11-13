from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def ignore_runtime_and_vcs(directory: str, names: list[str]) -> set[str]:
    ignored = {".git", ".pipeline_repo", "__pycache__", ".pytest_cache"}
    return {name for name in names if name in ignored}


def run(command: list[str], cwd: Path) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="dante-ci-runtime-") as tmp:
        stage = Path(tmp) / "repo"
        shutil.copytree(ROOT, stage, ignore=ignore_runtime_and_vcs)

        canonical_wrapper = ROOT / "src" / "dante_wgs_pipeline.ps1"
        before = canonical_wrapper.read_bytes()

        # This is the historical runtime-view materializer that the legacy static
        # regression suite was built around. It is intentionally allowed to mutate
        # only the disposable staging copy, never the canonical Actions checkout.
        run([sys.executable, "tests/ci_materialize_runtime.py"], stage)
        run(
            [
                sys.executable,
                "-m",
                "unittest",
                "discover",
                "-s",
                "tests",
                "-p",
                "test_*.py",
                "-v",
            ],
            stage,
        )

        after = canonical_wrapper.read_bytes()
        if before != after:
            raise RuntimeError("Canonical wrapper changed while isolated regression was running")

    print("ISOLATED PYTHON REGRESSION PASS | legacy runtime-view suite preserved; canonical source unchanged")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
