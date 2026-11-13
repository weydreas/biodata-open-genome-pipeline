from __future__ import annotations

import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

ALLOWED_RUNTIME_SOURCE_EXCEPTIONS = {
    "patients/README.md",
    "reference_data/GRCh37_bcftools_ploidy.tsv",
    "reference_data/CHESS_2.2_REFERENCE.txt",
}
RUNTIME_ONLY_PREFIXES = (
    "patients/",
    "reference_data/",
    "annotation_data/",
    "logs/",
    "results/",
    ".pipeline_repo/",
)
FORBIDDEN_EXACT = {
    ".active_patient.json",
    ".env",
}


def tracked_files() -> list[str]:
    proc = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
    )
    return [item.decode("utf-8") for item in proc.stdout.split(b"\0") if item]


def main() -> int:
    violations: list[str] = []
    for path in tracked_files():
        normalized = path.replace("\\", "/")
        if normalized in FORBIDDEN_EXACT:
            violations.append(normalized)
            continue
        if normalized in ALLOWED_RUNTIME_SOURCE_EXCEPTIONS:
            continue
        if normalized.startswith(RUNTIME_ONLY_PREFIXES):
            violations.append(normalized)

    if violations:
        formatted = "\n".join(f"  - {path}" for path in sorted(violations))
        raise RuntimeError(
            "Mutable runtime/data paths are tracked in the Git index:\n" + formatted
        )

    print(
        "GIT INDEX BOUNDARY PASS | mutable patient/reference/runtime paths are untracked; "
        f"immutable exceptions={len(ALLOWED_RUNTIME_SOURCE_EXCEPTIONS)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
