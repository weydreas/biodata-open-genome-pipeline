#!/usr/bin/env python3
"""Fail CI if the tracked public snapshot contains patient data or obvious credentials."""
from __future__ import annotations
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode('utf-8').split('\0')
DATA_EXT = ('.vcf', '.vcf.gz', '.fastq', '.fastq.gz', '.fq', '.fq.gz', '.bam', '.bai', '.cram', '.crai', '.bcf', '.sqlite', '.sqlite3', '.db')
DATA_DIRS = {'patients', 'annotation_data', 'reference_data', 'results', 'logs'}
TEXT_EXT = {'.py', '.ps1', '.cmd', '.md', '.txt', '.json', '.yml', '.yaml', '.toml', '.sh', '.ini'}
RE_SECRET = [
    re.compile(r'sk-(?:proj-)?[A-Za-z0-9_-]{20,}'),
    re.compile(r'ghp_[A-Za-z0-9]{20,}'),
    re.compile(r'github_pat_[A-Za-z0-9_]{15,}'),
    re.compile(r'AKIA[0-9A-Z]{16}'),
    re.compile(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
]
violations: list[str] = []
for rel in filter(None, FILES):
    normalized = rel.replace('\\', '/')
    parts = normalized.split('/')
    name = parts[-1]
    lower = normalized.lower()
    if (parts[0].lower() in DATA_DIRS or name.lower() == '.env'
            or name.lower().startswith('.env.') and name.lower() != '.env.example'
            or name.lower().startswith('.layout_migration_')
            or lower.endswith(DATA_EXT)):
        violations.append(f'tracked sensitive path: {rel}')
    path = ROOT / rel
    if path.is_file() and path.suffix.lower() in TEXT_EXT and path.stat().st_size < 2_000_000:
        content = path.read_text(encoding='utf-8', errors='replace')
        for pattern in RE_SECRET:
            if pattern.search(content):
                violations.append(f'possible credential pattern: {rel}')
                break
if violations:
    print('PUBLIC REPOSITORY SAFETY CHECK FAILED', file=sys.stderr)
    for message in violations:
        print(message, file=sys.stderr)
    raise SystemExit(1)
print(f'PUBLIC REPOSITORY SAFETY CHECK PASSED: {len(list(filter(None, FILES)))} tracked files')
