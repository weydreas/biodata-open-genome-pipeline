#!/usr/bin/env python3
"""Research a deliberately small list of prioritized genetic variants.

The agent uses the OpenAI Responses API with hosted web search. It does not
parse or upload a whole genome. Input must be a TSV/CSV shortlist prepared by a
deterministic annotation/prioritization step.

Required environment variable: OPENAI_API_KEY
Optional environment variable: OPENAI_MODEL (default: gpt-5.6)
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
import traceback
from pathlib import Path


ALLOWED_DOMAINS = [
    "pubmed.ncbi.nlm.nih.gov",
    "ncbi.nlm.nih.gov",
    "clinicalgenome.org",
    "gnomad.broadinstitute.org",
    "ensembl.org",
    "pharmgkb.org",
    "cpicpgx.org",
    "clinicaltrials.gov",
    "nature.com",
    "science.org",
    "cell.com",
    "nejm.org",
    "jamanetwork.com",
    "bmj.com",
]


def load_dotenv(path: Path) -> None:
    """Load a small local .env file without adding a python-dotenv dependency."""
    if not path.is_file():
        return
    for raw in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[7:].strip()
        if "=" not in line:
            continue
        name, value = line.split("=", 1)
        name, value = name.strip(), value.strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in {'"', "'"}:
            value = value[1:-1]
        if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", name):
            os.environ.setdefault(name, value)


def sniff_delimiter(path: Path) -> str:
    with path.open("rt", encoding="utf-8-sig", errors="replace") as handle:
        sample = handle.read(8192)
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;").delimiter
    except csv.Error:
        return "\t"


def read_rows(path: Path) -> list[dict[str, str]]:
    delimiter = sniff_delimiter(path)
    with path.open("rt", encoding="utf-8-sig", errors="replace", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter=delimiter)]


def first(row: dict[str, str], *names: str) -> str:
    lowered = {str(k).strip().lower(): str(v).strip() for k, v in row.items() if k is not None}
    for name in names:
        if lowered.get(name.lower()):
            return lowered[name.lower()]
    return ""


def variant_label(row: dict[str, str]) -> str:
    rsid = first(row, "rsid", "rsid_or_id", "id")
    if re.fullmatch(r"rs\d+", rsid, re.IGNORECASE):
        return rsid.lower()
    chrom = first(row, "chrom", "chromosome", "chr")
    pos = first(row, "position", "pos", "start")
    ref = first(row, "ref", "reference")
    alt = first(row, "alt", "alternative")
    if chrom and pos:
        return f"{chrom}:{pos}:{ref or '?'}:{alt or '?'}"
    gene = first(row, "gene", "symbol")
    return gene or "unidentified_variant"


def prompt_for(row: dict[str, str]) -> str:
    label = variant_label(row)
    compact = {k: v for k, v in row.items() if k and v and str(v).strip() not in {"", "."}}
    return f"""Conduct an evidence-controlled research review of the germline genetic finding {label}.

Provided structured annotations (treat as leads, verify important claims):
{json.dumps(compact, ensure_ascii=False, indent=2)}

Return a concise Russian-language Markdown report with these sections:
1. Variant identity: genome build, alleles, gene/transcript, consequence; explicitly flag any unresolved strand/build/allele ambiguity.
2. Database evidence: ClinVar classification and review status, ClinGen if relevant, population frequency including ancestry caveats, functional predictions only as supporting—not clinical—evidence.
3. Human evidence: phenotype associations from primary studies and meta-analyses; distinguish GWAS association, segregation, functional evidence and causal evidence.
4. Genotype interpretation: interpret the exact genotype/zygosity only if the risk/effect allele is verified. State effect size and direction where available.
5. Actionability: CPIC/PharmGKB or recognized clinical guidance if applicable; otherwise say no established clinical action. Do not infer a diagnosis from a variant.
6. Conflicts and uncertainty: contradictory findings, weak sample sizes, ancestry transferability, publication bias, linkage disequilibrium, multiple testing.
7. Verdict: pathogenic/likely pathogenic/VUS/benign risk allele/protective allele/not interpretable, evidence confidence (high/moderate/low), and what evidence would change the conclusion.

Rules:
- Search primary and authoritative sources. Prefer ClinVar/ClinGen, PubMed papers, gnomAD, CPIC and PharmGKB.
- Cite every material claim with a working source link.
- Never equate statistical association with disease causation.
- Never treat a computational predictor alone as clinical evidence.
- If the supplied genotype allele cannot be oriented with confidence, stop the genotype-specific conclusion and explain exactly what is missing.
- Do not give treatment instructions; identify only evidence-backed discussion points for a clinician.
"""


def safe_name(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")
    return value[:100] or "variant"


def response_sources(response_dict: dict) -> list[dict]:
    sources: list[dict] = []
    for item in response_dict.get("output", []):
        if item.get("type") != "web_search_call":
            continue
        action = item.get("action") or {}
        for source in action.get("sources") or []:
            if isinstance(source, dict):
                sources.append(source)
    return sources


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Research a small prioritized variant list via OpenAI web search")
    parser.add_argument("input", type=Path, help="TSV/CSV shortlist; never pass the full genome")
    parser.add_argument("--output-dir", type=Path, default=Path("variant_research"))
    parser.add_argument("--max-variants", type=int, default=20, help="Hard cap for this run (default: 20; max: 200)")
    parser.add_argument("--model", default=os.getenv("OPENAI_MODEL", "gpt-5.6"))
    parser.add_argument(
        "--reasoning",
        choices=["low", "medium", "high", "xhigh", "max"],
        default=os.getenv("OPENAI_REASONING", "high"),
    )
    parser.add_argument("--yes", action="store_true", help="Skip interactive privacy/cost confirmation")
    parser.add_argument("--delay", type=float, default=1.0, help="Seconds between requests")
    return parser.parse_args()


def main() -> int:
    load_dotenv(Path(__file__).resolve().parent / ".env")
    args = parse_args()
    if not args.input.is_file():
        raise FileNotFoundError(args.input)
    if not 1 <= args.max_variants <= 200:
        raise ValueError("--max-variants must be between 1 and 200")
    api_key = os.getenv("OPENAI_API_KEY", "")
    if not api_key or re.fullmatch(r"X+", api_key, re.IGNORECASE):
        raise RuntimeError(
            "OPENAI_API_KEY is missing or still a placeholder. Edit the .env file beside this script."
        )

    rows = read_rows(args.input)
    if not rows:
        raise RuntimeError("Input contains no data rows")
    selected = rows[: args.max_variants]
    print(f"Rows in input: {len(rows):,}")
    print(f"Variants selected for API research: {len(selected):,}")
    print("Only these selected rows—not the full genome—will be sent to the OpenAI API.")
    if not args.yes:
        answer = input("Continue? Type YES: ").strip()
        if answer != "YES":
            print("Cancelled. No API requests were made.")
            return 0

    try:
        from openai import OpenAI
    except ImportError as exc:
        raise RuntimeError("OpenAI SDK is not installed. Run: py -m pip install --upgrade openai") from exc

    args.output_dir.mkdir(parents=True, exist_ok=True)
    client = OpenAI()
    jsonl_path = args.output_dir / "research_results.jsonl"
    failures = 0

    with jsonl_path.open("a", encoding="utf-8") as jsonl:
        for index, row in enumerate(selected, start=1):
            label = variant_label(row)
            report_path = args.output_dir / f"{index:03d}_{safe_name(label)}.md"
            if report_path.exists() and report_path.stat().st_size > 100:
                print(f"[{index}/{len(selected)}] Cached: {label}")
                continue
            print(f"[{index}/{len(selected)}] Researching {label}...", flush=True)
            try:
                response = client.responses.create(
                    model=args.model,
                    reasoning={"effort": args.reasoning},
                    tools=[
                        {
                            "type": "web_search",
                            "filters": {"allowed_domains": ALLOWED_DOMAINS},
                        }
                    ],
                    tool_choice="required",
                    include=["web_search_call.action.sources"],
                    store=False,
                    input=prompt_for(row),
                )
                dumped = response.model_dump() if hasattr(response, "model_dump") else {}
                text = response.output_text
                sources = response_sources(dumped)
                report_path.write_text(
                    f"# {label}\n\n"
                    f"Модель: `{getattr(response, 'model', args.model)}`  \n"
                    f"Response ID: `{getattr(response, 'id', 'unknown')}`\n\n"
                    f"{text}\n\n## Технический список использованных источников\n\n"
                    + "\n".join(
                        f"- {source.get('title') or source.get('url', 'source')}: {source.get('url', '')}"
                        for source in sources
                    )
                    + "\n",
                    encoding="utf-8",
                )
                jsonl.write(
                    json.dumps(
                        {
                            "label": label,
                            "input": row,
                            "report": str(report_path),
                            "sources": sources,
                            "response_id": getattr(response, "id", None),
                            "model": getattr(response, "model", args.model),
                            "usage": dumped.get("usage"),
                            "store": False,
                        },
                        ensure_ascii=False,
                    )
                    + "\n"
                )
                jsonl.flush()
            except Exception as exc:
                failures += 1
                print(f"  ERROR for {label}: {type(exc).__name__}: {exc}", file=sys.stderr)
            if index < len(selected) and args.delay > 0:
                time.sleep(args.delay)

    print(f"Done. Reports: {args.output_dir.resolve()}")
    print(f"Failures: {failures}")
    return 1 if failures else 0


def run_cli() -> None:
    try:
        code = main()
    except SystemExit as exc:
        code = exc.code if isinstance(exc.code, int) else 1
    except KeyboardInterrupt:
        print("\nERROR: cancelled by user.", file=sys.stderr)
        code = 130
    except BaseException as exc:
        print(f"\nERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        traceback.print_exc()
        code = 1
    if code and sys.stdin and sys.stdin.isatty():
        try:
            input("\nPress Enter to close...")
        except (EOFError, KeyboardInterrupt):
            pass
    raise SystemExit(code)


if __name__ == "__main__":
    run_cli()
