#!/usr/bin/env python3
"""Build BWA index and create a duplicate-marked target-only BAM.

v8.2.3 alignment architecture:
- align the complete FASTQ input in resumable chunks;
- single-end input uses double-buffered file staging: while BWA processes chunk N,
  Python decompresses/stages chunk N+1 to hide FASTQ input latency;
- each completed chunk is atomically published as a coordinate-sorted full-genome BAM;
- single-end chunks do not need name collation/fixmate;
- paired-end chunks use samtools collate -> fixmate -> coordinate sort;
- after every input chunk is complete, samtools merge streams all coordinate-sorted
  checkpoints through one global markdup pass and only then retains target-overlapping reads.

This keeps the analytical contract of a whole-input alignment while making expensive BWA
work restartable from the last completed chunk and avoiding the previous whole-genome
query-name sort scratch explosion.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import BinaryIO


INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")
CHECKPOINT_SCHEMA = "reference-alignment-checkpoints-v1"
DEFAULT_CHUNK_RECORDS = 10_000_000
STREAM_WRITE_BUFFER = 8 * 1024 * 1024
GIB = 1024 ** 3

DOCKER_STAGE_LABEL = "biodata.dante.stage=15"
DOCKER_RUNTIME_LABEL = "biodata.dante.runtime=reference-alignment"


class TransientDockerError(RuntimeError):
    """Docker daemon/transport failure for which the current chunk may be retried."""


def docker_engine_ready(timeout_seconds: int = 15) -> bool:
    try:
        completed = subprocess.run(
            ["docker", "info", "--format", "{{.ServerVersion}}"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=timeout_seconds,
            check=False,
        )
        return completed.returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def wait_for_docker_engine(max_wait_seconds: int = 120) -> bool:
    deadline = time.monotonic() + max_wait_seconds
    while time.monotonic() < deadline:
        if docker_engine_ready():
            return True
        time.sleep(5)
    return False


def managed_container_args(name: str) -> list[str]:
    return [
        "--name", name,
        "--label", DOCKER_STAGE_LABEL,
        "--label", DOCKER_RUNTIME_LABEL,
    ]


def force_remove_container(name: str) -> None:
    if not name:
        return
    try:
        subprocess.run(
            ["docker", "rm", "-f", name],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=20,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        pass


def cleanup_stale_stage15_containers() -> int:
    """Remove only containers created/labeled by the v8 Stage 15 runtime."""
    try:
        completed = subprocess.run(
            [
                "docker", "ps", "-aq",
                "--filter", f"label={DOCKER_STAGE_LABEL}",
                "--filter", f"label={DOCKER_RUNTIME_LABEL}",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=20,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return 0
    if completed.returncode != 0:
        return 0
    ids = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    for container_id in ids:
        force_remove_container(container_id)
    if ids:
        print(f"[RECOVERY] Removed {len(ids)} orphan Stage 15 Docker container(s).", file=sys.stderr, flush=True)
    return len(ids)


def run(command: list[str]) -> None:
    completed = subprocess.run(command, check=False)
    if completed.returncode:
        raise RuntimeError(f"command exited with {completed.returncode}: {command[0]}")


def capture(command: list[str]) -> str:
    completed = subprocess.run(command, check=False, stdout=subprocess.PIPE, text=True)
    if completed.returncode:
        raise RuntimeError(f"command exited with {completed.returncode}: {command[0]}")
    return completed.stdout.strip()


def mount(source: Path, target: str) -> list[str]:
    return ["--mount", f"type=bind,source={source.resolve()},target={target}"]


def build_bwa_command(
    index_dir: Path,
    checkpoint_dir: Path,
    bwa_image: str,
    threads: int,
    read_group: str,
    r1_name: str,
    r2_name: str | None,
    container_name: str | None = None,
) -> list[str]:
    command = [
        "docker",
        "run",
        "--rm",
        *(managed_container_args(container_name) if container_name else []),
        *mount(index_dir, "/index"),
        *mount(checkpoint_dir, "/checkpoint"),
        "--entrypoint",
        "bwa",
        bwa_image,
        "mem",
        "-t",
        str(threads),
        "-R",
        read_group,
        "/index/GRCh37",
        f"/checkpoint/{r1_name}",
    ]
    if r2_name is not None:
        command.append(f"/checkpoint/{r2_name}")
    return command


def build_bwa_stdin_command(
    index_dir: Path,
    bwa_image: str,
    threads: int,
    read_group: str,
    container_name: str | None = None,
) -> list[str]:
    """Run single-end BWA MEM with FASTQ supplied on stdin.

    This avoids writing an uncompressed ``input_chunk_*.fastq`` to the checkpoint
    volume and then reading it back.  ``-`` is BWA's stdin input path.
    """
    return [
        "docker",
        "run",
        "--rm",
        *(managed_container_args(container_name) if container_name else []),
        "-i",
        *mount(index_dir, "/index"),
        "--entrypoint",
        "bwa",
        bwa_image,
        "mem",
        "-t",
        str(threads),
        "-R",
        read_group,
        "/index/GRCh37",
        "-",
    ]


def ensure_bwa_index(
    data_dir: Path,
    index_dir: Path,
    fasta_name: str,
    fasta_sha256: str,
    bwa_image: str,
    bwa_contract: str,
) -> Path:
    index_dir.mkdir(parents=True, exist_ok=True)
    final_prefix = index_dir / "GRCh37"
    state_path = index_dir / "GRCh37.index.state.json"
    try:
        state = json.loads(state_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        state = {}
    index_ready = all(
        Path(str(final_prefix) + suffix).is_file()
        and Path(str(final_prefix) + suffix).stat().st_size > 0
        for suffix in INDEX_SUFFIXES
    )
    if (
        index_ready
        and state.get("status") == "complete"
        and state.get("fasta_sha256") == fasta_sha256
        and state.get("bwa_contract") == bwa_contract
    ):
        return final_prefix

    part_prefix = index_dir / "GRCh37.part"
    for suffix in INDEX_SUFFIXES:
        part = Path(str(part_prefix) + suffix)
        if part.exists():
            part.unlink()
    run(
        [
            "docker",
            "run",
            "--rm",
            *mount(data_dir, "/data"),
            *mount(index_dir, "/index"),
            "--entrypoint",
            "bwa",
            bwa_image,
            "index",
            "-p",
            "/index/GRCh37.part",
            f"/data/fasta/{fasta_name}",
        ]
    )
    missing = [suffix for suffix in INDEX_SUFFIXES if not Path(str(part_prefix) + suffix).is_file()]
    if missing:
        raise RuntimeError(f"BWA index is incomplete; missing: {', '.join(missing)}")
    for suffix in INDEX_SUFFIXES:
        part = Path(str(part_prefix) + suffix)
        if part.stat().st_size <= 0:
            raise RuntimeError(f"BWA index component is empty: {part.name}")
    for suffix in INDEX_SUFFIXES:
        os.replace(Path(str(part_prefix) + suffix), Path(str(final_prefix) + suffix))
    state_part = Path(str(state_path) + ".part")
    state_part.write_text(
        json.dumps(
            {
                "status": "complete",
                "schema": "bwa-index-v2",
                "fasta_sha256": fasta_sha256,
                "bwa_contract": bwa_contract,
                "bwa_image_used": bwa_image,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    os.replace(state_part, state_path)
    return final_prefix


def format_duration(seconds: float | None) -> str:
    if seconds is None or seconds < 0 or seconds == float("inf"):
        return "unknown"
    total = max(0, int(seconds))
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def format_reads(value: int) -> str:
    if value >= 1_000_000_000:
        return f"{value / 1_000_000_000:.2f}B"
    if value >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if value >= 1_000:
        return f"{value / 1_000:.1f}k"
    return str(value)


def format_gib(value: int) -> str:
    return f"{value / GIB:.1f} GiB"


def stage15_progress(percent: float, detail: str) -> None:
    bounded = min(100.0, max(0.0, percent))
    print(f"[OVERALL] Stage 15 | {bounded:.2f}% | {detail}", file=sys.stderr, flush=True)


class AlignmentProgress:
    """Translate only completed BWA MEM batches into global read-count progress."""

    PROCESS_RE = re.compile(
        r"\[M::mem_process_seqs\]\s+Processed\s+([0-9,]+)\s+reads\b",
        re.IGNORECASE,
    )

    def __init__(self, total_sequences: int, completed_sequences: int = 0) -> None:
        self.total = max(1, total_sequences)
        self.base = max(0, min(completed_sequences, self.total))
        self.processed = self.base
        self.session_processed = 0
        self.started = time.monotonic()
        self.last_emit = 0.0

    def handle_line(self, line: bytes) -> None:
        text = line.decode("utf-8", errors="replace")
        match = self.PROCESS_RE.search(text)
        if not match:
            return
        batch = int(match.group(1).replace(",", ""))
        self.session_processed += batch
        global_processed = min(self.total, self.base + self.session_processed)
        self.processed = global_processed
        now = time.monotonic()
        if now - self.last_emit < 3.0 and global_processed < self.total:
            return
        self.last_emit = now
        elapsed = max(0.001, now - self.started)
        speed = self.session_processed / elapsed
        remaining = max(0, self.total - global_processed)
        eta = remaining / speed if speed > 0 else None
        ratio = global_processed / self.total
        overall = 8.0 + 82.0 * ratio
        stage15_progress(
            overall,
            f"ALIGNMENT {ratio * 100:.2f}% | {format_reads(global_processed)} / "
            f"{format_reads(self.total)} reads | {format_reads(int(speed))} reads/s | "
            f"ETA {format_duration(eta)}",
        )


def pump_bwa(stream, progress: AlignmentProgress) -> None:  # type: ignore[no-untyped-def]
    try:
        while True:
            line = stream.readline()
            if not line:
                return
            sys.stderr.buffer.write(line)
            sys.stderr.buffer.flush()
            progress.handle_line(line)
    finally:
        stream.close()


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(part, path)


def clean_transient_files(checkpoint_dir: Path) -> None:
    for pattern in (
        "input_chunk_*.fastq",
        "input_chunk_*.fastq.part",
        "chunk_*.part.bam",
        "collate_*.bam",
        "collate_*.tmp*",
        "coord_*.tmp*",
        "coord_*.bam",
    ):
        for path in checkpoint_dir.glob(pattern):
            if path.is_file():
                try:
                    path.unlink()
                except OSError:
                    pass
    tmp = checkpoint_dir / "tmp"
    if tmp.is_dir():
        for path in tmp.iterdir():
            if path.is_file():
                try:
                    path.unlink()
                except OSError:
                    pass


def checkpoint_identity(args: argparse.Namespace, read_mode: str) -> dict:
    return {
        "schema": CHECKPOINT_SCHEMA,
        "read_mode": read_mode,
        "r1_sha256": args.r1_sha256.lower(),
        "r2_sha256": (args.r2_sha256 or "").lower(),
        "r1_records": int(args.r1_records),
        "r2_records": int(args.r2_records or 0),
        "fasta_sha256": args.fasta_sha256.lower(),
        "bwa_contract": args.bwa_contract,
        "samtools_contract": args.samtools_contract,
        "checkpoint_format": "coordinate-bam-v1",
    }


def load_manifest(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, json.JSONDecodeError):
        return {}


def manifest_matches(manifest: dict, identity: dict) -> bool:
    """Match analytical checkpoint identity while allowing chunk-size policy changes.

    v7.5.11 manifests stored ``chunk_records`` inside identity.  Chunk size only
    controls checkpoint cadence; it does not change alignment semantics.  Ignore
    that legacy field so completed 5M-record checkpoints remain resumable when a
    later run switches to 20M checkpoints.
    """
    if manifest.get("schema") != CHECKPOINT_SCHEMA:
        return False
    stored = manifest.get("identity")
    if not isinstance(stored, dict):
        return False
    stored = dict(stored)
    stored.pop("chunk_records", None)
    expected = dict(identity)
    expected.pop("chunk_records", None)
    return stored == expected


def reset_stale_checkpoints(checkpoint_dir: Path, manifest_path: Path) -> None:
    stale = list(checkpoint_dir.glob("chunk_*.bam"))
    if stale or manifest_path.exists():
        print(
            "[CHECKPOINT] Existing checkpoint set is incompatible with the current FASTQ/reference/tool identity; "
            "clearing derived checkpoint files only.",
            file=sys.stderr,
            flush=True,
        )
    for path in stale:
        try:
            path.unlink()
        except OSError:
            pass
    try:
        manifest_path.unlink()
    except OSError:
        pass
    clean_transient_files(checkpoint_dir)


def valid_completed_chunks(manifest: dict, checkpoint_dir: Path) -> list[dict]:
    result: list[dict] = []
    expected_index = 0
    for item in manifest.get("completed_chunks", []):
        if not isinstance(item, dict) or int(item.get("index", -1)) != expected_index:
            break
        bam = checkpoint_dir / str(item.get("bam_name", ""))
        try:
            size = bam.stat().st_size
        except OSError:
            break
        if size <= 0 or size != int(item.get("bam_length", -1)):
            break
        result.append(item)
        expected_index += 1
    return result


def read_fastq_record(handle: BinaryIO) -> tuple[bytes, bytes, bytes, bytes] | None:
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline()
    plus = handle.readline()
    quality = handle.readline()
    if not sequence or not plus or not quality:
        raise RuntimeError("FASTQ became truncated while creating an alignment chunk")
    return header, sequence, plus, quality


def write_record(handle: BinaryIO, record: tuple[bytes, bytes, bytes, bytes]) -> None:
    for line in record:
        handle.write(line)


def open_fastq_binary(path: Path):  # type: ignore[no-untyped-def]
    if path.name.lower().endswith((".gz", ".bgz", ".bgzf")):
        return gzip.open(path, "rb")
    return path.open("rb")


def consume_chunk(
    r1_handle: BinaryIO,
    r2_handle: BinaryIO | None,
    count: int,
    r1_out: BinaryIO | None,
    r2_out: BinaryIO | None,
) -> int:
    consumed = 0
    for _ in range(count):
        record1 = read_fastq_record(r1_handle)
        if record1 is None:
            break
        record2 = read_fastq_record(r2_handle) if r2_handle is not None else None
        if r2_handle is not None and record2 is None:
            raise RuntimeError("R2 ended before R1 while creating paired alignment chunks")
        if r1_out is not None:
            write_record(r1_out, record1)
        if r2_out is not None and record2 is not None:
            write_record(r2_out, record2)
        consumed += 1
    return consumed



def stage_single_end_chunk(
    r1_handle: BinaryIO,
    count: int,
    output_path: Path,
) -> int:
    """Materialize one single-end FASTQ chunk atomically for BWA.

    The next chunk can be prepared by a background thread while BWA/samtools
    process the current chunk.  At most two uncompressed chunks are present
    transiently (current + prefetched next), so restartable checkpoints remain
    small while gzip decompression/I/O is hidden behind alignment.
    """
    part = output_path.with_name(output_path.name + ".part")
    for path in (part, output_path):
        try:
            path.unlink()
        except OSError:
            pass
    written = 0
    try:
        with part.open("wb", buffering=STREAM_WRITE_BUFFER) as out:
            written = consume_chunk(r1_handle, None, count, out, None)
        if written != count:
            raise RuntimeError(
                f"FASTQ ended while staging chunk {output_path.name}: "
                f"expected {count:,}, wrote {written:,}"
            )
        os.replace(part, output_path)
        return written
    except BaseException:
        try:
            part.unlink()
        except OSError:
            pass
        raise


def prefetch_single_end_chunk(
    r1_handle: BinaryIO,
    count: int,
    output_path: Path,
    result: dict,
) -> None:
    """Background wrapper for one-record-stream-ahead FASTQ staging."""
    try:
        result["written"] = stage_single_end_chunk(r1_handle, count, output_path)
    except BaseException as exc:
        result["error"] = exc


def clean_chunk_attempt_files(
    checkpoint_dir: Path,
    chunk_index: int,
    output_bam: Path,
) -> None:
    """Remove only failed-attempt BAM/sort scratch, never staged FASTQ/checkpoints."""
    part = output_bam.with_name(output_bam.stem + ".part.bam")
    try:
        part.unlink()
    except OSError:
        pass
    prefix = f"{chunk_index:05d}"
    tmp_dir = checkpoint_dir / "tmp"
    if tmp_dir.is_dir():
        for pattern in (f"coord_{prefix}*", f"collate_{prefix}*"):
            for path in tmp_dir.glob(pattern):
                if path.is_file():
                    try:
                        path.unlink()
                    except OSError:
                        pass


def checkpoint_bytes(chunks: list[dict]) -> int:
    return sum(max(0, int(item.get("bam_length", 0))) for item in chunks)


def disk_guard(
    checkpoint_dir: Path,
    chunks: list[dict],
    completed_sequences: int,
    total_sequences: int,
    reserve_gib: float,
) -> None:
    free = shutil.disk_usage(checkpoint_dir).free
    reserve = int(reserve_gib * GIB)
    stored = checkpoint_bytes(chunks)
    remaining_sequences = max(0, total_sequences - completed_sequences)
    predicted_remaining = 0
    transient = 8 * GIB
    if completed_sequences > 0 and stored > 0:
        bytes_per_sequence = stored / completed_sequences
        predicted_remaining = int(bytes_per_sequence * remaining_sequences)
        average_chunk = stored / max(1, len(chunks))
        transient = int(max(8 * GIB, average_chunk * 3.0))
    required = reserve + predicted_remaining + transient
    detail = (
        f"DISK {format_gib(free)} free | checkpoints {format_gib(stored)} | "
        f"reserve {reserve_gib:.0f} GiB"
    )
    if predicted_remaining:
        detail += f" | projected remaining checkpoints {format_gib(predicted_remaining)}"
    stage15_progress(8.0 + 82.0 * (completed_sequences / max(1, total_sequences)), detail)
    minimum_start = reserve + transient
    if free < minimum_start:
        raise SystemExit(
            "SAFE STOP: insufficient free disk for the next resumable alignment chunk. "
            f"Free {format_gib(free)}, safety requirement {format_gib(minimum_start)}. "
            "Completed checkpoint BAMs were preserved; free space and rerun item 15 to continue."
        )
    if predicted_remaining and free < required:
        raise SystemExit(
            "SAFE STOP: current checkpoint growth predicts that the volume will run out before Stage 15 can finish. "
            f"Free {format_gib(free)}, projected need including safety reserve {format_gib(required)}. "
            "Completed checkpoint BAMs were preserved; free/move storage and rerun item 15 to continue."
        )


def quickcheck_bam(root: Path, samtools_image: str, bam: Path) -> None:
    run(
        [
            "docker",
            "run",
            "--rm",
            *mount(root, "/work"),
            "--entrypoint",
            "samtools",
            samtools_image,
            "quickcheck",
            "-v",
            "/work/" + bam.relative_to(root).as_posix(),
        ]
    )


def quickcheck_checkpoint(checkpoint_dir: Path, samtools_image: str, bam: Path) -> None:
    run(
        [
            "docker",
            "run",
            "--rm",
            *mount(checkpoint_dir, "/checkpoint"),
            "--entrypoint",
            "samtools",
            samtools_image,
            "quickcheck",
            "-v",
            f"/checkpoint/{bam.name}",
        ]
    )


def process_single_end_stream_chunk(
    *,
    r1_handle: BinaryIO,
    expected_records: int,
    checkpoint_dir: Path,
    index_dir: Path,
    bwa_image: str,
    samtools_image: str,
    bwa_threads: int,
    sort_threads: int,
    sort_memory: str,
    read_group: str,
    chunk_index: int,
    output_bam: Path,
    total_sequences: int,
    completed_sequences: int,
) -> int:
    """Stream one single-end FASTQ chunk directly from gzip reader into BWA stdin."""
    part = output_bam.with_name(output_bam.stem + ".part.bam")
    try:
        part.unlink()
    except OSError:
        pass
    tmp_dir = checkpoint_dir / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{chunk_index:05d}"

    run_token = f"{os.getpid()}-{chunk_index:05d}"
    bwa_container = f"dante-stage15-bwa-{run_token}"
    samtools_container = f"dante-stage15-sort-{run_token}"
    bwa_command = build_bwa_stdin_command(index_dir, bwa_image, bwa_threads, read_group, bwa_container)
    samtools_pipeline = (
        "set -o pipefail; "
        f"samtools sort -@ {sort_threads} -m {sort_memory} "
        f"-T /checkpoint/tmp/coord_{prefix} -o /checkpoint/{part.name} -"
    )
    samtools_command = [
        "docker",
        "run",
        "--rm",
        *managed_container_args(samtools_container),
        "-i",
        *mount(checkpoint_dir, "/checkpoint"),
        "--entrypoint",
        "bash",
        samtools_image,
        "-lc",
        samtools_pipeline,
    ]

    progress = AlignmentProgress(total_sequences, completed_sequences)
    bwa = subprocess.Popen(
        bwa_command,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert bwa.stdin is not None and bwa.stdout is not None and bwa.stderr is not None
    stderr_thread = threading.Thread(target=pump_bwa, args=(bwa.stderr, progress), daemon=True)
    stderr_thread.start()
    samtools = subprocess.Popen(samtools_command, stdin=bwa.stdout)
    bwa.stdout.close()

    written = 0
    write_error: BaseException | None = None
    buffer = bytearray()
    try:
        while written < expected_records:
            record = read_fastq_record(r1_handle)
            if record is None:
                break
            buffer.extend(record[0])
            buffer.extend(record[1])
            buffer.extend(record[2])
            buffer.extend(record[3])
            written += 1
            if len(buffer) >= STREAM_WRITE_BUFFER:
                bwa.stdin.write(buffer)
                buffer.clear()
        if buffer:
            bwa.stdin.write(buffer)
            buffer.clear()
    except (BrokenPipeError, OSError) as exc:
        write_error = exc
    finally:
        try:
            bwa.stdin.close()
        except OSError:
            pass

    try:
        bwa_code = bwa.wait()
        stderr_thread.join()
    except BaseException:
        force_remove_container(bwa_container)
        force_remove_container(samtools_container)
        raise
    if bwa_code:
        try:
            samtools.terminate()
        except OSError:
            pass
        samtools_code = samtools.wait()
        force_remove_container(bwa_container)
        force_remove_container(samtools_container)
        message = (
            f"alignment chunk failed while streaming FASTQ (bwa={bwa_code}, samtools={samtools_code}, "
            f"records_sent={written:,}, writer_error={write_error!r})"
        )
        if bwa_code == 125 or samtools_code == 125:
            raise TransientDockerError(message)
        raise RuntimeError(message)
    try:
        samtools_code = samtools.wait()
    except BaseException:
        force_remove_container(bwa_container)
        force_remove_container(samtools_container)
        raise
    if write_error is not None:
        force_remove_container(bwa_container)
        force_remove_container(samtools_container)
        if isinstance(write_error, BrokenPipeError):
            raise TransientDockerError(f"FASTQ streaming pipe broke after {written:,} records: {write_error}")
        raise RuntimeError(f"FASTQ streaming failed after {written:,} records: {write_error}")
    if samtools_code:
        force_remove_container(bwa_container)
        force_remove_container(samtools_container)
        message = f"alignment chunk failed (bwa={bwa_code}, samtools={samtools_code})"
        if samtools_code == 125:
            raise TransientDockerError(message)
        raise RuntimeError(message)
    if written != expected_records:
        raise RuntimeError(
            f"FASTQ ended unexpectedly in chunk {chunk_index}: expected {expected_records:,}, streamed {written:,}"
        )
    if not part.is_file() or part.stat().st_size <= 0:
        raise RuntimeError(f"alignment chunk {chunk_index} did not create a non-empty BAM")
    quickcheck_checkpoint(checkpoint_dir, samtools_image, part)
    os.replace(part, output_bam)
    return written


def process_chunk(
    *,
    checkpoint_dir: Path,
    index_dir: Path,
    bwa_image: str,
    samtools_image: str,
    bwa_threads: int,
    sort_threads: int,
    sort_memory: str,
    read_group: str,
    chunk_index: int,
    r1_chunk: Path,
    r2_chunk: Path | None,
    output_bam: Path,
    total_sequences: int,
    completed_sequences: int,
) -> None:
    part = output_bam.with_name(output_bam.stem + ".part.bam")
    try:
        part.unlink()
    except OSError:
        pass
    tmp_dir = checkpoint_dir / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{chunk_index:05d}"
    r2_name = r2_chunk.name if r2_chunk is not None else None
    run_token = f"{os.getpid()}-{chunk_index:05d}"
    bwa_container = f"dante-stage15-bwa-{run_token}"
    samtools_container = f"dante-stage15-sort-{run_token}"
    bwa_command = build_bwa_command(
        index_dir,
        checkpoint_dir,
        bwa_image,
        bwa_threads,
        read_group,
        r1_chunk.name,
        r2_name,
        bwa_container,
    )

    if r2_chunk is None:
        # Single-end duplicate matching only needs coordinate, position and orientation.
        # No whole-genome name sort/collate is required.
        samtools_pipeline = (
            "set -o pipefail; "
            f"samtools sort -@ {sort_threads} -m {sort_memory} "
            f"-T /checkpoint/tmp/coord_{prefix} -o /checkpoint/{part.name} -"
        )
    else:
        # Paired-end fixmate requires mates to be grouped. collate is explicitly designed
        # for this and avoids the expensive full query-name sort used in older releases.
        samtools_pipeline = (
            "set -o pipefail; "
            f"samtools collate -O -u -@ {sort_threads} -T /checkpoint/tmp/collate_{prefix} - | "
            "samtools fixmate -m -u - - | "
            f"samtools sort -@ {sort_threads} -m {sort_memory} "
            f"-T /checkpoint/tmp/coord_{prefix} -o /checkpoint/{part.name} -"
        )
    samtools_command = [
        "docker",
        "run",
        "--rm",
        *managed_container_args(samtools_container),
        "-i",
        *mount(checkpoint_dir, "/checkpoint"),
        "--entrypoint",
        "bash",
        samtools_image,
        "-lc",
        samtools_pipeline,
    ]

    progress = AlignmentProgress(total_sequences, completed_sequences)
    bwa = subprocess.Popen(bwa_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert bwa.stdout is not None and bwa.stderr is not None
    stderr_thread = threading.Thread(target=pump_bwa, args=(bwa.stderr, progress), daemon=True)
    stderr_thread.start()
    samtools = subprocess.Popen(samtools_command, stdin=bwa.stdout)
    bwa.stdout.close()
    try:
        bwa_code = bwa.wait()
        stderr_thread.join()
    except BaseException:
        force_remove_container(bwa_container)
        force_remove_container(samtools_container)
        raise
    if bwa_code:
        try:
            samtools.terminate()
        except OSError:
            pass
        samtools_code = samtools.wait()
        force_remove_container(bwa_container)
        force_remove_container(samtools_container)
        message = f"alignment chunk failed (bwa={bwa_code}, samtools={samtools_code})"
        if bwa_code == 125 or samtools_code == 125:
            raise TransientDockerError(message)
        raise RuntimeError(message)
    try:
        samtools_code = samtools.wait()
    except BaseException:
        force_remove_container(bwa_container)
        force_remove_container(samtools_container)
        raise
    if samtools_code:
        force_remove_container(bwa_container)
        force_remove_container(samtools_container)
        message = f"alignment chunk failed (bwa={bwa_code}, samtools={samtools_code})"
        if samtools_code == 125:
            raise TransientDockerError(message)
        raise RuntimeError(message)
    if not part.is_file() or part.stat().st_size <= 0:
        raise RuntimeError(f"alignment chunk {chunk_index} did not create a non-empty BAM")
    quickcheck_checkpoint(checkpoint_dir, samtools_image, part)
    os.replace(part, output_bam)


def finalize_chunks(
    *,
    root: Path,
    checkpoint_dir: Path,
    data_dir: Path,
    fasta_name: str,
    targets: Path,
    output: Path,
    samtools_image: str,
    chunks: list[dict],
    merge_threads: int,
    index_threads: int,
) -> tuple[Path, int]:
    if not chunks:
        raise RuntimeError("no completed alignment chunks are available for finalization")
    part = output.with_name(output.stem + ".part" + output.suffix)
    part_index = Path(str(part) + ".bai")
    for path in (part, part_index):
        try:
            path.unlink()
        except OSError:
            pass
    targets_container = "/work/" + targets.relative_to(root).as_posix()
    part_container = "/work/" + part.relative_to(root).as_posix()
    chunk_names = [f"/checkpoint/{item['bam_name']}" for item in chunks]
    # Coordinate-sorted chunk BAMs contain the complete aligned genome. Merge preserves
    # coordinate order; one global markdup pass therefore sees duplicate candidates across
    # chunk boundaries. Only after that do we retain target-overlapping evidence.
    pipeline = (
        "set -o pipefail; "
        f"samtools merge -u -c -p -@ {merge_threads} -o - {' '.join(chunk_names)} | "
        f"samtools markdup -@ {merge_threads} -u - - | "
        f"samtools view -@ {merge_threads} -b -L '{targets_container}' -o '{part_container}' -"
    )
    run(
        [
            "docker",
            "run",
            "--rm",
            *mount(root, "/work"),
            *mount(checkpoint_dir, "/checkpoint"),
            *mount(data_dir, "/data"),
            "--entrypoint",
            "bash",
            samtools_image,
            "-lc",
            pipeline,
        ]
    )
    if not part.is_file() or part.stat().st_size <= 0:
        raise RuntimeError("final merge/markdup did not create a non-empty partial target BAM")
    quickcheck_bam(root, samtools_image, part)
    alignment_count = capture(
        [
            "docker",
            "run",
            "--rm",
            *mount(root, "/work"),
            "--entrypoint",
            "samtools",
            samtools_image,
            "view",
            "-c",
            part_container,
        ]
    )
    if not alignment_count.isdigit() or int(alignment_count) <= 0:
        raise RuntimeError("target BAM contains zero alignments; check FASTQ, GRCh37 reference and target panel")
    stage15_progress(97.0, f"{int(alignment_count):,} target alignments retained | build BAM index")
    run(
        [
            "docker",
            "run",
            "--rm",
            *mount(root, "/work"),
            "--entrypoint",
            "samtools",
            samtools_image,
            "index",
            "-@",
            str(index_threads),
            part_container,
        ]
    )
    if not part_index.is_file() or part_index.stat().st_size <= 0:
        raise RuntimeError("samtools did not create the BAM index")
    return part, int(alignment_count)


def remove_legacy_unresumable_sort_files(root: Path) -> None:
    legacy = root / "Dante_work" / "reference_sort_tmp"
    if not legacy.is_dir():
        return
    removed = 0
    for pattern in ("name*", "coordinate*"):
        for path in legacy.glob(pattern):
            if path.is_file():
                try:
                    path.unlink()
                    removed += 1
                except OSError:
                    pass
    if removed:
        print(
            f"[CHECKPOINT] Removed {removed} legacy v7.5.10 unresumable sort scratch files.",
            file=sys.stderr,
            flush=True,
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--index-dir", type=Path, required=True)
    parser.add_argument("--fasta-name", required=True)
    parser.add_argument("--fasta-sha256", required=True)
    parser.add_argument("--r1", type=Path, required=True)
    parser.add_argument("--r1-records", type=int, required=True)
    parser.add_argument("--r1-sha256", required=True)
    parser.add_argument("--r2", type=Path)
    parser.add_argument("--r2-records", type=int)
    parser.add_argument("--r2-sha256")
    parser.add_argument("--targets-bed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bwa-image", required=True)
    parser.add_argument("--bwa-contract", required=True)
    parser.add_argument("--samtools-image", required=True)
    parser.add_argument("--samtools-contract", required=True)
    parser.add_argument("--threads", type=int, default=4, help="legacy/general thread count")
    parser.add_argument("--bwa-threads", type=int)
    parser.add_argument("--sort-threads", type=int)
    parser.add_argument("--sort-memory", default="512M")
    parser.add_argument("--index-threads", type=int)
    parser.add_argument("--chunk-records", type=int, default=DEFAULT_CHUNK_RECORDS)
    parser.add_argument("--checkpoint-dir", type=Path)
    parser.add_argument("--disk-reserve-gb", type=float, default=15.0)
    args = parser.parse_args()

    root = args.root.resolve()
    data_dir = args.data_dir.resolve()
    index_dir = args.index_dir.resolve()
    r1 = args.r1.resolve()
    r2 = args.r2.resolve() if args.r2 else None
    targets = args.targets_bed.resolve()
    output = args.output.resolve()
    checkpoint_dir = (args.checkpoint_dir or (root / "Dante_work" / "reference_alignment_checkpoints")).resolve()
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    threads = max(1, args.threads)
    bwa_threads = max(1, args.bwa_threads or threads)
    sort_threads = max(1, args.sort_threads or min(threads, 4))
    index_threads = max(1, args.index_threads or threads)
    sort_memory = args.sort_memory.strip() or "512M"
    if not re.fullmatch(r"[1-9][0-9]*(?:[KMG])?", sort_memory, re.IGNORECASE):
        raise SystemExit(f"invalid --sort-memory value: {sort_memory}")
    if args.chunk_records <= 0:
        raise SystemExit("--chunk-records must be positive")
    if args.r1_records <= 0:
        raise SystemExit("--r1-records must be a positive exact FASTQ record count")
    if not re.fullmatch(r"[0-9a-fA-F]{64}", args.r1_sha256):
        raise SystemExit("--r1-sha256 must be an exact SHA-256")
    if r2 is None and args.r2_records not in (None, 0):
        raise SystemExit("--r2-records was supplied without --r2")
    if r2 is None and args.r2_sha256:
        raise SystemExit("--r2-sha256 was supplied without --r2")
    if r2 is not None and (args.r2_records is None or args.r2_records <= 0):
        raise SystemExit("paired-end input requires a positive --r2-records count")
    if r2 is not None and args.r1_records != args.r2_records:
        raise SystemExit(f"paired FASTQ record counts differ: R1={args.r1_records:,}, R2={args.r2_records:,}")
    if r2 is not None and not re.fullmatch(r"[0-9a-fA-F]{64}", args.r2_sha256 or ""):
        raise SystemExit("paired-end input requires an exact --r2-sha256")
    total_sequences = args.r1_records + (args.r2_records or 0)
    read_mode = "paired" if r2 is not None else "single"
    required_files = [r1, targets, data_dir / "fasta" / args.fasta_name]
    if r2 is not None:
        required_files.append(r2)
    for path in required_files:
        if not path.is_file():
            raise SystemExit(f"Required file is missing: {path}")
    for path in [targets, output.parent]:
        try:
            path.relative_to(root)
        except ValueError as exc:
            raise SystemExit(f"Path must be inside pipeline root: {path}") from exc

    cleanup_stale_stage15_containers()
    stage15_progress(6.0, "REFERENCE PREP | verify/build persistent BWA GRCh37 index")
    ensure_bwa_index(data_dir, index_dir, args.fasta_name, args.fasta_sha256, args.bwa_image, args.bwa_contract)
    stage15_progress(
        8.0,
        f"REFERENCE PREP complete | {format_reads(total_sequences)} input reads | "
        f"BWA {bwa_threads}t | sort {sort_threads}t x {sort_memory} | resumable {args.chunk_records:,}-record chunks",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    remove_legacy_unresumable_sort_files(root)
    clean_transient_files(checkpoint_dir)

    sample = re_safe(r1.stem.split("_R1")[0] or "Dante")
    read_group = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"
    identity = checkpoint_identity(args, read_mode)
    manifest_path = checkpoint_dir / "checkpoint_manifest.json"
    manifest = load_manifest(manifest_path)
    if manifest and not manifest_matches(manifest, identity):
        reset_stale_checkpoints(checkpoint_dir, manifest_path)
        manifest = {}
    if not manifest:
        manifest = {
            "schema": CHECKPOINT_SCHEMA,
            "status": "in_progress",
            "identity": identity,
            "completed_chunks": [],
        }
        atomic_json(manifest_path, manifest)

    chunks = valid_completed_chunks(manifest, checkpoint_dir)
    if len(chunks) != len(manifest.get("completed_chunks", [])):
        manifest["completed_chunks"] = chunks
        atomic_json(manifest_path, manifest)
    completed_mate_records = sum(int(item["records_per_mate"]) for item in chunks)
    completed_sequences = completed_mate_records * (2 if r2 is not None else 1)
    if chunks:
        stage15_progress(
            8.0 + 82.0 * (completed_sequences / total_sequences),
            f"RESUME | {len(chunks)} completed chunks | {format_reads(completed_sequences)} / {format_reads(total_sequences)} reads preserved",
        )

    if completed_mate_records < args.r1_records:
        with open_fastq_binary(r1) as r1_handle:
            r2_context = open_fastq_binary(r2) if r2 is not None else None
            try:
                r2_handle = r2_context
                consumed_records = 0
                chunk_index = 0

                # Scan exactly across already-completed checkpoint boundaries.
                # Each manifest entry is authoritative for its own record count, so
                # older 5M/20M checkpoints remain compatible with today's 10M policy.
                resume_scan_total = completed_mate_records
                resume_scan_started = time.monotonic()
                for stored_chunk in chunks:
                    stored_records = int(stored_chunk["records_per_mate"])
                    skipped = consume_chunk(r1_handle, r2_handle, stored_records, None, None)
                    if skipped != stored_records:
                        raise RuntimeError("FASTQ ended while scanning to a completed checkpoint boundary")
                    consumed_records += skipped
                    chunk_index += 1
                    elapsed_scan = max(0.001, time.monotonic() - resume_scan_started)
                    scan_speed = consumed_records / elapsed_scan
                    scan_remaining = max(0, resume_scan_total - consumed_records)
                    scan_eta = scan_remaining / scan_speed if scan_speed > 0 else None
                    stage15_progress(
                        8.0 + 82.0 * (completed_sequences / total_sequences),
                        f"RESUME SCAN | {format_reads(consumed_records)} / "
                        f"{format_reads(resume_scan_total)} previously completed records skipped | "
                        f"{format_reads(int(scan_speed))} rec/s | ETA {format_duration(scan_eta)}",
                    )
                if resume_scan_total:
                    stage15_progress(
                        8.0 + 82.0 * (completed_sequences / total_sequences),
                        f"RESUME SCAN complete | next alignment starts at "
                        f"{format_reads(completed_sequences)} / {format_reads(total_sequences)} reads",
                    )

                if r2 is None:
                    # Double-buffered single-end path:
                    #   gzip R1 -> staged current FASTQ -> BWA/sort
                    #                  while
                    #              background thread prepares next FASTQ
                    #
                    # This intentionally returns to file-backed BWA input because the
                    # direct Python->Docker stdin stream showed repeated CPU starvation.
                    # Only one chunk is prefetched, so disk overhead is bounded.
                    initial_records = min(args.chunk_records, args.r1_records - consumed_records)
                    current_chunk = checkpoint_dir / f"input_chunk_{chunk_index:05d}_R1.fastq"
                    stage15_progress(
                        8.0 + 82.0 * (completed_sequences / total_sequences),
                        f"DOUBLE BUFFER | stage first {format_reads(initial_records)}-read chunk",
                    )
                    current_written = stage_single_end_chunk(r1_handle, initial_records, current_chunk)

                    while consumed_records < args.r1_records:
                        expected_records = current_written
                        output_chunk = checkpoint_dir / f"chunk_{chunk_index:05d}.bam"
                        disk_guard(
                            checkpoint_dir,
                            chunks,
                            completed_sequences,
                            total_sequences,
                            args.disk_reserve_gb,
                        )
                        stage15_progress(
                            8.0 + 82.0 * (completed_sequences / total_sequences),
                            f"CHUNK {chunk_index + 1} | align {format_reads(expected_records)} reads | "
                            "double-buffered checkpoint",
                        )

                        # Start preparing the following chunk before BWA starts on the
                        # current staged file.  The gzip reader is touched only by this
                        # one background producer, never by BWA/samtools.
                        remaining_after_current = args.r1_records - (consumed_records + expected_records)
                        prefetch_thread = None
                        prefetch_result: dict = {}
                        next_chunk = None
                        next_expected = 0
                        if remaining_after_current > 0:
                            next_expected = min(args.chunk_records, remaining_after_current)
                            next_chunk = checkpoint_dir / f"input_chunk_{chunk_index + 1:05d}_R1.fastq"
                            prefetch_thread = threading.Thread(
                                target=prefetch_single_end_chunk,
                                args=(r1_handle, next_expected, next_chunk, prefetch_result),
                                daemon=True,
                            )
                            prefetch_thread.start()

                        max_chunk_attempts = 3
                        terminal_error: BaseException | None = None
                        for chunk_attempt in range(1, max_chunk_attempts + 1):
                            try:
                                process_chunk(
                                    checkpoint_dir=checkpoint_dir,
                                    index_dir=index_dir,
                                    bwa_image=args.bwa_image,
                                    samtools_image=args.samtools_image,
                                    bwa_threads=bwa_threads,
                                    sort_threads=sort_threads,
                                    sort_memory=sort_memory,
                                    read_group=read_group,
                                    chunk_index=chunk_index,
                                    r1_chunk=current_chunk,
                                    r2_chunk=None,
                                    output_bam=output_chunk,
                                    total_sequences=total_sequences,
                                    completed_sequences=completed_sequences,
                                )
                                break
                            except TransientDockerError as exc:
                                clean_chunk_attempt_files(checkpoint_dir, chunk_index, output_chunk)
                                if chunk_attempt >= max_chunk_attempts:
                                    terminal_error = RuntimeError(
                                        f"Docker transport failed repeatedly while processing chunk {chunk_index + 1}; "
                                        "completed checkpoints are preserved. Restart Docker Desktop if needed and rerun item 15."
                                    )
                                    terminal_error.__cause__ = exc
                                    break
                                stage15_progress(
                                    8.0 + 82.0 * (completed_sequences / total_sequences),
                                    f"DOCKER RECOVERY | transient engine failure in chunk {chunk_index + 1} | "
                                    f"retry {chunk_attempt + 1}/{max_chunk_attempts}; staged FASTQ + completed checkpoints preserved",
                                )
                                cleanup_stale_stage15_containers()
                                if not wait_for_docker_engine(120):
                                    terminal_error = RuntimeError(
                                        "Docker engine did not recover within 120 seconds; completed checkpoints are preserved. "
                                        "Restart Docker Desktop/WSL and rerun item 15."
                                    )
                                    terminal_error.__cause__ = exc
                                    break
                                stage15_progress(
                                    8.0 + 82.0 * (completed_sequences / total_sequences),
                                    f"DOCKER RECOVERY complete | retry chunk {chunk_index + 1} from staged FASTQ",
                                )
                            except BaseException as exc:
                                terminal_error = exc
                                break

                        # Never close/advance the gzip context while its prefetch thread
                        # is still using it.
                        if prefetch_thread is not None:
                            prefetch_thread.join()

                        if terminal_error is not None:
                            try:
                                current_chunk.unlink()
                            except OSError:
                                pass
                            raise terminal_error

                        written = expected_records
                        consumed_records += written
                        chunk_entry = {
                            "index": chunk_index,
                            "records_per_mate": written,
                            "input_sequences": written,
                            "bam_name": output_chunk.name,
                            "bam_length": output_chunk.stat().st_size,
                        }
                        chunks.append(chunk_entry)
                        completed_sequences += written
                        manifest["completed_chunks"] = chunks
                        manifest["status"] = "in_progress"
                        manifest["current_chunk_records"] = int(args.chunk_records)
                        atomic_json(manifest_path, manifest)
                        stage15_progress(
                            8.0 + 82.0 * (completed_sequences / total_sequences),
                            f"CHECKPOINT SAVED | chunk {chunk_index + 1} | "
                            f"{format_reads(completed_sequences)} / {format_reads(total_sequences)} reads",
                        )

                        try:
                            current_chunk.unlink()
                        except OSError:
                            pass

                        chunk_index += 1
                        if next_chunk is None:
                            break
                        if "error" in prefetch_result:
                            raise RuntimeError(
                                f"Next FASTQ chunk prefetch failed after checkpoint {chunk_index}: "
                                f"{prefetch_result['error']}"
                            ) from prefetch_result["error"]
                        prefetched = int(prefetch_result.get("written", 0))
                        if prefetched != next_expected:
                            raise RuntimeError(
                                f"Prefetched FASTQ chunk size mismatch: expected {next_expected:,}, got {prefetched:,}"
                            )
                        current_chunk = next_chunk
                        current_written = prefetched

                else:
                    # Paired-end keeps file-backed R1/R2 staging. Add the same Docker
                    # transient retry semantics used by the single-end path.
                    while consumed_records < args.r1_records:
                        expected_records = min(args.chunk_records, args.r1_records - consumed_records)
                        disk_guard(
                            checkpoint_dir,
                            chunks,
                            completed_sequences,
                            total_sequences,
                            args.disk_reserve_gb,
                        )
                        output_chunk = checkpoint_dir / f"chunk_{chunk_index:05d}.bam"
                        chunk_sequences = expected_records * 2
                        stage15_progress(
                            8.0 + 82.0 * (completed_sequences / total_sequences),
                            f"CHUNK {chunk_index + 1} | align {format_reads(chunk_sequences)} reads | checkpoint after completion",
                        )
                        r1_chunk = checkpoint_dir / f"input_chunk_{chunk_index:05d}_R1.fastq"
                        r2_chunk = checkpoint_dir / f"input_chunk_{chunk_index:05d}_R2.fastq"
                        with r1_chunk.open("wb") as out1, r2_chunk.open("wb") as out2:
                            written = consume_chunk(r1_handle, r2_handle, expected_records, out1, out2)
                        if written != expected_records:
                            raise RuntimeError(
                                f"FASTQ ended unexpectedly in chunk {chunk_index}: expected {expected_records:,}, wrote {written:,}"
                            )

                        max_chunk_attempts = 3
                        try:
                            for chunk_attempt in range(1, max_chunk_attempts + 1):
                                try:
                                    process_chunk(
                                        checkpoint_dir=checkpoint_dir,
                                        index_dir=index_dir,
                                        bwa_image=args.bwa_image,
                                        samtools_image=args.samtools_image,
                                        bwa_threads=bwa_threads,
                                        sort_threads=sort_threads,
                                        sort_memory=sort_memory,
                                        read_group=read_group,
                                        chunk_index=chunk_index,
                                        r1_chunk=r1_chunk,
                                        r2_chunk=r2_chunk,
                                        output_bam=output_chunk,
                                        total_sequences=total_sequences,
                                        completed_sequences=completed_sequences,
                                    )
                                    break
                                except TransientDockerError as exc:
                                    clean_chunk_attempt_files(checkpoint_dir, chunk_index, output_chunk)
                                    if chunk_attempt >= max_chunk_attempts:
                                        raise RuntimeError(
                                            f"Docker transport failed repeatedly while processing paired chunk {chunk_index + 1}; "
                                            "completed checkpoints are preserved."
                                        ) from exc
                                    cleanup_stale_stage15_containers()
                                    if not wait_for_docker_engine(120):
                                        raise RuntimeError(
                                            "Docker engine did not recover within 120 seconds; completed checkpoints are preserved."
                                        ) from exc
                        finally:
                            for path in (r1_chunk, r2_chunk):
                                try:
                                    path.unlink()
                                except OSError:
                                    pass
                            clean_transient_files(checkpoint_dir)

                        consumed_records += written
                        chunk_entry = {
                            "index": chunk_index,
                            "records_per_mate": written,
                            "input_sequences": written * 2,
                            "bam_name": output_chunk.name,
                            "bam_length": output_chunk.stat().st_size,
                        }
                        chunks.append(chunk_entry)
                        completed_sequences += written * 2
                        manifest["completed_chunks"] = chunks
                        manifest["status"] = "in_progress"
                        manifest["current_chunk_records"] = int(args.chunk_records)
                        atomic_json(manifest_path, manifest)
                        stage15_progress(
                            8.0 + 82.0 * (completed_sequences / total_sequences),
                            f"CHECKPOINT SAVED | chunk {chunk_index + 1} | "
                            f"{format_reads(completed_sequences)} / {format_reads(total_sequences)} reads",
                        )
                        chunk_index += 1
            finally:
                if r2_context is not None:
                    r2_context.close()

    if completed_sequences != total_sequences:
        raise RuntimeError(
            f"checkpoint accounting mismatch: completed={completed_sequences:,}, expected={total_sequences:,}"
        )
    stage15_progress(90.0, f"ALL ALIGNMENT CHUNKS COMPLETE | {len(chunks)} checkpoints | global merge + markdup")
    manifest["status"] = "alignment_complete"
    atomic_json(manifest_path, manifest)

    part, alignment_count = finalize_chunks(
        root=root,
        checkpoint_dir=checkpoint_dir,
        data_dir=data_dir,
        fasta_name=args.fasta_name,
        targets=targets,
        output=output,
        samtools_image=args.samtools_image,
        chunks=chunks,
        merge_threads=sort_threads,
        index_threads=index_threads,
    )
    part_index = Path(str(part) + ".bai")
    stage15_progress(99.0, "BAM index complete | atomically publish finalized alignment")
    os.replace(part, output)
    os.replace(part_index, Path(str(output) + ".bai"))
    manifest["status"] = "finalized_pending_orchestrator"
    manifest["target_alignment_count"] = alignment_count
    atomic_json(manifest_path, manifest)
    # Keep checkpoints until the PowerShell orchestrator has written and validated the
    # final BAM state file.  This closes the small crash window between Python publishing
    # the BAM and the orchestrator publishing provenance/state.
    stage15_progress(100.0, "DONE | target-only BAM and BAI finalized | checkpoints retained until state validation")
    return 0


def re_safe(value: str) -> str:
    return "".join(character if character.isalnum() or character in "_.-" else "_" for character in value)[:80]


if __name__ == "__main__":
    raise SystemExit(main())
