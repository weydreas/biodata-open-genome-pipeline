#!/usr/bin/env python3
"""Run a subprocess while emitting a durable activity/progress heartbeat.

Most Dante commands keep their normal console streams and receive a time-based
heartbeat.  AnnotSV can additionally be observed through its stable console
milestones and the read position of its internal ``*.part.tsv.tmp`` stream.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import threading
import time
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a command with periodic activity output")
    parser.add_argument("--label", required=True)
    parser.add_argument("--heartbeat-seconds", type=float, default=15.0)
    parser.add_argument("--watch-file", action="append", type=Path, default=[])
    parser.add_argument("--watch-dir", action="append", type=Path, default=[])
    parser.add_argument("--progress-profile", choices=("annotsv",))
    parser.add_argument("--container-name")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command and args.command[0] == "--":
        args.command = args.command[1:]
    if not args.command:
        parser.error("a command is required after --")
    if args.heartbeat_seconds < 1:
        parser.error("--heartbeat-seconds must be at least 1")
    if args.progress_profile == "annotsv" and not args.container_name:
        parser.error("--container-name is required for the annotsv progress profile")
    return args


def format_bytes(value: int) -> str:
    units = ("bytes", "KiB", "MiB", "GiB", "TiB")
    amount = float(value)
    for unit in units:
        if amount < 1024 or unit == units[-1]:
            if unit == "bytes":
                return f"{int(amount):,} {unit}"
            return f"{amount:,.2f} {unit}"
        amount /= 1024
    return f"{value:,} bytes"


def format_elapsed(seconds: float) -> str:
    total = max(0, int(seconds))
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def directory_snapshot(path: Path) -> tuple[int, int] | None:
    if not path.is_dir():
        return None
    files = 0
    size = 0
    try:
        for root, _dirs, names in os.walk(path):
            for name in names:
                files += 1
                try:
                    size += (Path(root) / name).stat().st_size
                except OSError:
                    pass
    except OSError:
        return None
    return files, size


def watch_text(files: list[Path], directories: list[Path]) -> str:
    details: list[str] = []
    for path in files:
        try:
            if path.is_file():
                details.append(f"{path.name}: {format_bytes(path.stat().st_size)}")
            else:
                details.append(f"{path.name}: waiting for output")
        except OSError:
            details.append(f"{path.name}: temporarily unavailable")
    for path in directories:
        snapshot = directory_snapshot(path)
        if snapshot is None:
            details.append(f"{path.name}: waiting for output")
        else:
            count, size = snapshot
            details.append(f"{path.name}: {count:,} files, {format_bytes(size)}")
    return " | " + " | ".join(details) if details else ""


class AnnotSVProgress:
    """Translate AnnotSV 3.5.x milestones and its active input offset to progress."""

    _MILESTONES = (
        ("downloading the configuration data", 3.0, "load configuration"),
        ("checking all these configuration data", 7.0, "validate configuration"),
        ("vcf to bed", 12.0, "convert VCF to BED"),
        ("checking the annotation data sources", 18.0, "check annotation sources"),
        ("listing arguments", 23.0, "prepare run"),
        ("searching for sv overlaps", 28.0, "find gene/regulatory overlaps"),
        ("listing of the annotations to be realized", 32.0, "build annotation plan"),
        ("annotation in progress", 35.0, "main annotation pass"),
        ("writing of ", 92.0, "write final table"),
    )

    def __init__(self, watch_files: list[Path], container_name: str) -> None:
        self._lock = threading.Lock()
        self._watch_file = watch_files[0] if watch_files else None
        self._temporary_file = Path(f"{self._watch_file}.tmp") if self._watch_file else None
        self._container_name = container_name
        self._percent = 0.0
        self._stage = "start AnnotSV"
        self._phase = "startup"
        self._annotation_started: float | None = None
        self._overlap_genes: int | None = None
        self._regulated_genes: int | None = None
        self._last_position: int | None = None
        self._last_sample_time: float | None = None
        self._speed: float | None = None

    def handle_line(self, raw_line: bytes) -> None:
        line = raw_line.decode("utf-8", errors="replace").strip()
        lowered = line.lower()
        now = time.monotonic()
        with self._lock:
            for marker, percent, stage in self._MILESTONES:
                if marker in lowered and percent >= self._percent:
                    self._percent = percent
                    self._stage = stage
                    if marker == "annotation in progress":
                        self._phase = "annotation"
                        self._annotation_started = now
                        self._last_position = None
                        self._last_sample_time = None
                        self._speed = None
                    elif marker == "writing of ":
                        self._phase = "writing"
                    break

            matched = re.search(r"\.\.\.(\d+) genes overlapped with an sv", lowered)
            if matched:
                self._overlap_genes = int(matched.group(1))
            matched = re.search(r"\.\.\.(\d+) genes regulated by", lowered)
            if matched:
                self._regulated_genes = int(matched.group(1))

    def _read_container_position(self) -> int | None:
        # AnnotSV (the container PID 1) reads outputFile.tmp sequentially during
        # OrganizeAnnotation.  /proc exposes the current byte offset without
        # modifying AnnotSV or its data files.
        script = (
            'for f in /proc/1/fd/*; do '
            'p=$(readlink "$f" 2>/dev/null) || continue; '
            'case "$p" in *.part.tsv.tmp) '
            'n=${f##*/}; '
            "awk '/^pos:/ {print $2; exit}' \"/proc/1/fdinfo/$n\"; exit 0;; "
            'esac; done'
        )
        creation_flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
        try:
            completed = subprocess.run(
                ["docker", "exec", self._container_name, "sh", "-c", script],
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                timeout=5,
                check=False,
                creationflags=creation_flags,
            )
        except (OSError, subprocess.TimeoutExpired):
            return None
        if completed.returncode != 0:
            return None
        try:
            return int(completed.stdout.decode("ascii", errors="ignore").strip())
        except ValueError:
            return None

    def _sample_annotation_stream(self, now: float) -> tuple[float, float | None] | None:
        with self._lock:
            phase = self._phase
            annotation_started = self._annotation_started
        if phase != "annotation" or annotation_started is None or now - annotation_started < 3:
            return None
        if self._temporary_file is None:
            return None
        try:
            total = self._temporary_file.stat().st_size
        except OSError:
            return None
        if total <= 0:
            return None
        position = self._read_container_position()
        if position is None:
            return None
        position = min(max(0, position), total)

        with self._lock:
            if self._phase != "annotation":
                return None
            if self._last_position is not None and position < self._last_position:
                # AnnotSV first performs a quick preparatory scan, then reopens
                # the file for the expensive annotation pass.
                self._last_position = None
                self._last_sample_time = None
                self._speed = None
            if (
                self._last_position is not None
                and self._last_sample_time is not None
                and position > self._last_position
                and now > self._last_sample_time
            ):
                current_speed = (position - self._last_position) / (now - self._last_sample_time)
                self._speed = current_speed if self._speed is None else (0.35 * current_speed + 0.65 * self._speed)
            self._last_position = position
            self._last_sample_time = now
            ratio = position / total
            self._percent = 35.0 + 55.0 * ratio
            eta = None
            if self._speed is not None and self._speed > 0 and position < total:
                eta = (total - position) / self._speed
            return ratio, eta

    def snapshot(self) -> str:
        now = time.monotonic()
        stream = self._sample_annotation_stream(now)
        with self._lock:
            percent = self._percent
            stage = self._stage
            phase = self._phase
            overlap_genes = self._overlap_genes
            regulated_genes = self._regulated_genes

        details: list[str] = [f"~{percent:.1f}%", stage]
        if stream is not None:
            ratio, eta = stream
            details.append(f"annotation stream {ratio * 100:.1f}%")
            if eta is not None:
                details.append(f"ETA {format_elapsed(eta)}")
        elif phase == "annotation" and (overlap_genes is not None or regulated_genes is not None):
            counts = []
            if overlap_genes is not None:
                counts.append(f"{overlap_genes:,} gene overlaps")
            if regulated_genes is not None:
                counts.append(f"{regulated_genes:,} regulatory overlaps")
            details.append(", ".join(counts))
        elif phase == "writing" and self._watch_file is not None:
            try:
                details.append(f"output {format_bytes(self._watch_file.stat().st_size)}")
            except OSError:
                pass
        return " | ".join(details)


def stop_child(process: subprocess.Popen[bytes]) -> None:
    if process.poll() is not None:
        return
    # Ctrl+C is delivered to every process attached to the same Windows console,
    # so normally the child is already stopping.  Give it a short grace period
    # before using terminate/kill as a final cleanup guard.
    try:
        process.wait(timeout=5)
        return
    except subprocess.TimeoutExpired:
        pass
    try:
        process.terminate()
        process.wait(timeout=5)
        return
    except (OSError, subprocess.TimeoutExpired):
        pass
    try:
        process.kill()
    except OSError:
        pass


def main() -> int:
    args = parse_args()
    started = time.monotonic()
    stopped = threading.Event()
    console_lock = threading.Lock()
    tracker = (
        AnnotSVProgress(args.watch_file, args.container_name)
        if args.progress_profile == "annotsv"
        else None
    )

    print(f"[START]  {args.label}", flush=True)
    if tracker is None:
        process = subprocess.Popen(args.command)
        output_thread = None
    else:
        process = subprocess.Popen(args.command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        def forward_output() -> None:
            assert process.stdout is not None
            while True:
                line = process.stdout.readline()
                if not line:
                    return
                tracker.handle_line(line)
                with console_lock:
                    try:
                        sys.stdout.buffer.write(line)
                        sys.stdout.buffer.flush()
                    except AttributeError:
                        sys.stdout.write(line.decode("utf-8", errors="replace"))
                        sys.stdout.flush()

        output_thread = threading.Thread(target=forward_output, name="dante-child-output", daemon=True)
        output_thread.start()

    def heartbeat() -> None:
        # Generic file/directory watches should report a state, not spam the same
        # state every N seconds.  Emit the first heartbeat and then only when the
        # watched output description changes (file appears, grows, directory changes).
        # Activities without a watch target keep the traditional periodic heartbeat.
        # Specialized progress trackers (e.g. AnnotSV) also keep periodic snapshots.
        last_watch_status: str | None = None
        while not stopped.wait(args.heartbeat_seconds):
            if process.poll() is not None:
                return
            elapsed = format_elapsed(time.monotonic() - started)
            if tracker is None:
                watched = watch_text(args.watch_file, args.watch_dir)
                if args.watch_file or args.watch_dir:
                    if watched == last_watch_status:
                        continue
                    last_watch_status = watched
                status = f"elapsed {elapsed}{watched}"
            else:
                status = f"{tracker.snapshot()} | elapsed {elapsed}"
            with console_lock:
                print(f"[ACTIVE] {args.label} | {status}", flush=True)

    thread = threading.Thread(target=heartbeat, name="dante-heartbeat", daemon=True)
    thread.start()
    try:
        return_code = process.wait()
    except KeyboardInterrupt:
        print(f"\n[PAUSED] {args.label}; derived .part files remain safe to inspect/restart.", file=sys.stderr, flush=True)
        stop_child(process)
        return 130
    finally:
        stopped.set()
        thread.join(timeout=6)
        if output_thread is not None:
            output_thread.join(timeout=2)

    elapsed = format_elapsed(time.monotonic() - started)
    if return_code == 0:
        progress = " | 100.0%" if tracker is not None else ""
        print(f"[DONE]   {args.label}{progress} | elapsed {elapsed}", flush=True)
    else:
        print(f"[FAILED] {args.label} | exit code {return_code} | elapsed {elapsed}", file=sys.stderr, flush=True)
    return return_code


def run_cli() -> None:
    try:
        code = main()
    except SystemExit as exc:
        code = exc.code if isinstance(exc.code, int) else 1
    except KeyboardInterrupt:
        code = 130
    except BaseException as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        code = 1
    raise SystemExit(code)


if __name__ == "__main__":
    run_cli()
