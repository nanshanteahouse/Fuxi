"""Performance monitoring — wall time, CPU, memory, GPU tracking."""

import json
import os
import subprocess as _sp
import threading
import time as _time
import unicodedata
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Optional

from core.utils._gpu import is_gpu_active


def _display_width(s: str) -> int:
    """Terminal display width accounting for CJK fullwidth and emoji characters."""
    w = 0
    for ch in s:
        ea = unicodedata.east_asian_width(ch)
        w += 2 if ea in ("W", "F") else 1
    return w


def _pad(s: str, width: int, align: str = "<") -> str:
    """Pad string to `width` display columns, respecting terminal char width."""
    dw = _display_width(s)
    if dw >= width:
        return s
    padding = width - dw
    if align == "<":
        return s + " " * padding
    elif align == ">":
        return " " * padding + s
    else:  # "^"
        left = padding // 2
        right = padding - left
        return " " * left + s + " " * right


@dataclass
class PerformanceReport:
    """Performance metrics for a pipeline step."""

    step: str = ""
    wall_sec: float = 0.0
    cpu_sec: float = 0.0
    peak_rss_mib: float = 0.0  # MiB = 1024² bytes
    avg_cpu_pct: float = 0.0
    gpu_mem_mb: float = -1.0
    n_cells: int = 0
    n_genes: int = 0
    checkpoint_mib: float = 0.0  # MiB = 1024² bytes
    exit_status: str = "completed"  # "completed" | "killed" | "failed" | "skipped"


@contextmanager
def monitor_performance(step_name: str = "", log=None, child_pid: Optional[int] = None):
    """Time CPU, memory, and GPU usage for a code block.

    Parameters
    ----------
    step_name : str
        Label for this performance measurement.
    log : logger, optional
        Logger to write results to.
    child_pid : int, optional
        If provided, the child process (pid) is tracked instead of the parent.
    """
    import psutil

    report = PerformanceReport(step=step_name)
    t0 = _time.time()
    if child_pid is not None:
        proc = psutil.Process(child_pid)
        tcpu0 = proc.cpu_times()
    else:
        proc = psutil.Process()
        tcpu0_parent = _time.process_time()
    cpu_samples: list[float] = []
    peak_rss = 0
    stop = threading.Event()

    def _sample():
        nonlocal peak_rss
        while not stop.is_set():
            try:
                m = proc.memory_info().rss
                if m > peak_rss:
                    peak_rss = m  # track new peak
                cpu_samples.append(proc.cpu_percent())
            except Exception:
                pass
            stop.wait(1.0)

    sampler = threading.Thread(target=_sample, daemon=True)
    sampler.start()
    try:
        yield report
    finally:
        stop.set()
        sampler.join(timeout=5)
        dt = _time.time() - t0
        report.wall_sec = round(dt, 1)
        if child_pid is not None:
            try:
                tcpu1 = proc.cpu_times()
                report.cpu_sec = round(
                    (tcpu1.user + tcpu1.system) - (tcpu0.user + tcpu0.system), 1
                )
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                report.cpu_sec = 0.0
        else:
            report.cpu_sec = round(_time.process_time() - tcpu0_parent, 1)
        # psutil memory_info().rss returns bytes → convert to MiB (1024²)
        report.peak_rss_mib = round(peak_rss / (1024 * 1024), 1)
        report.avg_cpu_pct = round(sum(cpu_samples) / max(len(cpu_samples), 1), 1)
        if is_gpu_active():
            try:
                out = (
                    _sp.check_output(
                        ["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"]
                    )
                    .decode()
                    .strip()
                    .split("\n")
                )
                report.gpu_mem_mb = sum(float(m) for m in out)
            except Exception:
                report.gpu_mem_mb = -1.0
        else:
            report.gpu_mem_mb = 0.0
        if log:
            log.info(
                "[perf] wall=%.1fs cpu=%.1fs mem=%.1fMiB cpu%%=%.1f%% gpu=%.0fMB",
                report.wall_sec,
                report.cpu_sec,
                report.peak_rss_mib,
                report.avg_cpu_pct,
                report.gpu_mem_mb,
            )


class PerformanceSummary:
    """Aggregated performance summary across pipeline steps."""

    def __init__(self) -> None:
        self.pipeline_info: dict = {}
        self.steps: list[PerformanceReport] = []

    def add_step(self, step_num: str, desc: str, perf: PerformanceReport) -> None:
        """Upsert a step entry by step_num — replace if rerun, append if new.

        Ensures ``perf_report.json`` survives ``--step N`` reruns and
        ``--resume`` without losing history: when the same step is rerun,
        the new measurements replace the previous entry instead of being
        appended as duplicates.
        """
        clone = PerformanceReport(
            step=f"{step_num} {desc}",
            wall_sec=perf.wall_sec,
            cpu_sec=perf.cpu_sec,
            peak_rss_mib=perf.peak_rss_mib,
            avg_cpu_pct=perf.avg_cpu_pct,
            gpu_mem_mb=perf.gpu_mem_mb,
            n_cells=perf.n_cells,
            n_genes=perf.n_genes,
            checkpoint_mib=perf.checkpoint_mib,
            exit_status=perf.exit_status,
        )
        key = str(step_num)
        for i, s in enumerate(self.steps):
            if s.step.split(" ", 1)[0] == key:
                self.steps[i] = clone
                return
        self.steps.append(clone)

    @classmethod
    def load_existing(cls, path: str) -> Optional["PerformanceSummary"]:
        """Load a previously saved ``perf_report.json``.

        Returns ``None`` if the file is missing or cannot be parsed, so the
        caller can fall back to a fresh summary. Used by the runner when
        resuming or running a single step (``--step N``) to preserve
        historical step metrics instead of overwriting them.
        """
        if not os.path.exists(path):
            return None
        try:
            with open(path) as f:
                data = json.load(f)
        except (OSError, json.JSONDecodeError):
            return None
        inst = cls()
        inst.pipeline_info = data.get("pipeline", {}) or {}
        for s in data.get("steps", []) or []:
            inst.steps.append(
                PerformanceReport(
                    step=s.get("step", ""),
                    wall_sec=s.get("wall_sec", 0.0),
                    cpu_sec=s.get("cpu_sec", 0.0),
                    peak_rss_mib=s.get("peak_rss_mib", 0.0),
                    avg_cpu_pct=s.get("avg_cpu_pct", 0.0),
                    gpu_mem_mb=s.get("gpu_mem_mb", -1.0),
                    n_cells=s.get("n_cells", 0),
                    n_genes=s.get("n_genes", 0),
                    checkpoint_mib=s.get("checkpoint_mib", 0.0),
                    exit_status=s.get("exit_status", "completed"),
                )
            )
        return inst

    def to_dict(self) -> dict:
        """Full JSON-ready dict with pipeline, steps, and summary keys."""
        steps_list: list[dict] = []
        for s in self.steps:
            steps_list.append(
                {
                    "step": s.step,
                    "wall_sec": s.wall_sec,
                    "cpu_sec": s.cpu_sec,
                    "peak_rss_mib": s.peak_rss_mib,
                    "avg_cpu_pct": s.avg_cpu_pct,
                    "gpu_mem_mb": s.gpu_mem_mb,
                    "n_cells": s.n_cells,
                    "n_genes": s.n_genes,
                    "checkpoint_mib": s.checkpoint_mib,
                    "exit_status": s.exit_status,
                }
            )

        if steps_list:
            total_wall = sum(s["wall_sec"] for s in steps_list)
            max_rss = max(s["peak_rss_mib"] for s in steps_list)
            max_rss_step = max(steps_list, key=lambda x: x["peak_rss_mib"])["step"]
        else:
            total_wall = 0.0
            max_rss = 0.0
            max_rss_step = ""

        return {
            "pipeline": self.pipeline_info,
            "steps": steps_list,
            "summary": {
                "n_steps": len(steps_list),
                "total_wall_sec": total_wall,
                "max_peak_rss_mib": max_rss,
                "max_peak_rss_step": max_rss_step,
            },
        }

    def save_json(self, path: str) -> None:
        """Write pretty-printed JSON atomically.

        Writes to ``<path>.tmp.<pid>`` then ``os.replace`` for an atomic
        same-fs rename, so an interrupted save never leaves a partially
        written ``perf_report.json`` (which would also break the next
        ``PerformanceSummary.load_existing`` call).
        """
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        tmp_path = f"{path}.tmp.{os.getpid()}"
        try:
            with open(tmp_path, "w") as f:
                json.dump(self.to_dict(), f, indent=2)
            os.replace(tmp_path, path)  # atomic same-fs rename
        except Exception:
            # Clean up tmp file on failure so it doesn't leak
            try:
                os.unlink(tmp_path)
            except OSError:
                pass
            raise

    def print_terminal_summary(
        self, n_jobs: int = 0, modality: str = "", config_path: str = ""
    ) -> None:
        """Print a bordered table with step details and memory reference estimate."""
        data = self.to_dict()
        steps_list = data["steps"]
        summary = data["summary"]

        if not steps_list:
            print("No steps recorded.")
            return

        # Column widths
        c_step = 6
        c_desc = 26
        c_wall = 8
        c_cpu = 8
        c_mem = 8
        c_cells = 7

        total_w = c_step + c_desc + c_wall + c_cpu + c_mem + c_cells + 7  # borders

        # Pipeline header
        if modality or config_path:
            meta_parts = []
            if modality:
                meta_parts.append(modality)
            if config_path:
                meta_parts.append(config_path)
            meta_line = "  |  ".join(meta_parts)
            print(f"╔{'═' * (total_w - 2)}╗")
            print(f"║  {_pad(meta_line, total_w - 4)}║")
            print(f"╚{'═' * (total_w - 2)}╝")

        def _h(w):
            return "─" * w

        _b = "│"

        top = f"┌{_h(c_step)}┬{_h(c_desc)}┬{_h(c_wall)}┬{_h(c_cpu)}┬{_h(c_mem)}┬{_h(c_cells)}┐"
        sep = f"├{_h(c_step)}┼{_h(c_desc)}┼{_h(c_wall)}┼{_h(c_cpu)}┼{_h(c_mem)}┼{_h(c_cells)}┤"
        header = (
            f"{_b}{'Step':<{c_step}}{_b}{'Description':<{c_desc}}{_b}"
            f"{'Wall':<{c_wall}}{_b}{'CPU':<{c_cpu}}{_b}{'Mem(MiB)':<{c_mem}}{_b}{'Cells':<{c_cells}}{_b}"
        )
        bot_sep = f"├{_h(c_step)}┴{_h(c_desc)}┼{_h(c_wall)}┼{_h(c_cpu)}┼{_h(c_mem)}┼{_h(c_cells)}┤"
        bot = f"└{_h(total_w - 2)}┘"

        print(top)
        print(header)
        print(sep)

        for s in steps_list:
            step_label, desc = (s["step"].split(" ", 1) + [""])[:2]
            cells_str = f"{s['n_cells'] / 1000:.1f}k" if s["n_cells"] else ""
            print(
                f"{_b}{_pad(step_label, c_step)}{_b}{_pad(desc, c_desc)}{_b}"
                f"{_pad(f'{s["wall_sec"]:.1f}s', c_wall, '>')}{_b}"
                f"{_pad(f'{s["cpu_sec"]:.1f}s', c_cpu, '>')}{_b}"
                f"{_pad(f'{s["peak_rss_mib"]:,.0f}', c_mem, '>')}{_b}"
                f"{_pad(cells_str, c_cells, '>')}{_b}"
            )

        print(bot_sep)

        # Total / Peak row (two merged cells)
        total_str = f"Total: {summary['total_wall_sec']:.1f}s wall"
        if n_jobs:
            total_str += f" ({n_jobs} jobs)"
        peak_str = f"Peak: {summary['max_peak_rss_mib']:,.0f} MiB"
        if summary.get("max_peak_rss_step"):
            peak_str += f" ({summary['max_peak_rss_step']})"
        left_w = c_step + 1 + c_desc  # merged left cell width
        right_w = c_wall + 1 + c_cpu + 1 + c_mem + 1 + c_cells  # merged right cell width
        print(f"{_b}{_pad(total_str, left_w)}{_b}{_pad(peak_str, right_w)}{_b}")

        # Memory reference estimates
        n_genes_val = next((s["n_genes"] for s in steps_list if s["n_genes"]), 0)
        n_cells_total = sum(s["n_cells"] for s in steps_list) if steps_list else 0
        if n_cells_total and n_genes_val:
            mem_per_1k = summary["max_peak_rss_mib"] / (n_cells_total / 1000)
            mem_line = f"[mem] Memory reference: ~{mem_per_1k:.2f} MiB per 1k cells at {n_genes_val:,} genes"
            print(f"{_b}{_pad(mem_line, total_w - 2)}{_b}")
            est = self._estimate_memory(summary["max_peak_rss_mib"], n_cells_total, n_genes_val)
            parts = []
            for k, v in est.items():
                parts.append(f"{k} x {n_genes_val:,}: ~{v:.1f} GiB")
            est_line = "    -> " + "  ".join(parts)
            print(f"{_b}{_pad(est_line, total_w - 2)}{_b}")
            print(f"{_b}{_pad('(linear estimate, actual varies)', total_w - 2)}{_b}")

        print(bot)

    @staticmethod
    def _estimate_memory(peak_rss_mib: float, n_cells: int, n_genes: int) -> dict[str, float]:
        """Estimate memory requirements at scale based on current measurement.

        Returns estimated GiB for 50k/100k/200k/500k cells at same gene density.
        """
        if n_cells == 0 or n_genes == 0:
            return {}
        per_cell_gene_mib = peak_rss_mib / (n_cells * n_genes)
        return {
            "50k": round(50_000 * n_genes * per_cell_gene_mib / 1024, 1),
            "100k": round(100_000 * n_genes * per_cell_gene_mib / 1024, 1),
            "200k": round(200_000 * n_genes * per_cell_gene_mib / 1024, 1),
            "500k": round(500_000 * n_genes * per_cell_gene_mib / 1024, 1),
        }


@contextmanager
def timed_substep(name: str, log=None, *, log_level=None):
    """Time a sub-step within a pipeline step.

    Surfaces hidden bottlenecks inside a single pipeline step (PCA, Harmony,
    UMAP sweep, batch_diag, ...) by emitting clean log lines of the form
    ``[substep] <name> took <N.N>s``. Replaces manual timestamp archaeology
    when triaging step wall time.

    Usage::

        with timed_substep("Harmony", log=log):
            sc.external.pp.harmony_integrate(adata, batch_key)

    Optionally bind the yielded dict to read the wall time programmatically::

        with timed_substep("PCA", log=log) as t:
            sc.pp.pca(adata)
        log.info("PCA used %s sec", t["wall_sec"])

    Args:
        name: sub-step label (logged verbatim).
        log: optional logger with a ``.log(level, fmt, *args)`` method
            (standard ``logging.Logger``). If None, no log line is emitted
            but the timing is still captured on the yielded dict.
        log_level: optional explicit level (default ``logging.INFO``).

    Yields:
        dict with ``name`` and ``wall_sec`` (the latter filled on exit).
    """
    import logging as _logging

    t0 = _time.time()
    info = {"name": name, "wall_sec": 0.0}
    try:
        yield info
    finally:
        info["wall_sec"] = round(_time.time() - t0, 1)
        if log is not None:
            level = log_level if log_level is not None else _logging.INFO
            log.log(level, "[substep] %s took %.1fs", name, info["wall_sec"])


def record_memory_skip(step: str, operation: str, reason: str, cfg=None, log=None) -> None:
    """Record that a pipeline step skipped an operation due to memory_policy.

    Appends one JSON line to ``<results_dir>/memory_skips.jsonl`` so users can
    audit silent skips across a pipeline run. The file is append-only across
    ``--step N`` / ``--resume`` invocations — a single grep gives the full
    picture without re-parsing step logs.

    Each line is a JSON object::

        {"timestamp": "2026-07-25T14:35:01", "step": "03_integrate",
         "operation": "regress_out", "reason": "memory_policy=balanced ..."}

    Usage from a step::

        if _skip_regress:
            log.info("memory_policy=%s — skipping regress_out", mem_policy)
            record_memory_skip(
                step="03_integrate",
                operation="regress_out",
                reason=f"memory_policy={mem_policy} would dense-allocate",
                cfg=cfg, log=log,
            )

    Args:
        step: short step identifier (e.g. ``"03_integrate"``).
        operation: what was skipped (e.g. ``"regress_out"``, ``"PCA dense"``).
        reason: human-readable reason, ideally including the policy value.
        cfg: resolved config — needs a ``results_dir`` attribute. If missing,
            the skip is logged but not recorded to disk.
        log: optional logger for diagnostics about the recording itself.
    """
    import datetime as _dt

    record = {
        "timestamp": _dt.datetime.now().isoformat(timespec="seconds"),
        "step": step,
        "operation": operation,
        "reason": reason,
    }

    results_dir = getattr(cfg, "results_dir", None) if cfg is not None else None
    if not results_dir:
        if log is not None:
            log.warning("[memory-skip] could not record %s/%s: no results_dir", step, operation)
        return

    path = os.path.join(results_dir, "memory_skips.jsonl")
    try:
        os.makedirs(results_dir, exist_ok=True)
        with open(path, "a") as f:
            f.write(json.dumps(record) + "\n")
    except OSError as e:
        if log is not None:
            log.warning("[memory-skip] failed to write %s: %s", path, e)
