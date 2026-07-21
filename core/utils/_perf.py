"""Performance monitoring — wall time, CPU, memory, GPU tracking."""

import json
import os
import subprocess as _sp
import threading
import time as _time
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Optional


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
        """Convert step_num + desc into a PerformanceReport clone and append."""
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
        )
        self.steps.append(clone)

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
        """Write pretty-printed JSON, creating parent dirs if needed."""
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

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
            print(f"║  {meta_line:<{total_w - 4}}║")
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
                f"{_b}{step_label:<{c_step}}{_b}{desc:<{c_desc}}{_b}"
                f"{f'{s["wall_sec"]:.1f}s':>{c_wall}}{_b}"
                f"{f'{s["cpu_sec"]:.1f}s':>{c_cpu}}{_b}"
                f"{f'{s["peak_rss_mib"]:,.0f}':>{c_mem}}{_b}"
                f"{cells_str:>{c_cells}}{_b}"
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
        print(f"{_b}{total_str:<{left_w}}{_b}{peak_str:<{right_w}}{_b}")

        # Memory reference estimates
        n_genes_val = next((s["n_genes"] for s in steps_list if s["n_genes"]), 0)
        n_cells_total = sum(s["n_cells"] for s in steps_list) if steps_list else 0
        if n_cells_total and n_genes_val:
            mem_per_1k = summary["max_peak_rss_mib"] / (n_cells_total / 1000)
            mem_line = f" \U0001f4d0 Memory reference: ~{mem_per_1k:.2f} MiB per 1k cells at {n_genes_val:,} genes"
            print(f"{_b}{mem_line:<{total_w - 2}}{_b}")
            est = self._estimate_memory(summary["max_peak_rss_mib"], n_cells_total, n_genes_val)
            parts = []
            for k, v in est.items():
                parts.append(f"{k} × {n_genes_val:,}: ~{v:.1f} GiB")
            est_line = "    \u2192 " + "  ".join(parts)
            print(f"{_b}{est_line:<{total_w - 2}}{_b}")
            print(f"{_b}{'(linear estimate, actual varies)':<{total_w - 2}}{_b}")

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
