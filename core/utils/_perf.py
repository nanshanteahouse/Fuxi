"""Performance monitoring — wall time, CPU, memory, GPU tracking."""

import logging
import subprocess as _sp
import threading
import time as _time
from contextlib import contextmanager
from dataclasses import dataclass


@dataclass
class PerformanceReport:
    """Performance metrics for a pipeline step."""
    step: str = ""
    wall_sec: float = 0.0
    cpu_sec: float = 0.0
    peak_rss_mb: float = 0.0
    avg_cpu_pct: float = 0.0
    gpu_mem_mb: float = -1.0


@contextmanager
def monitor_performance(step_name: str = "", log=None):
    """Time CPU, memory, and GPU usage for a code block."""
    import psutil
    report = PerformanceReport(step=step_name)
    t0 = _time.time(); tcpu0 = _time.process_time()
    proc = psutil.Process()
    cpu_samples: list[float] = []
    peak_rss = 0
    stop = threading.Event()

    def _sample():
        nonlocal peak_rss
        while not stop.is_set():
            try:
                m = proc.memory_info().rss
                if m > peak_rss: peak_rss = m
                cpu_samples.append(proc.cpu_percent())
            except Exception:
                pass
            stop.wait(1.0)

    sampler = threading.Thread(target=_sample, daemon=True)
    sampler.start()
    try:
        yield report
    finally:
        stop.set(); sampler.join(timeout=5)
        dt = _time.time() - t0
        report.wall_sec = round(dt, 1)
        report.cpu_sec = round(_time.process_time() - tcpu0, 1)
        report.peak_rss_mb = round(peak_rss / 1e6, 1)
        report.avg_cpu_pct = round(sum(cpu_samples) / max(len(cpu_samples), 1), 1)
        try:
            out = _sp.check_output(["nvidia-smi","--query-gpu=memory.used",
                "--format=csv,noheader,nounits"]).decode().strip().split("\n")
            report.gpu_mem_mb = sum(float(m) for m in out)
        except Exception:
            report.gpu_mem_mb = -1.0
        if log:
            log.info("[perf] wall=%.1fs cpu=%.1fs mem=%.1fMB cpu%%=%.1f%% gpu=%.0fMB",
                     report.wall_sec, report.cpu_sec, report.peak_rss_mb,
                     report.avg_cpu_pct, report.gpu_mem_mb)
