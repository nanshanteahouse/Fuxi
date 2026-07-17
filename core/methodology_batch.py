#!/usr/bin/env python3
"""
methodology_batch.py — Batch methodology pattern extraction for all registered papers.

Scans projects/papers/*/ for papers with existing insights.yaml but missing
methodology_patterns, then runs methodology extraction in parallel using
ThreadPoolExecutor with configurable concurrency.

Usage:
    python core/methodology_batch.py                          # all papers
    python core/methodology_batch.py --dry-run                # preview only
    python core/methodology_batch.py --workers 4              # custom concurrency
    python core/methodology_batch.py --retry-failed           # retry previously failed papers

Framework code (tracked by git). Paper-specific data goes to projects/papers/ (gitignored).
"""

import os
import sys
import time
import argparse
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field

# Ensure repo root is on sys.path
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import yaml
from urllib.error import HTTPError, URLError

from core.paper_insights import PaperInsights, _build_cfg_from_env
from core.paper_converter import PmcXmlSource

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
_DEFAULT_WORKERS = 8
_PAPERS_ROOT = Path(_REPO_ROOT) / "projects" / "papers"


@dataclass
class BatchStats:
    """Track batch processing statistics."""
    total: int = 0
    success: int = 0
    skipped: int = 0
    failed: int = 0
    failed_papers: list[str] = field(default_factory=list)
    start_time: float = 0.0


@dataclass
class PaperJob:
    """A single paper to process."""
    dir_name: str
    xml_path: Path
    insights_path: Path


def _build_cfg():
    """Build AI config from environment."""
    return _build_cfg_from_env()

def discover_jobs() -> list[PaperJob]:
    """Find papers that have XML source and insights.yaml but lack methodology_patterns.

    Returns list of PaperJob ready for processing.
    """
    jobs: list[PaperJob] = []

    if not _PAPERS_ROOT.is_dir():
        logger.warning("Papers root not found: %s", _PAPERS_ROOT)
        return jobs

    for paper_dir in sorted(_PAPERS_ROOT.iterdir()):
        if not paper_dir.is_dir():
            continue

        insights_path = paper_dir / "insights.yaml"
        if not insights_path.is_file():
            continue

        # Check if methodology already extracted
        with open(insights_path) as f:
            data = yaml.safe_load(f)
        if data and "methodology_patterns" in data:
            logger.debug("Skip (already has methodology): %s", paper_dir.name)
            continue

        # Find XML source
        xml_files = list(paper_dir.glob("*.xml"))
        if not xml_files:
            logger.debug("Skip (no XML): %s", paper_dir.name)
            continue

        jobs.append(PaperJob(
            dir_name=paper_dir.name,
            xml_path=xml_files[0],
            insights_path=insights_path,
        ))

    return jobs


def process_paper(job: PaperJob, cfg) -> str:
    """Process a single paper: extract methodology patterns and write insights.yaml.

    Returns status string: "OK", "SKIP", or "FAIL: <reason>".
    """
    try:
        source = PmcXmlSource(xml_path=str(job.xml_path))
        result = PaperInsights().run(
            source=source,
            cfg=cfg,
            output_path=str(job.insights_path),
            force=True,
            extract_methodology=True,
        )
        if result == "SKIPPED":
            return "SKIP"
        return "OK"
    except (HTTPError, URLError) as e:
        return f"FAIL: HTTP/URL error: {e}"
    except RuntimeError as e:
        return f"FAIL: Runtime error: {e}"
    except Exception as e:
        return f"FAIL: {type(e).__name__}: {e}"


def run_batch(workers: int = _DEFAULT_WORKERS, dry_run: bool = False) -> BatchStats:
    """Run methodology extraction for all discovered papers in parallel.

    Args:
        workers: Number of parallel threads.
        dry_run: If True, only list jobs without processing.

    Returns:
        BatchStats with processing results.
    """
    jobs = discover_jobs()
    stats = BatchStats(total=len(jobs))
    stats.start_time = time.time()

    print(f"\n=== Methodology Batch Extraction ===")
    print(f"Papers to process: {len(jobs)}")
    print(f"Workers: {workers}")
    print(f"Dry run: {dry_run}")
    print()

    if dry_run:
        for i, job in enumerate(jobs, 1):
            print(f"  [{i:3d}/{len(jobs)}] {job.dir_name[:70]}")
        print(f"\n  Dry run complete. {len(jobs)} papers would be processed.")
        return stats

    if not jobs:
        print("  No papers to process.")
        return stats

    cfg = _build_cfg()

    with ThreadPoolExecutor(max_workers=workers) as executor:
        # Submit all jobs
        future_map = {
            executor.submit(process_paper, job, cfg): job
            for job in jobs
        }

        # Collect results as they complete
        for i, future in enumerate(as_completed(future_map), 1):
            job = future_map[future]
            elapsed = time.time() - stats.start_time
            rate = i / elapsed if elapsed > 0 else 0

            try:
                status = future.result()
            except Exception as e:
                status = f"FAIL: {type(e).__name__}: {e}"

            if status == "OK":
                stats.success += 1
            elif status.startswith("SKIP"):
                stats.skipped += 1
            else:
                stats.failed += 1
                stats.failed_papers.append(f"{job.dir_name}: {status}")

            # Progress with ETA
            eta = (len(jobs) - i) / rate if rate > 0 else 0
            print(
                f"\r  [{i:3d}/{len(jobs)}] "
                f"OK:{stats.success} FAIL:{stats.failed} SKIP:{stats.skipped} "
                f"| {rate:.1f}/s | ETA: {eta:.0f}s | "
                f"{job.dir_name[:50]}...",
                end="",
                flush=True,
            )

    print()  # newline after progress

    return stats


def print_summary(stats: BatchStats) -> None:
    """Print final summary report."""
    elapsed = time.time() - stats.start_time
    print(f"\n{'='*60}")
    print(f"Batch Complete")
    print(f"{'='*60}")
    print(f"  Total:    {stats.total}")
    print(f"  Success:  {stats.success}")
    print(f"  Failed:   {stats.failed}")
    print(f"  Skipped:  {stats.skipped}")
    print(f"  Elapsed:  {elapsed:.1f}s ({elapsed/stats.total:.1f}s/paper avg)")
    if stats.failed_papers:
        print(f"\n  Failed papers:")
        for fp in stats.failed_papers[:10]:
            print(f"    - {fp}")
        if len(stats.failed_papers) > 10:
            print(f"    ... and {len(stats.failed_papers) - 10} more")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch methodology pattern extraction for all registered papers."
    )
    parser.add_argument(
        "--dry-run", "-n", action="store_true", default=False,
        help="Preview papers to process without running extraction."
    )
    parser.add_argument(
        "--workers", "-w", type=int, default=_DEFAULT_WORKERS,
        help=f"Number of parallel workers (default: {_DEFAULT_WORKERS})."
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", default=False,
        help="Enable debug logging."
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s [%(name)s] %(message)s",
        stream=sys.stderr,
    )
    logging.getLogger("httpx").setLevel(logging.WARNING)
    logging.getLogger("httpcore").setLevel(logging.WARNING)
    logging.getLogger("openai").setLevel(logging.WARNING)

    stats = run_batch(workers=args.workers, dry_run=args.dry_run)
    if not args.dry_run:
        print_summary(stats)


if __name__ == "__main__":
    main()
