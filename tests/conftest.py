"""Pytest configuration — adds repo root to sys.path for clean imports.

Any test file under ``tests/`` can then do::

    from core.config.schema import CFG
    from core.utils import resolve_config
    from core.annotation.scoring import score_cluster_against_kb
    from rna.utils.evidence_fusion import fuse_all_clusters
    from core.cluster.evaluation import find_pareto_frontier
    from core.paper.registry import load_master_registry

without manual sys.path manipulation.
"""

import os
import sys

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)


import pytest  # noqa: E402  # after sys.path setup
import scanpy as sc  # noqa: E402

_FIXTURE_DIR = os.path.join(_REPO_ROOT, "tests", "fixtures")


@pytest.fixture(scope="session")
def synthetic_adata():
    """5k well-separated synthetic fixture, cached per session."""
    path = os.path.join(_FIXTURE_DIR, "synthetic_rna.h5ad")
    if not os.path.exists(path):
        pytest.skip(f"Fixture not found: {path}. Run tests/fixtures/generate_synthetic.py first.")
    return sc.read(path)


@pytest.fixture(scope="session")
def synthetic_overlapping_adata():
    path = os.path.join(_FIXTURE_DIR, "synthetic_overlapping.h5ad")
    if not os.path.exists(path):
        pytest.skip(f"Fixture not found: {path}")
    return sc.read(path)


@pytest.fixture(scope="session")
def synthetic_severely_overlapping_adata():
    path = os.path.join(_FIXTURE_DIR, "synthetic_severely_overlapping.h5ad")
    if not os.path.exists(path):
        pytest.skip(f"Fixture not found: {path}")
    return sc.read(path)


@pytest.fixture(scope="session")
def synthetic_100k_adata():
    path = os.path.join(_FIXTURE_DIR, "synthetic_100k.h5ad")
    if not os.path.exists(path):
        pytest.skip(f"Fixture not found: {path}")
    return sc.read(path)


@pytest.fixture(scope="session")
def shekhar_5k_adata():
    path = os.path.join(_FIXTURE_DIR, "shekhar_5k_subsample.h5ad")
    if not os.path.exists(path):
        pytest.skip(f"Fixture not found: {path}")
    return sc.read(path)


@pytest.fixture(scope="session")
def shekhar_20k_adata():
    path = os.path.join(_FIXTURE_DIR, "shekhar_20k_subsample.h5ad")
    if not os.path.exists(path):
        pytest.skip(f"Fixture not found: {path}")
    return sc.read(path)
