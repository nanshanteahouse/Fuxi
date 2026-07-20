"""Pytest configuration — adds repo root to sys.path for clean imports.

Any test file under ``tests/`` can then do::

    from core.config.schema import CFG
    from core.utils import resolve_config
    from core.annotation.scoring import score_cluster_against_kb
    from rna.utils.evidence_fusion import fuse_all_clusters
    from core.cluster.evaluation import find_pareto_frontier

without manual sys.path manipulation.
"""

import os
import sys

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
