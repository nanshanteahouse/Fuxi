"""End-to-end integration tests for 05_annotate_major quality report paths.

Verifies that each annotation path (unified KB, AI, score_genes) produces
the expected 05_annotation_quality.json output.

These tests are marked ``integration`` (slower, opt-in via ``-m integration``).
"""

import importlib.util
import json
import os
import sys
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from anndata import AnnData

# ── Ensure repo root is on sys.path (conftest.py also does this) ──────
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

# ── Load the 05_annotate_major module via file path ───────────────────
_STEP_PATH = os.path.join(_REPO_ROOT, "rna", "steps", "05_annotate_major.py")
_spec = importlib.util.spec_from_file_location("rna.steps._05_quality_report_paths", _STEP_PATH)
assert _spec is not None and _spec.loader is not None, f"Could not load {_STEP_PATH}"
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

_update_quality_report_pass_rate = _mod._update_quality_report_pass_rate
ai_annotate = _mod.ai_annotate
score_genes_mode = _mod.score_genes_mode


# ── Shared helpers ────────────────────────────────────────────────────


def _make_synthetic_adata(with_raw: bool = False) -> AnnData:
    """Create a synthetic AnnData with leiden clusters and optional .raw.

    Fixture duplicated from test_05_quality_report.py to keep test modules
    independent.
    """
    rng = np.random.RandomState(42)
    n_cells = 50
    n_genes = 20

    x = rng.poisson(lam=1.0, size=(n_cells, n_genes)).astype(np.float32)
    adata = AnnData(x)
    adata.var_names = [f"GENE_{i}" for i in range(n_genes)]
    adata.obs["leiden"] = rng.choice(["0", "1", "2", "3"], n_cells)
    adata.obsm["X_umap"] = rng.randn(n_cells, 2)

    if with_raw:
        n_raw = n_genes + 10
        raw_x = rng.poisson(lam=1.0, size=(n_cells, n_raw)).astype(np.float32)
        raw = AnnData(raw_x)
        raw.var_names = [f"GENE_{i}" for i in range(n_raw)]
        adata.raw = raw

    return adata


def _make_mock_std() -> MagicMock:
    """Create a minimal mock StandardOntology for marker validation.

    Minimal mock KB hierarchy duplicated from test_05_quality_report.py to
    keep test modules independent.
    """
    std = MagicMock()

    def _validate(adata, **kwargs):
        clusters = sorted(adata.obs["leiden"].unique(), key=lambda x: int(x))
        results = []
        for c in clusters:
            if c in ("0", "1"):
                status = "PASS"
            elif c == "2":
                status = "MARGINAL"
            else:
                status = "FAIL"
            results.append({"cluster": c, "status": status})
        return results

    std.validate.side_effect = _validate
    return std


def _apply_validation(adata: AnnData, std: MagicMock) -> None:
    """Apply mock standardizer validation to set marker_validation column.

    Mirrors the pattern used in main() for all three annotation paths.
    """
    validation_results = std.validate(adata)
    validation_map = {r["cluster"]: r["status"] for r in validation_results}
    adata.obs["marker_validation"] = (
        adata.obs["leiden"].astype(str).map(lambda c: validation_map.get(c, "NO_ONTOLOGY"))
    )


def _setup_rgg_result(adata: AnnData) -> None:
    """Populate adata.uns['rank_genes_groups'] in scanpy's recarray format.

    Required by ai_annotate to read per-cluster marker genes via
    ``sc.get.rank_genes_groups_df``.
    """
    groups = sorted(adata.obs["leiden"].unique(), key=lambda x: int(x))
    n_top = 5
    dtype = [(str(g), "O") for g in groups]

    names_rows = []
    scores_rows = []
    pvals_rows = []
    pvals_adj_rows = []
    lfc_rows = []

    rng = np.random.RandomState(99)
    for i in range(n_top):
        names_rows.append(tuple(f"MARKER_{g}_{i}" for g in groups))
        scores_rows.append(tuple(float(rng.uniform(0.5, 5.0)) for _ in groups))
        pvals_rows.append(tuple(float(rng.uniform(1e-10, 0.05)) for _ in groups))
        pvals_adj_rows.append(tuple(float(rng.uniform(1e-8, 0.1)) for _ in groups))
        lfc_rows.append(tuple(float(rng.randn()) for _ in groups))

    adata.uns["rank_genes_groups"] = {
        "names": np.rec.array(names_rows, dtype=dtype),
        "scores": np.rec.array(scores_rows, dtype=dtype),
        "pvals": np.rec.array(pvals_rows, dtype=dtype),
        "pvals_adj": np.rec.array(pvals_adj_rows, dtype=dtype),
        "logfoldchanges": np.rec.array(lfc_rows, dtype=dtype),
        "params": {
            "groupby": "leiden",
            "method": "wilcoxon",
            "use_raw": True,
        },
    }


# ── Tests ─────────────────────────────────────────────────────────────


@pytest.mark.integration
def test_e2e_unified_path(tmp_path: os.PathLike) -> None:
    """Unified KB annotation path produces 05_annotation_quality.json.

    Given: synthetic adata with cell_type and marker_validation columns
           (simulating unified_annotate + std.validate output).
    When:  _update_quality_report_pass_rate is called.
    Then:  05_annotation_quality.json is created in cfg.table_dir.
    """
    adata = _make_synthetic_adata()

    # Simulate cell_type assignment from unified KB annotation
    leiden_map = {"0": "T cell", "1": "B cell", "2": "NK cell", "3": "Monocyte"}
    adata.obs["cell_type"] = adata.obs["leiden"].map(leiden_map).astype("category")

    # Apply marker validation (as in main() unified path)
    std = _make_mock_std()
    _apply_validation(adata, std)

    cfg = MagicMock()
    cfg.table_dir = str(tmp_path)

    _update_quality_report_pass_rate(adata, cfg)

    quality_path = os.path.join(cfg.table_dir, "05_annotation_quality.json")
    assert os.path.exists(quality_path), f"Expected quality report at {quality_path}"


@pytest.mark.integration
@patch("core.ai.caller.ai_query")
def test_e2e_ai_path(mock_ai_query: MagicMock, tmp_path: os.PathLike) -> None:
    """AI annotation path produces JSON with strict_pass_rate field.

    Given: synthetic adata + mocked ai_query returning valid annotation.
    When:  ai_annotate runs successfully followed by quality report update.
    Then:  05_annotation_quality.json contains the strict_pass_rate field.
    """
    adata = _make_synthetic_adata(with_raw=False)
    cfg = MagicMock()
    cfg.ai = MagicMock()
    cfg.ai.max_tokens = 4096
    cfg.ai.model = "test-model"
    cfg.tissue = "retina"
    cfg.species = "human"
    cfg.table_dir = str(tmp_path)
    cfg.figure_dir = str(tmp_path)
    logger = MagicMock()

    # Pre-populate rank_genes_groups so that ai_annotate can read it
    # (rank_genes_groups will be mocked to no-op, preserving pre-set data)
    _setup_rgg_result(adata)

    groups = sorted(adata.obs["leiden"].unique(), key=lambda x: int(x))
    ai_return = json.dumps(
        {
            str(g): {
                "cell_type": {
                    "0": "T cell",
                    "1": "B cell",
                    "2": "NK cell",
                    "3": "Monocyte",
                }[str(g)],
                "state": "active",
                "subtype": "CD8+",
                "confidence": "high",
                "reasoning": "marker expression",
            }
            for g in groups
        }
    )
    mock_ai_query.return_value = ai_return

    with (
        patch.object(_mod.sc.tl, "rank_genes_groups"),
        patch.object(_mod, "safe_plot"),
    ):
        ann_result = ai_annotate(adata, cfg, logger)

    assert ann_result is not None, "ai_annotate should return annotations dict"

    # Apply marker validation (as in main() AI path)
    std = _make_mock_std()
    _apply_validation(adata, std)

    _update_quality_report_pass_rate(adata, cfg)

    quality_path = os.path.join(cfg.table_dir, "05_annotation_quality.json")
    assert os.path.exists(quality_path), f"Expected quality report at {quality_path}"
    with open(quality_path) as f:
        quality = json.load(f)
    assert "strict_pass_rate" in quality, "Quality report should include strict_pass_rate"


@pytest.mark.integration
def test_e2e_score_genes_path(tmp_path: os.PathLike) -> None:
    """Score_genes annotation path produces 05_annotation_quality.json.

    Given: synthetic adata with .raw + marker_dict configured.
    When:  score_genes_mode runs followed by quality report update.
    Then:  05_annotation_quality.json is created in cfg.table_dir.
    """
    adata = _make_synthetic_adata(with_raw=True)
    cfg = MagicMock()
    cfg.table_dir = str(tmp_path)
    cfg.figure_dir = str(tmp_path)
    cfg.marker = MagicMock()
    cfg.marker.marker_dict = {
        "T cell": ["GENE_0", "GENE_1"],
        "B cell": ["GENE_2", "GENE_3"],
        "NK cell": ["GENE_4", "GENE_5"],
    }
    cfg.marker.subcluster_types = []
    logger = MagicMock()

    with patch.object(_mod, "safe_plot"):
        score_genes_mode(adata, cfg, logger)

    # Apply marker validation (as in main() score_genes path)
    std = _make_mock_std()
    _apply_validation(adata, std)

    _update_quality_report_pass_rate(adata, cfg)

    quality_path = os.path.join(cfg.table_dir, "05_annotation_quality.json")
    assert os.path.exists(quality_path), f"Expected quality report at {quality_path}"
