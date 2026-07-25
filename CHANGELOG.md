# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
(via commit convention — versions not yet tagged).

## [Unreleased]

### Cluster-overhaul plan (plan #47)

A comprehensive rewrite of the clustering pipeline: 8 bug fixes, 6 new
features, significant GPU-accelerated performance improvements, and a
breaking config template cleanup.  The grid search engine now supports
three dispatch modes, progressive funnel search for 1M-scale datasets,
target-K binary search, and per-entry multi-metric enrichment that
powers a 3-tier resolution recommendation.

#### Bug Fixes

- **B1/B2 — stability NaN.** `_compute_stability` returned NaN when
  `n_seeds <= 1` or when all seeds failed.  Downstream multi-metric
  scoring now detects NaN entries and degrades gracefully: a low-variance
  guard (range < 0.01) silently removes the metric from the composite.
  *Before:* composite score was NaN, grid search aborted.  *After:*
  metric is dropped, remaining weights are renormalised, search proceeds.

- **C1/C6 — coherence NaN.** `_compute_cluster_coherence` returned NaN
  for empty `per_cell_scores` or when no cell type had valid (non-None)
  entries.  The multi-metric selector now degrades from 5-metric to
  3-metric (silhouette + stability + coherence) when `splitting_gain`
  and `kb_annotatable_rate` are absent, and down to silhouette-only
  when all auxiliary metrics are NaN.  *Before:* a single NaN in any
  entry poisoned the entire composite.  *After:* metric-level degradation
  with per-metric low-variance guards.

- **D1/D2/D4 — DE-gated selection for subtype data.** Added
  `_select_de_gated` for tissue subtypes where silhouette scores are
  flat across resolutions.  Finds the highest resolution where every
  cluster still expresses >= `multi_metric_de_gate_threshold` (default
  25) one-vs-rest DE genes.  *Before:* multi-metric on subtype data
  chose arbitrary mid-resolution parameters with no biological sanity
  check.  *After:* DE-gated selection ensures each cluster is
  transcriptionally distinct, following the Shekhar 2016 inverted
  pattern.

- **Dead code removed.** Removed the `if not results_summary:
  log.critical(...); sys.exit(1)` block from both `rna/steps/04_cluster_umap.py`
  and `spatial/steps/04_cluster.py`.  *Before:* step 04 hard-exited on
  empty grid results without cleanup.  *After:* empty results propagate
  normally; the pipeline runner handles the error state.

- **getattr default fix.** `on_non_counts` was an orphaned module-level
  variable in `core/config/schema.py` instead of a proper field inside
  `ScrubletSettings`.  *Before:* config validation ignored the setting.
  *After:* `on_non_counts` is a validated `ScrubletSettings` field
  defaulting to `"skip_warn"`.

- **Dead fields removed from config templates.** 10 config templates
  slimmed by ~177 lines total.  Fields managed by `global.yaml` defaults
  were removed: `cluster_selection_method`, `umap_selection_method`,
  `param_grid_min_dist`, `param_grid_spread`, `umap_min_dist`,
  `umap_spread`, all DE method params (`method`, `n_genes`,
  `pval_cutoff`, `logfc_cutoff`), `execution.*`, `integration.diagnose*`,
  `integration.gini_*_threshold`, `qc.min_cells_per_gene`.
  `qc.max_pct_mito` replaced with a commented guidance value.
  *Before:* templates duplicated defaults, creating drift risk.
  *After:* single source of truth in `global.yaml` for global defaults.

#### Features

- **`enrich_grid_results()`** — new public function in
  `core/cluster/evaluation.py` that takes a bare
  `results_summary` (n_neighbors, resolution, cluster_key) and enriches
  each entry **in place** with `stability_score`, `cluster_coherence`,
  `kb_annotatable_rate`, and `splitting_gain`.  Also persists the full
  enriched grid at `adata.uns["cluster_grid_results"]`.  Designed as
  the single enrichment entry point consumed by funnel mode, target
  mode, and the standard grid search path.

- **Funnel mode** — progressive grid search for datasets > 100k cells.
  Stratified KMeans++ subsample preserves rare clusters, runs full grid
  on the subsample, ranks by composite score, re-validates top-K
  candidates on full data.  Controlled by `funnel_enabled` (default
  `True` for datasets > `funnel_threshold` = 100,000),
  `funnel_subsample_size` (50,000), `funnel_top_k` (3).  Lineage
  recorded at `adata.uns["funnel_lineage"]`.

- **Target mode** — binary-search resolution to hit a target cluster
  count.  3-run median per candidate guards against Leiden seed
  variance.  Controlled by `target_n_clusters` (default `None` =
  disabled), `target_search_max_iters` (10).  Convergence within ±1
  cluster.

- **KB suggestion metric** — `kb_annotatable_rate` added to the
  5-metric composite scoring.  For each resolution, computes the
  fraction of clusters that can be confidently annotated via the tissue
  knowledge base.  Weight defaults to 0.1 in `multi_metric_weights`.

- **Three-mode dispatch.** `cluster_selection_method` accepts
  `"multi_metric"` (default, 5-metric composite with auto-degrade),
  `"pareto_elbow"` (Pareto frontier + normalized elbow), `"silhouette"`
  (max silhouette), or `None` (manual via `best_resolution` +
  `best_n_neighbors`).  Each mode emits a human-readable reason string
  for the run log.

- **3-tier resolution recommendation.** The multi-metric selector now
  emits coarse / balanced / fine resolution tiers in the log (not stored
  in the h5ad).  Coarse = lowest resolution with composite > 70% of max,
  balanced = argmax composite, fine = highest resolution with stability
  > 85% of max stability.

#### Performance

- **Early termination in `_compute_stability`.** When the last 2
  consecutive pairwise ARIs are both 1.0, the seed loop breaks early
  (clusters are already deterministic).  Saves up to 60% of stability
  compute on well-separated data.

- **`plot_per_combo`** — new flag (default `False`) to skip the
  per-(n_neighbors, resolution) individual UMAP scatter plot loop.
  Summary grid figure, best-params UMAP, and batch-colored UMAP still
  generate regardless.  *Before:* ~2h matplotlib cost at 1M-cell scale
  for cosmetic subplots.  *After:* ~80% plot wall time saved when
  `plot_per_combo: false`.

- **`umap_plot_mode`** — `"auto"` (default) renders full scatter under
  50k cells and subsamples above it; `"full"`, `"subsample"`, and
  `"skip"` for explicit control.  Replaces the hard-coded 100k-cell
  rendering that made large datasets stall for hours.

- **GPU acceleration** (`perf(cluster): 74e9240`).  `_compute_stability`
  routed through `gpu_leiden` (cuGraph, strips flavor/directed/
  n_iterations automatically).  `select_best_umap_params` routed through
  `gpu_neighbors` + `gpu_umap`.  Comparison figure reuses saved UMAP
  coords instead of 3× cold spectral recomputes (~2h+ → seconds on
  Li2026 1M cells).  See also the broader GPU integration below.

- **KNN/UMAP reuse.** Step 04 skips KNN rebuild when `n_neighbors` and
  `use_rep` match the existing graph.  UMAP sweep warm-starts from the
  previous embedding when parameter ranges overlap.  ~6x faster on 1M
  cells.

- **Trustworthiness metric for UMAP selection.** `umap_selection_metric`
  (default `"trustworthiness"`) uses `sklearn.manifold.t_sne.trustworthiness`
  as the quality score; `"convex_hull"` remains available for backward
  compatibility.

#### Config Changes

**Breaking** — fields removed from per-dataset config templates (now in
`global.yaml` defaults):

| Removed Field | Default value |
|---|---|
| `clustering.cluster_selection_method` | `"multi_metric"` |
| `clustering.umap_selection_method` | `"convex_hull"` |
| `clustering.param_grid_min_dist` | `[0.1, 0.3, 0.5]` |
| `clustering.param_grid_spread` | `[1.0]` |
| `clustering.umap_min_dist` | `0.3` |
| `clustering.umap_spread` | `1.0` |
| `de.method` | `"wilcoxon"` |
| `de.n_genes` | `50` |
| `de.pval_cutoff` | `0.05` |
| `de.logfc_cutoff` | `0.25` |
| `execution.*` (n_jobs, random_seed, memory_policy) | global defaults |
| `integration.diagnose` / `diagnose_report` | `true` |
| `integration.gini_batch_threshold` | `0.3` |
| `integration.gini_biology_threshold` | `0.6` |
| `qc.min_cells_per_gene` | `3` |
| `qc.max_pct_mito` | commented guidance |

Existing configs that still set these fields will **not** break —
`ClusteringSettings` uses `extra="ignore"` — but the values will be
ignored in favor of `global.yaml`.  Remove the dead fields to silence
any warnings.

**New opt-in keys:**

| Key | Type | Default | Description |
|---|---|---|---|
| `clustering.funnel_enabled` | `bool` | `True` | Progressive funnel grid search for large data |
| `clustering.funnel_threshold` | `int` | `100000` | Min n_obs to trigger funnel mode |
| `clustering.funnel_subsample_size` | `int` | `50000` | Subsample size for funnel grid |
| `clustering.funnel_top_k` | `int` | `3` | Top-K candidates to re-validate on full data |
| `clustering.target_n_clusters` | `int` or `None` | `None` | Target cluster count (binary search) |
| `clustering.target_search_max_iters` | `int` | `10` | Max binary search iterations |
| `clustering.plot_per_combo` | `bool` | `False` | Skip per-combo UMAP subplots |
| `clustering.umap_plot_mode` | `str` | `"auto"` | UMAP scatter render mode |
| `clustering.umap_plot_max_cells` | `int` | `50000` | Full-render threshold for auto mode |
| `clustering.multi_metric_weights` | `dict` | `{silhouette: 0.2, stability: 0.2, coherence: 0.3, splitting_gain: 0.2, kb_rate: 0.1}` | Composite metric weights |
| `clustering.stability_n_seeds` | `int` | `12` | Seeds for stability evaluation |
| `clustering.multi_metric_de_gate_threshold` | `int` | `25` | Min DE genes for DE-gated selection |
| `execution.device` | `str` | `"cpu"` | `"auto"`, `"cpu"`, or `"gpu"` |

#### GPU Integration Note

The cluster-overhaul builds on the broader GPU integration landed in
`feat(gpu): 237d57d` and `perf(cluster): 74e9240`.  For users with
NVIDIA GPUs and the `[rapids]` extra installed, setting
`execution.device: gpu` activates:

- GPU-accelerated PCA (8.7x on Li2026 1M cells)
- GPU-accelerated Harmony (26.7x)
- GPU-accelerated Leiden grid (20.9x)
- GPU-accelerated UMAP rebuild (7.9x via P0-3 warm-start)
- GPU-routed `_compute_stability` (cuGraph Leiden, ~2h CPU → minutes)
- GPU-routed UMAP sweep + comparison figure reuse (eliminates 3x
  spectral recompute, ~2h+ per cosmetic pass)

All dispatchers fall back to scanpy CPU paths gracefully when RAPIDS is
missing or probing fails.  See `core/utils/_gpu.py` for the full
dispatch layer.
