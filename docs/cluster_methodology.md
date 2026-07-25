# Clustering Methodology

## Overview

Step 04 of the Fuxi pipeline performs unsupervised clustering to partition cells into discrete groups. It wraps a full grid scan over `n_neighbors` x `resolution`, evaluates each combination against five quantitative metrics, and selects the best parameters automatically. The pipeline also optimizes UMAP visualization parameters (min_dist, spread) after the best clustering is found.

### What step 04 does (high level)

1. Builds a nearest-neighbor graph at each `n_neighbors` value in the grid.
2. Runs Leiden clustering at each `resolution` value.
3. Computes 5 evaluation metrics for each (n_neighbors, resolution) pair.
4. Selects the best pair via a weighted composite score (default) or other methods.
5. Optimizes UMAP min_dist / spread for visualization.
6. Saves the final checkpoint with the best clustering labels and UMAP embedding.

### The enrichment pipeline

After raw grid search, each entry is enriched with multi-metric scores:

```
grid_search_clustering(adata, param_grid)
  -> list of {n_neighbors, resolution, n_clusters, silhouette_score}

enrich_grid_results(adata, results, cfg)
  -> adds: stability_score, cluster_coherence, splitting_gain, kb_annotatable_rate
```

The enriched results are then evaluated by `select_best_params()`.

---

## Five Metrics

### 1. Silhouette score (weight: 0.15)

**What it measures:** How well each cell fits its own cluster versus the nearest neighbor cluster. Ranges [-1, 1], higher is better.

**How it's computed:** `sklearn.metrics.silhouette_score` on the PCA/integrated embedding. For datasets larger than `SILHOUETTE_SAMPLE_THRESHOLD`, it's computed on a random 10,000-cell subsample to bound cost.

**Why it's included:** It's the most widely understood clustering quality metric. It provides a baseline that other metrics can refine.

**Limitations:** Tends to prefer compact, spherical clusters. Can be flat across resolutions in subtype-level data where all partitions are similarly mediocre.

### 2. Stability score (weight: 0.20)

**What it measures:** Reproducibility of the clustering across random seeds. Ranges [0, 1], higher is better.

**How it's computed:** Leiden clustering is run `n_seeds` times (default 12, configurable as `stability_n_seeds`) with different random seeds. Pairwise adjusted Rand index (ARI) between all label sets is computed; the mean pairwise ARI is the stability score. Has early termination: if the last 2 consecutive ARIs both equal 1.0, the computation stops.

**Why it's included:** A clustering that changes drastically with seed choice is not trustworthy. Stability guards against spurious partitions.

### 3. Cluster coherence (weight: 0.45), DOMINANT

**What it measures:** Biological interpretability of clusters. Ranges [0, 1], higher is better.

**How it's computed:**

For each cluster, per-cell scores from known marker genes (from `cfg.marker.marker_dict`) are averaged across cell types. A cluster is "coherent" when:

- Its top cell-type marker score exceeds `min_expression` (float or percentile threshold).
- The ratio of the best score to the second-best score exceeds `dominance_threshold` (default 2.5).

Coherence = (number of coherent clusters) / (number of valid clusters). Clusters smaller than `min_cluster_size` are excluded from the denominator.

**Why coherence dominates (weight 0.45):**

Coherence is the only metric that directly measures biological interpretability. A high silhouette score can come from splitting a homogeneous population into arbitrary chunks. A high stability score just means Leiden converged consistently. But coherence tells you: "does each cluster correspond to a single known cell type?"

The weight was set empirically by testing on retina datasets where ground-truth cell types were known. Coherence consistently outperformed other metrics in selecting the resolution that matched manual expert annotation.

**Edge case handling:**

- If all entries have `cluster_coherence < 0.1`, the metric is auto-degraded and the system falls back to silhouette + stability. This catches cases where the marker dict is mismatched to the tissue.
- If no marker dict is provided (`marker_dict` empty) or `adata.raw` is missing, coherence is skipped entirely.

### 4. Splitting gain (weight: 0.15)

**What it measures:** How many new clusters are created per unit resolution increase. Lower is generally better (smooth transitions).

**How it's computed:**

```
splitting_gain(r_i) = max(0, (n_clusters(r_i) - n_clusters(r_{i-1})) / (r_i - r_{i-1}))
```

Computed within each n_neighbors group across the resolution sweep.

**Why it's included:** Prevents the selector from rewarding resolution jumps that create many tiny clusters from a single noisy split. A plateau where increasing resolution adds few clusters is preferred over a cliff where one step creates 10 new clusters.

### 5. KB annotatable rate (weight: 0.05)

**What it measures:** What fraction of clusters can be annotated using the tissue knowledge base. Ranges [0, 1], higher is better.

**How it's computed:** For each cluster, the best-matching cell type (by marker score) is identified. If its mean score exceeds 0.5, the cluster is considered "annotatable." The rate is the fraction of clusters that meet this threshold.

**Why it's included:** A clustering where every cluster matches a known cell type is more interpretable. The low weight (0.05) reflects that this metric is only available when a tissue KB is loaded (`cfg.tissue_kb` set), and that it shouldn't overly influence the result.

### Weight summary

| Metric | Default weight | When available |
|--------|:-------------:|---------------|
| silhouette | 0.15 | Always |
| stability | 0.20 | Always |
| cluster_coherence | 0.45 | Requires marker_dict |
| splitting_gain | 0.15 | Requires >= 2 resolutions per n_neighbors |
| kb_annotatable_rate | 0.05 | Requires tissue_kb + marker_dict |

**Degradation ladder:** If no entry has `cluster_coherence`, the system degrades to 3-metric (silhouette + stability + splitting_gain), then to 2-metric (silhouette + stability), then to silhouette-only if all other metrics fail.

---

## Selection Methods

The selection method is controlled by `cfg.clustering.cluster_selection_method`.

### multi_metric (default)

Computes a weighted composite score for each (n_neighbors, resolution) pair using rank-based normalization. Each metric is independently rank-normalized to [0, 1] via `scipy.stats.rankdata`. The weighted sum of normalized scores is the composite. The pair with the highest composite is selected.

**3-tier resolution recommendation:** The system logs three recommended resolution levels:

| Tier | Selection criterion |
|------|-------------------|
| **coarse** | Lowest resolution with composite > 0.7 of max composite |
| **balanced** | The composite argmax (the actual selection) |
| **fine** | Highest resolution with stability > 0.85 of max stability |

This gives users a quick sense of the resolution landscape without re-running.

### pareto_elbow

Computes the Pareto frontier in (n_clusters, silhouette_score) space, then picks the point closest to the ideal point (k_min=0, s_max=1) on normalized axes. Useful when you want Pareto efficiency between cluster count and silhouette quality.

### silhouette

Simple baseline: pick the (n_neighbors, resolution) pair with the highest silhouette score. Not recommended for production use.

### de_gated

For subtype-level data (detected by `_detect_granularity`). Automatically triggers when the granularity classifier determines the data is "subtype" rather than "tissue." Selects the highest resolution where every cluster has at least `de_gate_threshold` (default 25) differentially expressed genes (one-vs-rest Wilcoxon, padj < 0.05, log2FC > 1.0). Follows the Shekhar 2016 merge.clusters.DE pattern (inverted).

### None (manual)

Uses `best_resolution` and optionally `best_n_neighbors` from config. When `best_n_neighbors=0`, picks the best silhouette at the given resolution. When a user sets both, does exact match.

### Granularity detection

Before selection, `_detect_granularity()` classifies the data as "tissue" (multiple distinct cell types) or "subtype" (FACS-enriched / similar cells). It computes the coefficient of variation (CV) of silhouette scores within the median n_neighbors group. Low CV + low max n_clusters signals subtype data, which triggers DE-gated selection automatically.

---

## UMAP Parameter Selection

After the best (n_neighbors, resolution) is selected, the pipeline sweeps over `param_grid_min_dist` x `param_grid_spread` to find optimal UMAP layout parameters.

Three metrics are available, controlled by `cfg.clustering.umap_selection_metric`:

### trustworthiness (default)

**What it measures:** How well the UMAP embedding preserves local neighborhoods from the high-dimensional space. Ranges [0, 1], higher is better.

**How it's computed:** `sklearn.manifold.t_sne.trustworthiness` between the PCA embedding and the UMAP coordinates with `n_neighbors=15`.

**Why it's the default:** Trustworthiness is a principled dimensionality reduction quality metric. It directly measures whether cells that are neighbors in the original space remain neighbors in the 2D projection. This is what UMAP is supposed to do. Default changed from convex_hull to trustworthiness in v2.1.

### convex_hull (legacy)

**What it measures:** The area of the convex hull of the 2D UMAP embedding. Larger area means better spread of cells.

**How it's computed:** `scipy.spatial.ConvexHull(umap_coords).volume` (area for 2D).

**Why it's legacy:** Convex hull area rewards spreading cells out, but doesn't measure whether the layout preserves the data structure. A layout that scatters cells randomly can have a large convex hull area but terrible local structure. Kept for backward compatibility.

### fixed (skip sweep)

When `umap_selection_metric="fixed"`, the sweep is skipped entirely. The pipeline uses `cfg.clustering.umap_min_dist` and `cfg.clustering.umap_spread` as-is. This saves substantial wall time on large datasets (approximately 3,700s on 1M cells) when you already know good UMAP parameters.

### Performance notes

- The KNN graph is reused from the grid search when `n_neighbors` matches (deduplicates O(n log n) graph building).
- UMAP warm-start: the first (min_dist, spread) combo uses spectral/PAGA initialization; subsequent combos warm-start from the previous embedding. This saves approximately 40 min per combo on 1M-cell datasets.
- Sweep results (including UMAP coordinates) are saved to `sweep_results` so the comparison figure can render without re-running UMAP.

---

## Three-Mode Grid Scan Dispatch

The pipeline dispatches to one of three grid-scan modes based on dataset size and config:

### 1. Target mode

**Trigger:** `cfg.clustering.target_n_clusters` is set to a value (e.g., 20).

**Behavior:**
- For each value in `param_grid_n_neighbors`, binary-search resolution to hit the target cluster count.
- Resolution is searched in [0.1, 5.0] with 3-run median at each step (seeds 42, 123, 456) to guard against Leiden non-monotonicity.
- Converges when |actual_k - target_k| <= 1.
- Results are enriched with multi-metric scores.
- Best used when you have a strong prior about the number of cell types.

**Example config:**
```yaml
clustering:
  target_n_clusters: 20
  param_grid_n_neighbors: [15, 20, 30]
```

### 2. Funnel mode

**Trigger:** `funnel_enabled=True` (default) AND `adata.n_obs > funnel_threshold` (default 100,000). Also requires `target_n_clusters` to be None.

**Behavior:**
1. **Stratified subsample** the full dataset to `funnel_subsample_size` (default 50,000) via KMeans++ per batch. Preserves rare clusters (forces >= 5 cells per cell type if available).
2. Run the full grid search on the subsample.
3. Rank candidates by composite score, keep top-K (default 3).
4. Re-validate top-K candidates on the full dataset (Leiden + silhouette + enrichment).
5. Pick the best from full validation.

**Why it matters:** Full grid search on 1M+ cells would run Leiden 18-30 times x 5 stability seeds = 90-150 Leiden runs. Funnel reduces this by testing the parameter landscape on a subsample and only validating promising candidates at full scale.

**Funnel lineage:** The subsample index and composite score changes are stored at `adata.uns["funnel_lineage"]` for debugging.

### 3. Full-grid mode

**Trigger:** Default. When target_n_clusters is None AND (funnel is disabled OR n_obs <= threshold).

**Behavior:**
- Exhaustive grid search over all (n_neighbors, resolution) pairs.
- Each n_neighbors group rebuilds the KNN graph once, then runs Leiden at all resolutions.
- Results are enriched and evaluated.

### Dispatch decision tree

```
if target_n_clusters is set:
  -> TARGET mode (binary-search resolution)

elif funnel_enabled AND n_obs > funnel_threshold:
  -> FUNNEL mode (subsample -> rank -> re-validate)

else:
  -> FULL GRID mode (exhaustive)
```

---

## KB Suggestion: Strategy C

### What it does

When `target_n_clusters` is not set in the config but `cfg.tissue` is provided, the pipeline queries the tissue knowledge base for a suggested cluster count. The suggestion is logged as an informational message (not automatically applied).

### Strategy C: species-priority with all-source median fallback

The KB maintains a database of cluster counts from published studies across tissues and species, stored in per-tissue `sources/*.yaml` files under `core/kb/<tissue>/`.

The algorithm (`suggest_target_n_clusters` in `core/kb/__init__.py`):

1. **Species-filtered match:** If species is provided, filter entries matching both tissue AND species.
   - >= 2 matches: return median cluster count.
   - 1 match: log a warning (single source), return that count.
   - 0 matches: fall through to step 2.

2. **All-source median:** Return the median of all entries for this tissue (any species).
   - 0 matches: return None (no suggestion).
   - 1 match: log a warning (single source).

**Why species filtering first:** Cluster counts vary substantially by species for the same tissue (e.g., human retina ~16 subtypes, mouse retina ~39 subtypes). Species-priority ensures the suggestion is relevant. The all-source median fallback provides a reasonable guess when the specific species has no published data.

### Example output

```
KB suggests target_n_clusters=39 for tissue=retina species=mouse.
Set clustering.target_n_clusters=39 in config to activate target mode.
```

---

## Adaptive n_neighbors

### Per-dataset-size scaling

The config field `param_grid_n_neighbors_adaptive: bool = True` enables automatic scaling of the n_neighbors grid based on dataset size. When enabled, the grid values are computed proportionally to `adata.n_obs`.

**Key principle:** Larger datasets require more neighbors to capture the same local structure. A fixed grid of [15, 20, 30] that works for 10,000 cells would be too dense for 1,000,000 cells.

**Scaling logic (planned):** The adaptive scaling uses a square-root heuristic relative to a reference dataset size. When `param_grid_n_neighbors_adaptive` is False, the explicit `param_grid_n_neighbors` list is used as-is.

### Adaptive resolution expansion for small datasets

When `multi_metric_adaptive_resolution=True` (default) and `adata.n_obs < 3000`, the pipeline automatically extends the resolution grid to include higher resolutions (3.0, 5.0) that would otherwise be omitted. This prevents under-clustering on small datasets where the default resolution range [0.3, 2.0] may not produce enough clusters.

---

## Lab Server vs GPU

The clustering step can run on CPU or GPU. The device is selected via `cfg.execution.device`.

| Setting | Behavior |
|---------|----------|
| `"auto"` | Use GPU if RAPIDS is available, else CPU |
| `"cpu"` | Force CPU execution |
| `"gpu"` | Force GPU execution (raises error if unavailable) |

### What gets accelerated on GPU

| Operation | CPU library | GPU library | Typical speedup |
|-----------|-------------|-------------|:---------------:|
| KNN graph | `scanpy.pp.neighbors` | `rapids_singlecell` | 2-5x |
| Leiden clustering | `scanpy.tl.leiden` (igraph) | `cuGraph` | 20x (18 combos: 8520s -> 408s) |
| UMAP | `scanpy.tl.umap` (umap-learn) | `cuml.UMAP` | 1.3-8x (varies with warm-start) |
| PCA | `scanpy.pp.pca` | `cuml.PCA` | 8.7x |
| Harmony | `harmonypy` | `rapids_singlecell` | 26.7x |

### CPU-only tasks

Some parts of the enrichment pipeline remain CPU-bound and see no GPU acceleration:

- Multi-metric enrichment (stability Leiden on GPU was added in commit 74e9240).
- DE gene testing (Wilcoxon rank-sum).
- Matplotlib figure generation.
- Marker score computation.

### Configuration for each environment

**Lab server (no GPU):**
```yaml
execution:
  device: cpu
clustering:
  funnel_threshold: 100000  # or lower for slower CPUs
```

**Workstation with NVIDIA GPU:**
```yaml
execution:
  device: auto  # or gpu
clustering:
  funnel_threshold: 100000
```

**Detailed GPU performance notes** are in [`notes/features/2026-07-25_gpu_rapids_integration.md`](../notes/features/2026-07-25_gpu_rapids_integration.md). Key findings from the Li2026 1.1M-cell benchmark:

- Step 04 total wall: CPU ~10h -> GPU ~3h (3.3x).
- The dominant gain was in the Leiden grid search (20.9x) and final UMAP rebuild (7.9x).
- UMAP warm-start (P0-3) was the most impactful CPU optimization, reducing per-combo UMAP from 47min to 7min.
- The comparison figure fix (commit 74e9240) saved approximately 2.4h by reusing sweep coordinates instead of recomputing UMAP.

---

## Config Reference

### Clustering settings (schema default values)

```yaml
clustering:
  # Grid
  param_grid_n_neighbors: [15, 20, 30]
  param_grid_n_neighbors_adaptive: true
  param_grid_resolutions: [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]

  # Selection
  cluster_selection_method: multi_metric
  multi_metric_weights:
    silhouette: 0.2
    stability: 0.2
    cluster_coherence: 0.3
    splitting_gain: 0.2
    kb_annotatable_rate: 0.1
  stability_n_seeds: 12
  multi_metric_adaptive_resolution: true
  multi_metric_coverage_ratio_threshold: 1.5
  multi_metric_coherence_dominance: 1.5
  multi_metric_granularity_cv_threshold: 0.05
  multi_metric_granularity_min_clusters: 10
  multi_metric_de_gate_threshold: 25

  # Target mode
  target_n_clusters: null
  target_search_max_iters: 10

  # Funnel mode
  funnel_enabled: true
  funnel_threshold: 100000
  funnel_subsample_size: 50000
  funnel_top_k: 3

  # UMAP
  umap_selection_method: convex_hull
  umap_selection_metric: trustworthiness
  param_grid_min_dist: [0.1, 0.3, 0.5]
  param_grid_spread: [1.0]
  umap_min_dist: 0.3
  umap_spread: 1.0
  umap_plot_mode: auto
  umap_plot_max_cells: 50000
  plot_per_combo: false
```

### Execution settings

```yaml
execution:
  device: auto         # auto | cpu | gpu
  n_jobs: -1
  random_seed: 42
```

---

## Example Workflow

### Default workflow (tissue-level, 50k cells)

```yaml
# config.yaml
clustering:
  cluster_selection_method: multi_metric
```

Pipeline output:

```
[Grid search] 3 n_neighbors x 6 resolutions = 18 combinations
[Enrichment]  computing stability, coherence, splitting_gain, kb_rate...
[Selection]  composite=0.7814 sil=0.4231(n=0.67) stab=0.923(n=0.83)
             coherence=0.889 split_gain=0.167 kb_rate=0.933 k=16
             weights=sil:0.13,stab:0.17,coh:0.43,split:0.17,kb:0.10
             3-tier: coarse: r=0.30 (k=8) / balanced: r=0.80 (k=16)
             / fine: r=1.50 (k=24)
```

### Target mode (known cluster count)

```yaml
clustering:
  target_n_clusters: 20
  param_grid_n_neighbors: [15, 20, 25, 30]
```

### Funnel mode (1M cells, auto-triggered)

No config change needed. The funnel will trigger automatically when `n_obs > 100,000`.

### Manual override

```yaml
clustering:
  cluster_selection_method: null
  best_resolution: 1.0
  best_n_neighbors: 20
```
