# Fuxi (伏羲) Architecture

## Overview

Fuxi is a modular, multi-omics single-cell analysis pipeline that unifies
scRNA-seq (Scanpy), scATAC-seq (Snapatac2), and spatial transcriptomics
(Squidpy) under a shared core infrastructure.  The pipeline follows a
**layered decision architecture** in which each layer contributes
increasingly confident evidence toward a final analytical result.

---

## 1. Evidence Fusion Pattern

The **Evidence Fusion Pattern** is Fuxi's canonical design for combining
multiple sources of biological evidence into a single high-confidence
decision.  It is a 5-layer decision engine with strict priority ordering.

### Canonical implementation

| Module | Role |
|--------|------|
| ``rna/utils/marker_scoring.py`` | Layer 1–2: statistical + similarity scoring |
| ``rna/utils/evidence_fusion.py`` | Layer 3–5: rule engine + fusion + AI arbitration |

### The 5 layers (evaluated in priority order)

```
Layer 0 — Expert rules         [hard-coded biological knowledge]
Layer 1 — High-confidence      [score >= 0.7]
Layer 2 — Medium-confidence    [0.5 <= score < 0.7]
Layer 3 — Low-confidence        [0.25 <= score < 0.5]
Layer 4 — Unknown / fallback   [all else]
```

#### Layer 0 — Expert Rules
Priority: **Highest**.  If the expert-rule engine matches a cluster
(e.g. "these 5 marker genes are present"), the rule's cell type is
accepted immediately, subject to a quality gate: if independent marker
scoring finds zero KB overlap, confidence is downgraded from ``rule``
to ``low``.

```python
# core/evidence_fusion.py — Tier 0
if expert_rule_result is not None:
    rule_score, rule_n = _resolve_score(marker_scores, expert_rule_result)
    if rule_score < 0.25 and rule_n == 0:
        conf = 'low'  # downgrade: rule triggered by noise
    else:
        conf = 'rule'
    return FusionDecision(expert_rule_result, conf, ...)
```

#### Layers 1–3 — Marker Scoring
When no expert rule fires, the pipeline falls through to marker-scoring
evidence.  The scoring engine uses:

1. **Hypergeometric (Fisher's exact) test** — enrichment of KB marker
   genes in the cluster's top-N DE genes, optionally with
   consensus-weighted counts (gold markers count 3x, high 2x, etc.).
2. **Cosine similarity** — vector-space overlap between cluster DE genes
   and the KB marker set.
3. **Confidence multiplier** — cell types with <= 2 markers get 0.5x,
   3–5 markers get 0.8x, > 5 markers get 1.0x.
4. **Phylogenetic weight** — when ``target_class`` is set, scores are
   multiplied by a taxonomic-distance weight (same class + same order = 1.0,
   same class + different order = 0.8, different class = 0.6–0.9).

```python
# rna/utils/marker_scoring.py — score construction
final_score = hypergeometric_score * conf_mult
if target_class:
    p_weight = phylogenetic_weight(source_cls, target_class, ...)
    final_score *= p_weight
if neg_penalty:
    final_score *= 0.5
```

#### Layer 4 — Unknown / Diagnostic
When no layer produces a confident result, the cluster is classified as
``Unknown`` with a diagnostic category explaining why:

| Category | Meaning |
|----------|---------|
| ``no_kb_match`` | No KB cell type had any marker overlap |
| ``low_quality_data`` | Cluster failed low-quality detection |
| ``ambiguous`` | Multiple cell types with score >= 0.25 |
| ``weak_signal`` | Best score > 0 but < 0.25 |
| ``true_unknown`` | Could not determine by any method |

### Data flow

```
  Cluster DE genes
         │
         ▼
  ┌──────────────────┐
  │  marker_scoring   │  Fisher exact + cosine sim + phylogenetic weight
  │  (statistical)    │  → {type_key: Score}
  └────────┬─────────┘
           │
           ▼
  ┌──────────────────┐
  │  apply_expert    │  Hard-coded biological rules
  │  _rules          │  → (cell_type or None, [alternative_rules])
  └────────┬─────────┘
           │
           ▼
  ┌──────────────────┐
  │  fuse_evidence    │  Tiered decision (0 → 1 → 2 → 3 → 4)
  │  (orchestrator)   │  → FusionDecision
  └────────┬─────────┘
           │
           ▼
  ┌──────────────────┐
  │  AI arbitration   │  Optional: compare with LLM suggestion
  │  (optional)       │  → ai_agreed flag in FusionDecision
  └──────────────────┘
```

### Caller

The fusion engine is invoked from ``rna/steps/05_annotate_major.py``
(and ``rna/annotation_engine.py``) via:

```python
decision = fuse_evidence(
    marker_scores=marker_scores,
    expert_rule_result=rule_result,
    ai_suggestion=ai_results.get(cluster),
    unconstrained=CFG.ai.unconstrained_annotation,
)
```

### Quality metadata

``fuse_all_clusters()`` can return quality metadata on request:

```python
decisions, quality = fuse_all_clusters(
    all_scores, all_rules, ..., return_quality=True
)
# quality = {
#     "annotated_by_rule": 3,
#     "annotated_by_scoring": 18,
#     "unknown": 1,
#     "ambiguity": 0,
#     "ai_agreed": 21,
#     "total": 22,
#     "diagnostic_summary": {"no_kb_match": 1},
# }
```

### Design principles

1. **Strict priority** — higher layers never override lower ones; they
   can only downgrade confidence (quality gate for rules).
2. **Diagnostic transparency** — every ``Unknown`` decision carries a
   machine-readable ``DiagnosticInfo`` explaining why, enabling automated
   downstream audits.
3. **AI as arbiter, not oracle** — the LLM is consulted only at medium/low
   confidence tiers and records agreement; it never overrides the
   statistical result (unless ``unconstrained_annotation`` is active).
4. **Extensibility** — ``DECISION_TIERS`` is a simple list of
   ``(name, predicate)`` tuples; new tiers can be inserted without
   changing the core fusion loop.

---

## 2. Module structure

```
fuxi/
├── core/                    # Shared infrastructure
│   ├── config.py            # Unified Config + nested modality configs
│   ├── utils.py             # safe_write, safe_plot, resolve_config, ...
│   ├── ai_caller.py         # LLM calls with retry + caching
│   ├── ai_prompts.py        # Annotation / interpretation templates
│   ├── run_pipeline.py      # CLI dispatcher (--modality rna|atac|spatial)
│   ├── dataset_schema.py    # dataset.yaml Python model
│   ├── dataset_detector.py  # Auto-detect modality from file patterns
│   ├── path_validation.py   # Safe path traversal guards
│   ├── clustering.py        # Shared grid-search clustering interface
│   │   ├── anatomy.py          # Anatomical adjacency loading & CCI filtering
│   └── preprocess/          # Preprocessing pipeline
│       ├── preprocessor.py      # Orchestrator
│       ├── metadata_parser.py   # Phase 4: dataset.yaml generation
│       ├── matrix_loader.py     # Phase 5: config generation
│       ├── format_detector.py   # File format classification
│       ├── archive_extractor.py # Archive extraction (tar.gz, zip, ...)
│       └── superseries_detector.py
├── rna/
│   ├── steps/               # 12 pipeline steps (scanpy)
│   ├── utils/
│   │   ├── evidence_fusion.py    # 5-layer decision engine
│   │   ├── marker_scoring.py     # Statistical scoring + expert rules
│   │   └── cluster_evaluation.py # Silhouette / Pareto-elbow
│   ├── annotation_engine.py # Cross-module annotation API
│   └── ortholog.py          # Gene name conversion (Ensembl ↔ symbol)
├── atac/
│   └── steps/               # 10 pipeline steps (snapatac2)
├── spatial/
│   └── steps/               # 11 pipeline steps (squidpy)
├── projects/                # Dataset-specific configs
├── tests/                   # Unified test suite
└── templates/               # Config templates
```

## 3. Config system

The ``Config`` dataclass (``core/config.py``) uses **nested modality
configs** for clean separation:

```python
cfg = Config()
cfg.modality = "rna"
cfg.rna.n_top_genes = 4000       # RNA-specific
cfg.atac.min_fragments = 1000    # ATAC-specific
cfg.spatial.library_id = ""      # Spatial-specific
```

Backward-compatible access via ``__getattr__``:

```python
cfg.n_top_genes  # → resolves to cfg.rna.n_top_genes when modality == "rna"
```

Common fields (directly on ``Config``):
- ``data_dir``, ``results_dir``, ``h5ad_dir``, etc.
- ``tissue``, ``species``, ``expression_type``
- ``ai`` (``AIConfig``)
- ``n_jobs``, ``random_seed``, ``h5ad_compression``, etc.

## 4. Preprocessing pipeline

The preprocessing pipeline (``core/preprocess/``) phases:

| Phase | Module | Output |
|-------|--------|--------|
| 0 | ``preprocessor.py`` | Validate input directory |
| 1 | ``archive_extractor.py`` | Extract tar.gz / zip / gz / bz2 |
| 2 | ``superseries_detector.py`` | Detect SuperSeries structure |
| 3 | ``format_detector.py`` | Classify files by format + infer modality |
| 4 | ``metadata_parser.py`` | Generate ``dataset.yaml`` |
| 5 | ``matrix_loader.py`` | Generate ``config_GSE_ID.py`` |
| 6 | ``preprocessor.py`` | Summary report |
