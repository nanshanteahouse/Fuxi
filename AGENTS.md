# Fuxi (伏羲)

Unified single-cell multi-omics pipeline — scRNA-seq (Scanpy), scATAC-seq (Snapatac2), Spatial (Squidpy). Python 3.14+ on WSL2.

## Development conventions

**Commit message.** Use [Conventional Commits](https://www.conventionalcommits.org/):
```
<type>(<scope>): <subject>
```
Types: `feat` / `fix` / `docs` / `refactor` / `perf` / `test` / `chore`
Scope: semantic name (e.g. `enrichment`), NOT step number (e.g. `09_enrichment`)
Subject: imperative, lowercase, ≤72 chars. Body explains *why*, not *what*.

**Commit scope.** Before each commit, review what changed in this session (`git diff --stat`, `git status`). Unless explicitly instructed otherwise, only stage and commit files that belong to *this session's work*. Do not revert, delete, or touch unrelated files.

**Config access.** Use nested topic paths: `CFG.hvg.n_top_genes`, `CFG.clustering.cluster_selection_method`. `.py` configs are rejected — use `.yaml`.

**Core scripts.** Step scripts under `rna/steps/`, `atac/steps/`, `spatial/steps/`, `bulk/steps/` must not be edited in place. Copy to `projects/{modality}/{GSE_ID}/` first.
**Ad-hoc scripts.** One-off / dataset-specific analysis scripts under `adhoc/`. Not part of the pipeline, no compatibility guarantee — use once and discard.

**Code organization.** Soft caps: 500 LOC for core modules and step scripts, 400 for utility modules. Algorithm engines under `*_utils/` at 500 LOC.
## Running methods

### Running modes

Both modes use `--step N` to run one step at a time. The difference is whether the Agent pauses for user input:

**Auto mode** — Execute steps sequentially with default settings. After each step, report progress and continue. User can interrupt at any point (Ctrl+C). Suitable for familiar datasets or batch reproduction.

**Interactive mode** — Execute `--step N` one at a time. After each step, present results, ask questions, offer options, and wait for confirmation before proceeding. Suitable for exploratory analysis or new datasets.
**TUI mode** — Launch the unified terminal interface: `python -m core.tui`. Keyboard-navigable dashboard with registry browser, pipeline runner, config editor, and results viewer. Ideal for project exploration and batch management.
**MCP mode** — Start an AI-accessible server that lets LLM agents (Claude Desktop, VS Code Copilot, custom agents) query the registry, check pipeline status, and trigger downloads/preprocessing/pipeline runs through the Model Context Protocol. See `docs/mcp_guide_zh-CN.md` for setup instructions.




```bash
# List steps
python core/run_pipeline.py --modality rna --list
python core/run_pipeline.py --modality atac --list
python core/run_pipeline.py --modality bulk --list

# Run full pipeline / single step / resume
python core/run_pipeline.py --modality rna --config projects/rna/<GSE_ID>/config_<GSE_ID>.yaml
python core/run_pipeline.py --modality atac --step 0 --config ...
python core/run_pipeline.py --modality rna --resume --config ...
# Spatial
python core/run_pipeline.py --modality spatial --step 0 --config ...

# Bulk RNA-seq
python core/run_pipeline.py --modality bulk --list
python core/run_pipeline.py --modality bulk --config projects/bulk/<GSE_ID>/config_<GSE_ID>.yaml
python core/run_pipeline.py --modality bulk --step 2 --config ...

# Step range / subclustering / subset
python core/run_pipeline.py --modality rna --steps 0-2 --config ...
python core/run_pipeline.py --modality rna --step 6 --cell-type "Müller Glia" --config ...
# Subset: filter cells by sample/obs criteria (Step 00)
python core/run_pipeline.py --modality rna --config config_pcw8.yaml
# Config: CFG.downsample.sample_keep=["SCR205"] or CFG.downsample.obs_filter="stage=='PCW8'"

# Paper tools
python core/paper/insights.py --pmid <PMID>       # AI paper interpretation
python core/paper/insights.py --pmid <PMID> --methodology  # + methodology patterns
python adhoc/migration_scripts/methodology_batch.py          # batch methodology backfill
python -m core.paper.registry report               # print summary
python -m core.paper.registry verify              # check registry consistency
python -m core.paper.registry register --gse GSE123456  # register GSE → PMID link
python -m core.paper.registry find-orphans        # list orphan datasets
python core/pipeline/reproduce.py --all --dry-run       # preview reproducibility
python core/pipeline/reproduce.py <paper_dir>           # reproduce a single paper

# MCP server
python -m core.ai.mcp_server                         # stdio mode (for AI agents)
python -m core.ai.mcp_server --http 8080             # HTTP mode (for remote clients)
```

### Adding a dataset

**Automated (recommended):**
```bash
python core/preprocess/preprocessor.py --gse <GSE_ID> --data-root $FUXI_DATA_ROOT --download
```
**Manual:** Generate a starter config from the format specs: `python -m core.config scaffold --list` (see format keys), then `python -m core.config scaffold --format <KEY> --out projects/{modality}/{GSE_ID}/config_<GSE_ID>.yaml`. Edit the file afterwards — do not edit `core/preprocess/config_specs.py` and re-render by hand; the specs are the source of truth.

**Schema → config workflow.** `core/config/schema.py` is the single source of truth for config fields/defaults; starter configs are rendered on demand from `core/preprocess/config_specs.py` specs (no committed template files). When adding a schema field: optionally add it to the relevant spec(s) in `config_specs.py`; `validate_specs()` fails CI on any dotted path that is not a real schema field.

### Key paths

| Module | Location |
|--------|----------|
| RNA steps | `rna/steps/` (13 steps: 00_load → 12_cell_interaction) |
| ATAC steps | `atac/steps/` (13 steps: 00_load → 12_integrate) |
| Spatial steps | `spatial/steps/` (11 steps: 00_load → 10_cell_interaction) |
| Bulk steps | `bulk/steps/` (5+1 steps: 00_load → 04_exploratory, optional 05_batch) | PyDESeq2 |
| Shared core | `core/` — sub-packages: ai/, annotation/, cluster/, config/, interaction/, kb/, paper/, pipeline/, preprocess/, utils/ |
| Paper tools | `core/paper/` (insights.py, registry.py, converter.py, cross_paper.py) |
| Methodology tools | `core/paper/insights.py --methodology`, `adhoc/migration_scripts/methodology_batch.py` |
| Knowledge base | `core/kb/` (tissue ontologies, marker validation, adjacency) |
| Ad-hoc scripts | `adhoc/` (one-off migration, ortholog processing, dataset-specific analysis) |
| Brainstorming | `projects/notebook/` (methodology_ideas, keywords, etc.) |
| TUI | `core/tui/` (7 backends, 6 screens, 4 widgets) — `python -m core.tui` |
| MCP server | `core/ai/mcp_server.py` (10 tools: registry + pipeline + execution) — `python -m core.ai.mcp_server` |

### Key design patterns

**Step dispatch.** `core/run_pipeline.py` runs each step as a separate `subprocess.run()` — never import steps directly. Steps self-identify checkpoint files; skip if output exists. `--resume` scans for first missing checkpoint.

**Config loading.** `resolve_config()` loads `.yaml` via `yaml.safe_load()` + `Config.model_validate()`. Each call returns a new Config instance. Path resolution via `model_post_init` hook.

**Three annotation modes (RNA Step 05):**
1. **Unified KB** (if `CFG.tissue_kb` set): marker scoring → evidence fusion → optional AI fallback
2. **AI** (if `CFG.ai.enabled`): LLM with StandardOntology normalization
3. **Score_genes** (fallback): `sc.tl.score_genes()` with `CFG.marker.marker_dict`

**Cluster selection (Step 04).** Controlled by `CFG.clustering.cluster_selection_method`:
- `multi_metric` (RNA default, MMACS v2): 5-metric composite, tissue or subtype path
- `pareto_elbow` (ATAC/Spatial default): Pareto frontier over n_clusters vs silhouette
- `silhouette` / `None` (manual)
UMAP sweep auto-selects `min_dist`/`spread` via convex hull area.

**snRNA-seq adaptation.** Auto-detected from GEO keywords → `CFG.qc.is_nuclei=True`. Tightens mito threshold to `max_pct_mito_nuclei` (default 3.0%) and MAD multiplier to 1.5×.

**Subset filtering.** Step 00 supports `CFG.downsample.sample_keep` and `CFG.downsample.obs_filter`. Output dirs auto-append `_subset`.

### Dataset & Paper lookup

When user asks about a dataset (GSE), paper (PMID), or reproduction:
use `python -m core.paper.registry status --gse <ID>` or `--pmid <ID>`.
Read `docs/agent/registry_lookup.md` for the full decision tree, fuzzy search,
registration commands, and programmatic API reference.
### Methodology Pattern Analysis

Extract methodological fingerprints from papers (5 dimensions: archetype,
strategy, narrative, toolbox, contribution).
See `docs/agent/methodology_patterns.md` for commands and schema reference.

### Critical conventions

- Use `.venv/bin/python` for all Python commands
- Source `.env` first: `set -a && source .env && set +a`
- Steps run as **subprocesses** — never import step scripts directly
- Every step must prepend repo root to `sys.path`
- Config: `.yaml` + `resolve_config()` → Pydantic v2; `.py` rejected
- `data_root()` requires `FUXI_DATA_ROOT` env var
- Import pattern: `from core.utils import ...`, `from core.ai.caller import ...`

## Session report generation

When user explicitly asks for a report (写报告, 写总结, write a report, etc.),
generate a markdown report under `notes/`. Read `docs/agent/report_generation.md`
for directory mapping, naming rules, and edge cases.
Mandatory conventions (see `docs/agent/report_generation.md` → Notes 库维护惯例 section):
- Every new report carries a 状态头 (status header): `> 状态：现行` + `> 后继/关联` lines (full-width colon).
- When a report completes a backlog item, check it off in `notes/BACKLOG.md` (set done + 备注 link) in the same session.
- Write a two-week work log under `notes/logs/`; move superseded reports to `notes/archive/` (never delete) and update `notes/INDEX.md` in the same session.


## Notes are private — never commit

The entire `notes/` directory is **gitignored**. Notes contain internal
thinking, incomplete analysis, and experimental reports — these are for
internal reference only and must **never** appear in the public repository.

- **NEVER** run `git add notes/` or `git add -f notes/` under any circumstance.
- **NEVER** include notes/ files in commits, Pull Requests, or pushes.
- Before every commit, verify with `git status` that no notes/ paths are staged.
- If a commit hook or `git add -A` accidentally stages notes/ files, unstage
  them before committing.

## Pre-push hooks & SSH keep-alive

Pre-push: pyright type check + smoke test (~50-80s on WSL2).
If `git push` fails with SSH disconnect, see `docs/agent/ssh_troubleshooting.md`.
