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

**Config access.** Use nested topic paths: `CFG.hvg.n_top_genes`, `CFG.clustering.cluster_selection_method`. `.py` configs are rejected — use `.yaml`.

**Core scripts.** Step scripts under `rna/steps/`, `atac/steps/`, `spatial/steps/` must not be edited in place. Copy to `projects/{modality}/{GSE_ID}/` first.

## Running methods

### Running modes

Both modes use `--step N` to run one step at a time. The difference is whether the Agent pauses for user input:

**Auto mode** — Execute steps sequentially with default settings. After each step, report progress and continue. User can interrupt at any point (Ctrl+C). Suitable for familiar datasets or batch reproduction.

**Interactive mode** — Execute `--step N` one at a time. After each step, present results, ask questions, offer options, and wait for confirmation before proceeding. Suitable for exploratory analysis or new datasets.


```bash
# List steps
python core/run_pipeline.py --modality rna --list
python core/run_pipeline.py --modality atac --list

# Run full pipeline / single step / resume
python core/run_pipeline.py --modality rna --config projects/rna/<GSE_ID>/config_<GSE_ID>.yaml
python core/run_pipeline.py --modality atac --step 0 --config ...
python core/run_pipeline.py --modality rna --resume --config ...

# Paper tools
python core/paper_insights.py --pmid <PMID>       # AI paper interpretation
python core/paper_registry.py --build              # build paper→GSE→config index
python core/paper_registry.py --verify             # check registry consistency
python core/run_reproduce.py --all --dry-run       # preview reproducibility
python core/run_reproduce.py <paper_dir>           # reproduce a single paper
```

### Key paths

| Module | Location |
|--------|----------|
| RNA steps | `rna/steps/` (13 steps: 00_load → 12_cell_interaction) |
| ATAC steps | `atac/steps/` (10 steps: 00_load → 09_integrate) |
| Spatial steps | `spatial/steps/` (11 steps: 00_load → 10_cell_interaction) |
| Shared core | `core/` (config, utils, ai_caller, preprocess) |
| Paper tools | `core/paper_insights.py`, `core/paper_registry.py`, `core/run_reproduce.py` |
| Project configs | `projects/{modality}/{GSE_ID}/config_*.yaml` |
| Config templates | `templates/config_templates/*.yaml` |

### Dataset lookup & unified registry

When user requests to analyze or reproduce a dataset, check the **master registry** first:

1. **Quick summary** — `python -m core.registry report`
   (counts papers, datasets, links, orphans).
2. **Consistency check** — `python -m core.registry verify`
   (detects stale config paths, missing dirs, orphan datasets).
3. **Find orphans** — `python -m core.registry find-orphans`
   (datasets or supplement dirs with no linked paper).
4. **Programmatic query** — python import:
   ```python
   from core.registry import load_master_registry
   reg = load_master_registry()
   reg.get_dataset_links("41578023")   # paper → datasets
   reg.get_paper_links("GSE118614")    # dataset → papers
   reg.get_paper(paper_id="41578023")  # by PMID/slug
   reg.find_orphans()
   ```
5. **If paper found** — Read `projects/papers/<paper_dir>/insights.yaml` for metadata.
   Cross-validate with NCBI: `python core/paper_insights.py --pmid <PMID>`.
6. **If not found** — Ask user whether to download and register:
   ```bash
   python core/paper_insights.py --pmid <PMID>
   python core/paper_registry.py --build
   ```


### Critical conventions

- Use `.venv/bin/python` for all Python commands
- Source `.env` first: `set -a && source .env && set +a`
- Steps run as **subprocesses** — never import step scripts directly
- Every step must prepend repo root to `sys.path`
- Config: `.yaml` + `resolve_config()` → Pydantic v2; `.py` rejected
- `data_root()` requires `FUXI_DATA_ROOT` env var
- Import pattern: `from core.utils import ...`, `from core.ai_caller import ...`
