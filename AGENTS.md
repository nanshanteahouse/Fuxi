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

Pipeline supports two modes, chosen by the Agent based on user preference:

**Auto mode** — Run full step range with all default settings. Suitable for familiar datasets or batch reproduction:
```bash
python core/run_pipeline.py --modality rna --resume --config <config>.yaml
```

**Interactive mode** — Agent runs `--step N` one at a time. After each step, present results to the user, ask questions, offer options, and wait for confirmation before proceeding. Suitable for exploratory analysis or new datasets.


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

### Paper workflow

1. **Check registry** — `python core/paper_registry.py --verify` or grep `projects/papers/registry.yaml`
2. **Read insights** — `projects/papers/<paper_dir>/insights.yaml`
3. **Cross-validate** — `python core/paper_insights.py --pmid <PMID>`


### Critical conventions

- Use `.venv/bin/python` for all Python commands
- Source `.env` first: `set -a && source .env && set +a`
- Steps run as **subprocesses** — never import step scripts directly
- Every step must prepend repo root to `sys.path`
- Config: `.yaml` + `resolve_config()` → Pydantic v2; `.py` rejected
- `data_root()` requires `FUXI_DATA_ROOT` env var
- Import pattern: `from core.utils import ...`, `from core.ai_caller import ...`
