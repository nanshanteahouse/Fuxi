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

**Core scripts.** Step scripts under `rna/steps/`, `atac/steps/`, `spatial/steps/` must not be edited in place. Copy to `projects/{modality}/{GSE_ID}/` first.
**Ad-hoc scripts.** One-off / dataset-specific analysis scripts under `adhoc/`. Not part of the pipeline, no compatibility guarantee — use once and discard.

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
python core/paper/insights.py --pmid <PMID>       # AI paper interpretation
python core/paper/insights.py --pmid <PMID> --methodology  # + methodology patterns
python adhoc/migration_scripts/methodology_batch.py          # batch methodology backfill
python -m core.paper.registry report               # print summary
python -m core.paper.registry verify              # check registry consistency
python -m core.paper.registry register --gse GSE123456  # register GSE → PMID link
python -m core.paper.registry find-orphans        # list orphan datasets
python core/pipeline/reproduce.py --all --dry-run       # preview reproducibility
python core/pipeline/reproduce.py <paper_dir>           # reproduce a single paper
```

### Key paths

| Module | Location |
|--------|----------|
| RNA steps | `rna/steps/` (13 steps: 00_load → 12_cell_interaction) |
| ATAC steps | `atac/steps/` (14 steps: 00_load → 13_integrate) |
| Spatial steps | `spatial/steps/` (11 steps: 00_load → 10_cell_interaction) |
| Shared core | `core/` — sub-packages: ai/, annotation/, cluster/, config/, interaction/, kb/, paper/, pipeline/, preprocess/, utils/ |
| Paper tools | `core/paper/` (insights.py, registry.py, converter.py, cross_paper.py) |
| Methodology tools | `core/paper/insights.py --methodology`, `adhoc/migration_scripts/methodology_batch.py` |
| Knowledge base | `core/kb/` (tissue ontologies, marker validation, adjacency) |
| Ad-hoc scripts | `adhoc/` (one-off migration, ortholog processing, dataset-specific analysis) |
| Brainstorming | `projects/notebook/` (methodology_ideas, keywords, etc.) |

### Dataset & Paper lookup

When user requests to analyze a dataset, look up a paper, or reproduce results,
use targeted queries first rather than scanning the entire registry.

#### 1. Targeted status check (preferred for exact IDs)

```bash
# Check a single GSE — shows registration, data, config, pipeline status
python -m core.paper.registry status --gse GSE123456

# Check a paper and all its linked datasets
python -m core.paper.registry status --pmid 31493975
```

This is the **primary entry point** for any exact GSE or PMID query.
One command replaces piecing together 3+ separate checks, and would have
prevented the "not registered" misdiagnosis described earlier.

#### 2. Fuzzy paper search (vague description)

```bash
# Search by keyword in title, author, journal, year, PMID, slug
python -m core.paper.registry list-papers --query "retina development"
python -m core.paper.registry list-papers --query "author:Norrie"
python -m core.paper.registry list-papers --query "2024"
```

Use `list-papers` when the user gives a vague description (no exact PMID/GSE).
Once narrowed to a candidate, use `status --pmid` for the full picture.

#### 3. Global registry commands

```bash
python -m core.paper.registry report          # summary counts
python -m core.paper.registry verify          # consistency check
python -m core.paper.registry find-orphans    # orphan datasets
```

#### 4. Decision tree from `status` output

| Status shows | Next step |
|---|---|
| Data downloaded + config exists | `python core/run_pipeline.py --modality <mod> --config <path>` |
| Data downloaded + no config | Generate config via `core/preprocess/` or copy from `templates/config_templates/` |
| Data not downloaded | Download data (GEO) then repeat status check |
| Not registered | `register --gse <GSE>` or `register --pmid <PMID>` first |
| PMID not in registry | `python core/paper/insights.py --pmid <PMID>` then register |

#### 5. Registration

```bash
python -m core.paper.registry register --gse GSE123456              # GSE → PMID auto-link
python -m core.paper.registry register --gse GSE123456 --dry-run    # preview first
python -m core.paper.registry register --pmid 31493975              # paper → GSE linkage
```

#### 6. Programmatic query (for scripts / advanced use)

```python
from core.paper.registry import load_master_registry
reg = load_master_registry()
reg.get_dataset_links("41578023")     # paper → datasets
reg.get_paper_links("GSE123456")      # dataset → papers
reg.get_paper(paper_id="41578023")    # by PMID/slug
reg.find_orphans()

### Methodology Pattern Analysis

Extract and compare the methodological "fingerprint" of any single-cell/omics paper —
domain-agnostic across retina, cancer, immunology, and beyond.

The framework captures 5 dimensions:
- **archetype** — what type of study (atlas, development, perturbation, disease_comp, multiomic, ...)
- **strategy** — experimental design (species, tissue, cell count, modalities, validation)
- **narrative** — how the story is told (argument_structure: bottom_up / top_down / comparative / ...)
- **toolbox** — computational pipeline (framework, algorithms, AI usage)
- **contribution** — what the paper leaves to the field (novelty, reusable assets, field impact)

```bash
# Single paper
python core/paper/insights.py --pmid <PMID> --methodology

# Batch backfill all registered papers
python adhoc/migration_scripts/methodology_batch.py
python adhoc/migration_scripts/methodology_batch.py --dry-run    # preview first
python adhoc/migration_scripts/methodology_batch.py --workers 8  # custom concurrency

# Read methodology_patterns from insights.yaml
python -c "import yaml; d=yaml.safe_load(open('projects/papers/<paper>/insights.yaml')); print(d.get('methodology_patterns', {}).get('archetype', {})"

# Schema reference
# templates/schemas/methodology_patterns.yaml
```

### Critical conventions

- Use `.venv/bin/python` for all Python commands
- Source `.env` first: `set -a && source .env && set +a`
- Steps run as **subprocesses** — never import step scripts directly
- Every step must prepend repo root to `sys.path`
- Config: `.yaml` + `resolve_config()` → Pydantic v2; `.py` rejected
- `data_root()` requires `FUXI_DATA_ROOT` env var
- Import pattern: `from core.utils import ...`, `from core.ai.caller import ...`

## Session report generation

When the user **explicitly asks** to summarize the session or write a report
("写报告", "写总结", "生成报告", "总结一下", "记录一下", "write a report",
"summarize"), generate a markdown report under `notes/` covering:

- What was found / the problem
- How it was solved
- Results and outcomes
- Gaps, caveats, or future work

### Topic → directory mapping

Analyze the session's dominant subject and place the report accordingly:

```
Dominant subject                 → Directory        Naming format
──────────────────────────────────────────────────────────────────
Bug diagnosis / system audit    → audit/            YYYY-MM-DD_<topic>.md
Architecture change / migration → engineering/      YYYY-MM-DD_<topic>.md
Feature implementation          → features/         YYYY-MM-DD_<topic>.md
Technical research / lit review → research/         YYYY-MM-DD_<topic>.md
Knowledge-base update           → kb/               YYYY-MM-DD_<topic>.md
Paper insights / supplements    → supplements/      YYYY-MM-DD_<topic>.md
Reproduction verification       → reproduction/     YYYY-MM-DD_<topic>.md
Reference docs / indices        → reference/        YYYY-MM-DD_<topic>.md
Work log / weekly summary       → logs/             recent_work_summary_YYYY-MM-DD.md
```

Boundary decisions:
- **audit vs engineering** — if both bug diagnosis AND fix implementation are
  present, choose by dominant purpose: evaluation → audit, implementation →
  engineering.
- **features vs engineering** — small-scoped local changes → features;
  system-level architecture changes → engineering.

### Naming rules

- **`YYYY-MM-DD_<topic>.md`** — date is for sorting only; `<topic>` is a short
  slug identifying the report content (e.g. `atac_pipeline_rewrite`). Use the
  report completion date.
- **`recent_work_summary_YYYY-MM-DD.md`** — the date *describes* the covered
  time period. This format is used only for work-log entries under `logs/`.

### Edge cases

- If the session spans multiple topics, pick the dominant one. When truly
  unsure, ask the user.
- If the user asks for a report on a past event (not the current session),
  produce it with the event date, not today's date.
- After writing any report, update `notes/INDEX.md` to add the new entry
  under the corresponding directory section.
- If the notes/ topic classification itself changes, update this mapping.

## Notes are private — never commit

The entire `notes/` directory is **gitignored**. Notes contain internal
thinking, incomplete analysis, and experimental reports — these are for
internal reference only and must **never** appear in the public repository.

- **NEVER** run `git add notes/` or `git add -f notes/` under any circumstance.
- **NEVER** include notes/ files in commits, Pull Requests, or pushes.
- Before every commit, verify with `git status` that no notes/ paths are staged.
- If a commit hook or `git add -A` accidentally stages notes/ files, unstage
  them before committing.
