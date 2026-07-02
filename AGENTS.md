# Fuxi (伏羲) — Agent Quick Reference

> Full knowledge base: [CLAUDE.md](CLAUDE.md)

## One-liners

```bash
python core/run_pipeline.py --modality rna --list
python core/run_pipeline.py --modality atac --list
python core/run_pipeline.py --modality rna --config projects/rna/<GSE_ID>/config_<GSE_ID>.py
```

## Key paths

| Module | Location |
|--------|----------|
| Shared core | `core/` (config, utils, ai_caller, ai_prompts, run_pipeline, preprocess) |
| RNA steps | `rna/steps/` (12 scripts) |
| ATAC steps | `atac/steps/` (10 scripts) |
| Project configs | `projects/{modality}/{GSE_ID}/config_*.py` |

## Critical conventions

- Steps run as **subprocesses** via `run_pipeline.py` — never imported directly
- Every step script must add repo root to `sys.path`: `sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))`
- Config loaded dynamically: `CFG = resolve_config(args.config)`
- `data_root()` requires `FUXI_DATA_ROOT` env var (no hardcoded defaults)
- Import pattern: `from core.utils import ...`, `from core.ai_caller import ...`

## Commit message discipline

- Commit message 中的每条 `Cx`/`Mx`/`Nx`/`mx`/`Sx` 声明必须在 `git show <commit> --stat` 或 diff 中**直接验证**
- 若某 commit message 中包含的声明确实未实现，应在该 commit 的 `git notes` 中标记 `UNFIXED`
