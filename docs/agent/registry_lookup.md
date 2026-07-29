# Dataset & Paper Lookup

How to look up datasets, papers, and registry entries. Use targeted queries first rather than scanning the entire registry.

## 1. Targeted status check (preferred for exact IDs)

```bash
# Check a single GSE — shows registration, data, config, pipeline status
python -m core.paper.registry status --gse GSE123456

# Check a paper and all its linked datasets
python -m core.paper.registry status --pmid 31493975
```

This is the **primary entry point** for any exact GSE or PMID query.
One command replaces piecing together 3+ separate checks.

## 2. Fuzzy paper search (vague description)

```bash
# Search by keyword in title, author, journal, year, PMID, slug
python -m core.paper.registry list-papers --query "retina development"
python -m core.paper.registry list-papers --query "author:Norrie"
python -m core.paper.registry list-papers --query "2024"
```

Use `list-papers` when the user gives a vague description (no exact PMID/GSE).
Once narrowed to a candidate, use `status --pmid` for the full picture.

## 3. Global registry commands

```bash
python -m core.paper.registry report          # summary counts
python -m core.paper.registry verify          # consistency check
python -m core.paper.registry find-orphans    # orphan datasets
```

## 4. Decision tree from `status` output

| Status shows | Next step |
|---|---|
| Data downloaded + config exists | `python core/run_pipeline.py --modality <mod> --config <path>` |
| Data downloaded + no config | Generate config via `core/preprocess/` or copy from `templates/config_templates/` |
| Data not downloaded | Download data (GEO) then repeat status check |
| Not registered | `register --gse <GSE>` or `register --pmid <PMID>` first |
| PMID not in registry | `python core/paper/insights.py --pmid <PMID>` then register |

## 5. Registration

```bash
python -m core.paper.registry register --gse GSE123456              # GSE → PMID auto-link
python -m core.paper.registry register --gse GSE123456 --dry-run    # preview first
python -m core.paper.registry register --pmid 31493975              # paper → GSE linkage
```

## 6. Programmatic query (for scripts / advanced use)

```python
from core.paper.registry import load_master_registry
reg = load_master_registry()
reg.get_dataset_links("41578023")     # paper → datasets
reg.get_paper_links("GSE123456")      # dataset → papers
reg.get_paper(paper_id="41578023")    # by PMID/slug
reg.find_orphans()
```
