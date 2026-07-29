# Methodology Pattern Analysis

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
