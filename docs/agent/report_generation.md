# Session Report Generation

When the user **explicitly asks** to summarize the session or write a report
("写报告", "写总结", "生成报告", "总结一下", "记录一下", "write a report",
"summarize"), generate a markdown report under `notes/` covering:

- What was found / the problem
- How it was solved
- Results and outcomes
- Gaps, caveats, or future work

## Topic → directory mapping

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

## Boundary decisions

- **audit vs engineering** — if both bug diagnosis AND fix implementation are
  present, choose by dominant purpose: evaluation → audit, implementation →
  engineering.
- **features vs engineering** — small-scoped local changes → features;
  system-level architecture changes → engineering.

## Naming rules

- **`YYYY-MM-DD_<topic>.md`** — date is for sorting only; `<topic>` is a short
  slug identifying the report content (e.g. `atac_pipeline_rewrite`). Use the
  report completion date.
- **`recent_work_summary_YYYY-MM-DD.md`** — the date *describes* the covered
  time period. This format is used only for work-log entries under `logs/`.

## Edge cases

- If the session spans multiple topics, pick the dominant one. When truly
  unsure, ask the user.
- If the user asks for a report on a past event (not the current session),
  produce it with the event date, not today's date.
- After writing any report, update `notes/INDEX.md` to add the new entry
  under the corresponding directory section.
- If the notes/ topic classification itself changes, update this mapping.
