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

## Notes 库维护惯例（Mandatory conventions）

写报告只是第一步。以下 5 条为**强制**惯例，每份报告生成时逐条执行，
以保证 notes/ 作为可检索工作库长期有序（索引不腐、过时可见、待办可勾）。

### 规则 1 — 状态头必填

每份新报告在 `# 标题` 之后紧跟状态头（blockquote 格式，**全角冒号** `：`）：

```
> 状态：现行
```

若本报告取代了旧报告，旧报告的状态头改为：

```
> 状态：已被 <后继报告名> 取代
> 后继：<后继报告名>
> 关联：<演化链上的其他报告，逗号分隔>
```

被移入归档区的旧报告标记为：

```
> 状态：已归档→notes/archive/
> 后继：<后继报告名>
```

示例（演化链末端 CURRENT 报告）：

```
> 状态：现行
> 关联：2026-06-30_cluster_selection_pareto_elbow.md, 2026-07-12_multi_metric_clustering_mmacs.md
```

### 规则 2 — BACKLOG 勾选

报告完成某条待办项时，**同一 session 内**更新 `notes/BACKLOG.md`：把对应条目
状态改为 `done`，并在 `备注` 追加实施报告链接。示例（勾选前后）：

```
- [ ] at-001 写 private marker gate 失败测试 — 来源: research/2026-08-03_transition_annotation_todo.md | 状态: open
- [ ] at-001 写 private marker gate 失败测试 — 来源: research/2026-08-03_transition_annotation_todo.md | 状态: done | 备注: archive/2026-08-03_transition_annotation_p0.md
```

勾选后 `grep -c '来源: ' notes/BACKLOG.md` 仍应等于条目总数（125），保证格式
不因勾选而破坏。

### 规则 3 — 日志周期

每两周写一篇工作日志 `notes/logs/recent_work_summary_YYYY-MM-DD.md`，沿用既有
模板（如 `notes/logs/recent_work_summary_2026-07-17.md`）：H1 + `> 日期/涵盖时间范围/
提交数/核心主题` 元数据 blockquote + TOC + 编号主题段 + 提交表 + 页脚。`提交数` 取
`git log --since=<起始日> --until=<截止日> --oneline | wc -l`（只读窗口计数）。

### 规则 4 — 归档策略

被取代的报告**移入 `notes/archive/`（永不删除）**，文件名保持不变，`notes/INDEX.md`
保留指针指向归档位置与后继。只做 `mv`，不合并、不重写正文；归档文件仍需状态头
（见规则 1）。

### 规则 5 — INDEX 更新

每份报告完成后**同一 session 内**更新 `notes/INDEX.md`：在对应目录段新增一行条目，
并刷新该段计数。索引条目必须与磁盘一一对应——无幽灵条目（只在索引、不在磁盘）、
无缺漏（在磁盘、不在索引）。

### 规则 6 — 链导读维护

每个主题演化链（见 `notes/INDEX.md`「演化链」节）在 `notes/reference/chains/` 下有
对应链导读 `chain_<id>.md`。当某条链**新增成员报告**时（新报告与链内既有报告构成
取代/后继关系），**同一 session 内**更新对应导读：

- 时间线表追加该成员行（basename/日期/关键贡献）；
- 关键决策点表补充新决策（来源标注为成员 basename + 章节）；
- 状态头 `> 关联：` 行追加新成员 basename；
- 若链尾变化（新成员成为权威指针），更新「现状 + 权威指针」段。

新建主题链时（首份报告开启新演化链），在同一 session 创建 `chain_<id>.md` 并给
`notes/INDEX.md` 演化链节添加对应条目与 📖 导读链接。导读**不写 commit hash**
（提交可能经 rebase 重写），引用一律用文件名 + 章节锚定。

---

### 边界说明

- **supplements/ 排除**：`notes/supplements/`（论文附表）为使用者单独维护区域——
  不增删、不改动、不纳入索引逐条登记，根 INDEX 仅保留入口。
- **notes/ 永不提交**：整个 `notes/` 目录被 gitignore，属私有内容；**严禁** `git add notes/`
  或以任何形式把 notes/ 文件带进 commit / PR / push（详见 AGENTS.md「Notes are private — never commit」）。
