# Fuxi 分析管线使用指南

> 适用于：**单细胞组学研究人员** | 无需编程背景即可上手

---

## 目录

1. [管线能做什么？](#1-管线能做什么)
2. [前置准备](#2-前置准备)
3. [快速开始：运行第一条管线](#3-快速开始运行第一条管线)
4. [scRNA-seq 管线详解](#4-scrna-seq-管线详解)
5. [scATAC-seq 管线详解](#5-scatac-seq-管线详解)
6. [空间转录组管线详解](#6-空间转录组管线详解)
7. [结果文件说明](#7-结果文件说明)
8. [常用运行技巧](#8-常用运行技巧)
9. [配置文件详解](#9-配置文件详解)
10. [常见问题（FAQ）](#10-常见问题faq)

---

## 1. 管线能做什么？

当你从 GEO 等公共数据库下载好单细胞数据、运行完预处理脚本得到配置文件后，分析管线将自动完成从原始数据到生物学结论的**全流程计算**：

| 阶段 | 做了什么 | 生物学意义 |
|------|---------|-----------|
| 🔬 数据加载 | 自动识别并读取 6 种常见单细胞数据格式 | 统一为内部格式，屏蔽格式差异 |
| 🧹 质量控制 | 去除双细胞、死细胞、低质量细胞 | 保证下游分析基于可靠数据 |
| 🔗 批次整合 | 归一化 + 高变基因 + PCA + 批次校正 | 消除技术差异，保留生物学信号 |
| 🗺️ 聚类与可视化 | 多参数网格搜索 + UMAP 降维 | 发现细胞亚群，呈现数据结构 |
| 🏷️ 细胞注释 | KB 打分 / AI 大模型 / 标记基因三种模式自动注释 | 将聚类编号转化为有生物学意义的细胞类型 |
| 🔍 差异分析 | 标记基因 + 阶段比较 + 时序趋势三层 DE | 鉴定各类细胞的特征基因和发育动态 |
| 🌳 轨迹推断 | PAGA + 扩散伪时间 | 重建细胞分化/发育路径 |
| 🧬 通路富集 | GO/KEGG 过表达分析 + GSEA | 揭示细胞类型的生物学功能 |
| 🧭 GRN 调控网络 | 伪细胞聚合 + decoupler TF 活性推断 | 鉴定驱动各类细胞身份的转录因子 |

**简单来说：配置文件就绪 → 一条命令 → 从原始数据到论文级图表全部自动产出。**

---

## 2. 前置准备

### 2.1 安装环境

```bash
# Linux / WSL
cd /path/to/Fuxi
python -m venv .venv
source .venv/bin/activate
pip install -r requirements/rna.txt  # 或 requirements.txt（全部模态）
```

### 2.2 设置环境变量

```bash
# 必需：数据根目录（存放所有下载的原始数据）
export FUXI_DATA_ROOT=/data/geo_datasets

# WSL 用户注意：必须关闭 HDF5 文件锁
export HDF5_USE_FILE_LOCKING=FALSE

# 可选：AI 注释所需的 API Key
export LLM_API_KEY=sk-your-api-key-here
```

### 2.3 确认配置文件就绪

运行管线前，你需要确保数据集目录下已有两个文件：

```
projects/{模态}/{数据集ID}/
├── dataset.yaml          # 数据集元信息清单
└── config_{数据集ID}.py   # 管线配置文件
```

这两个文件通常由**预处理脚本**自动生成。如果你还没有配置文件，请先参考《Fuxi 预处理脚本使用指南》。

---

## 3. 快速开始：运行第一条管线

### 3.1 查看可用步骤

```bash
# 查看 scRNA-seq 管线的所有步骤
python core/run_pipeline.py --modality rna --list

# 查看 scATAC-seq 管线的所有步骤
python core/run_pipeline.py --modality atac --list

# 查看空间转录组管线的所有步骤
python core/run_pipeline.py --modality spatial --list
```

你会看到类似这样的输出：

```
Fuxi — RNA-seq pipeline step list
============================================================
  [00] Load raw data → 00_raw.h5ad
  ...
  [11] GRN regulatory network analysis (decoupler) → 11_grn.h5ad
```

### 3.2 一键运行全流程

```bash
# scRNA-seq 全流程（12 步，从数据加载到 GRN 分析）
python core/run_pipeline.py --modality rna --config projects/rna/{数据集ID}/config_{数据集ID}.py

# scATAC-seq 全流程（10 步）
python core/run_pipeline.py --modality atac --config projects/atac/{数据集ID}/config_{数据集ID}.py

# 空间转录组全流程（10 步）
python core/run_pipeline.py --modality spatial --config projects/spatial/{数据集ID}/config_{数据集ID}.py
```

运行过程中，终端会实时显示每一步的进度和耗时：

```
============================================================
[run] [RNA] Step [00]: Load raw data → 00_raw.h5ad
============================================================
[run] Step [00] completed (took 45.2s).

============================================================
[run] [RNA] Step [02]: Scrublet doublet detection (per sample) → 01_doublet.h5ad
============================================================
[run] Step [02] completed (took 120.8s).
...
============================================================
[run] Fuxi RNA-seq pipeline execution finished.
============================================================
[run] Step timing summary:
  [00]    45.2s  Load raw data → 00_raw.h5ad
  ...
  [Total] 1845.3s  10 steps total
```

### 3.3 检查点与断点续跑

管线采用**检查点机制**：每步完成后会保存中间结果文件。如果中途因故中断，可以使用 `--resume` 从断点继续，已完成的步骤会自动跳过：

```bash
python core/run_pipeline.py --modality rna --resume --config projects/rna/{数据集ID}/config_{数据集ID}.py
```

### 3.4 运行单个步骤

如果你想单独运行或重跑某一步：

```bash
# 只运行第 06 步（细胞注释）
python core/run_pipeline.py --modality rna --step 5 --config projects/rna/{数据集ID}/config_{数据集ID}.py

# 运行步骤 02 到 05
python core/run_pipeline.py --modality rna --steps 2-5 --config projects/rna/{数据集ID}/config_{数据集ID}.py

# 运行步骤 00, 02, 04（跳着跑）
python core/run_pipeline.py --modality rna --steps 0,2,4 --config projects/rna/{数据集ID}/config_{数据集ID}.py
```

---

## 4. scRNA-seq 管线详解

scRNA-seq 管线包含 12 个步骤（编号 00-11），数据依次流转：

```
原始数据 → 00_load → 01_doublet → 02_qc
         → 03_integrate → 04_cluster → 05_annotate
         → 06_subcluster → 07_markers → 08_trajectory
         → 09_enrichment → 10_exploratory → 11_grn
```

### Step 00：数据加载

**输入**：原始数据文件（自动识别格式） | **输出**：`00_raw.h5ad`

管线自动识别并加载以下格式之一：

| 格式 | 典型文件 | 适用场景 |
|------|---------|---------|
| `10X_h5` | `*filtered_feature_bc_matrix.h5` | Cell Ranger 标准输出 |
| `10X_mtx` | `matrix.mtx.gz` + `barcodes.tsv.gz` + `features.tsv.gz` | Cell Ranger 原始输出 |
| `csv_matrix` | 基因×细胞计数矩阵（CSV/TSV/MTX） | 自定义实验、Smart-seq2 等 |
| `h5ad` | `*.h5ad` | 已预处理的数据 |

> **R 格式（`.rds` / `.qs`）**：管线不原生支持。可使用 [r2h5ad](https://github.com/nanshanteahouse/r2h5ad) 转为 h5ad 后加载。

加载过程中自动完成：
- **样本/阶段映射**：根据 barcode 后缀（如 `-1`、`-2`）自动标注每个细胞的样本来源和发育阶段
- **多文件合并**：如果某个数据集有多个 H5 文件，自动合并为一个 AnnData 对象
- **格式兼容**：自动处理 legacy 2 列 genes.tsv 与标准 3 列 features.tsv 的转换

### Step 01：双细胞检测（Scrublet）

**输入**：`00_raw.h5ad` | **输出**：`01_doublet.h5ad`

按样本独立运行 Scrublet 算法，检测因液滴中包裹了两个细胞而产生的"双细胞"假象。

- 大样本（>15,000 细胞）串行处理以避免内存溢出
- 小样本通过并行加速
- **非 count 格式数据自动跳过**：当 `expression_type` 为 `TPM`/`FPKM`/`CPM`/`log1p_counts` 时，自动禁用 Scrublet（因为负二项分布假设对归一化后数据不成立）
- 如果某样本的 Scrublet 运行失败，该样本的所有细胞标记为非双细胞（优雅降级，不阻塞管线）

输出：`doublet_scores`（双细胞概率分数）和 `predicted_doublet`（是/否）两列。

### Step 02：质量控制（QC）

**输入**：`01_doublet.h5ad` | **输出**：`02_qc.h5ad`<br>
**诊断图**：`{figure_dir}/02_qc/` — `nFeature_distribution.png`、`nCount_vs_nFeature.png`、`pct_mito_distribution.png`

两种模式，由 `use_adaptive_thresholds` 控制：

| 模式 | 配置 | 行为 |
|------|------|------|
| 硬阈值（默认） | `use_adaptive_thresholds=False` | 使用 config 中的固定阈值 |
| MAD 自适应 | `use_adaptive_thresholds=True` | 每个指标用 Median ± N × MAD 计算自适应上/下界，硬阈值作为安全地板/天花板 |

过滤维度（全部表示为 `(lo, hi)` 的阈值字典）：

1. **去除双细胞**：剔除 `predicted_doublet=True` 的细胞
2. **基因数过滤**：去除检测基因数过少（空液滴）或过多（双细胞漏网）的细胞（默认 500-7500）
3. **线粒体过滤**：去除线粒体基因占比 > 20% 的细胞（死细胞/受损细胞）
4. **nCount 上限**（仅 raw_counts）：对 `total_counts` 做上界过滤，TPM/FPKM/CPM 下自动跳过
5. **复杂度过滤**（仅 raw_counts）：`log10(基因数)/log10(总UMI)` 下界过滤，TPM/FPKM/CPM 下自动跳过

**始终生成 3 张诊断图**，标注实际使用的阈值线（硬阈值或 MAD），无需人工看图决策即可提供可溯源的审计追踪。

### Step 03：归一化与批次整合

**输入**：`02_qc.h5ad` | **输出**：`03_integrated.h5ad`

这是管线中最关键的整合步骤，将不同样本/批次的数据对齐到同一个分析空间：

1. **高变基因筛选（HVG）**：识别信息量最大的基因（默认约 4000 个）。自动尝试多种方法（`seurat_v3` → `cell_ranger` → `seurat` → 手动方差筛选），确保在任何数据上都能成功
2. **回归协变量**（可选）：去除 `total_counts` 和 `pct_counts_mt` 等技术因素对表达量的影响
3. **归一化**：每个细胞的总计数归一化到 10,000，然后 log1p 变换
4. **PCA 降维**：主成分分析（默认 100 维），生成肘部图
5. **Harmony 批次校正**：消除不同样本/批次间的技术差异（可选，通过 `CFG.use_harmony` 控制）

> 💡 完整基因表达矩阵保存在 `.raw` 中，下游标记基因计算和差异分析均使用 `.raw`，确保不会因 HVG 筛选丢失信息。

### Step 04：聚类与 UMAP

**输入**：`03_integrated.h5ad` | **输出**：`04_clustered.h5ad`

执行**多参数网格搜索**，自动找到最优聚类方案：

- 遍历多种 `n_neighbors` 参数（默认 [15, 20, 30]）和 Leiden 分辨率（默认 [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]）
- 每种组合计算一次 UMAP 降维和 Leiden 聚类
- 以 **Pareto 拐点算法**（Pareto 前沿 + 归一化拐点检测）为标准，在聚类数和 Silhouette 分数之间自动选出最优参数
- 生成所有参数组合的 UMAP 对比图（PDF）和网格搜索汇总表（CSV）

> 💡 这一步的计算量较大，因为它尝试了 3×6=18 种参数组合。但这是值得的——你不用手动反复尝试不同参数。

### Step 05：细胞类型注释

**输入**：`04_clustered.h5ad` | **输出**：`05_annotated.h5ad`

这是管线的核心步骤，将"Cluster 0, Cluster 1, ..."转化为"Rod Photoreceptor, Bipolar Cell, Müller Glia, ..."等有生物学意义的标签。管线支持三种注释模式，按优先级自动选择：

#### 模式一：KB 知识库模式（最高准确度）

如果你研究的是已有知识库支持的组织（如视网膜），只需在配置中设置：
```python
CFG.tissue_kb = "retina"
```

管线将自动执行一个精密的多层决策流程：
1. **标记基因计算**：对每个聚类，用 Wilcoxon 秩和检验找出高表达基因
2. **知识库打分**：将各聚类的标记基因与知识库中各细胞类型的已知标记基因进行超几何检验 + 余弦相似度双重打分
3. **专家规则匹配**：应用优先级排序的确定性匹配规则（如"同时高表达 RHO 和 PDE6A → Rod Photoreceptor"）
4. **证据融合**：5 层决策引擎综合 marker 分数、专家规则、层级结构等信息，给出带置信度的共识注释
5. **AI 兜底**（可选）：对低置信度聚类，自动调用大语言模型（LLM）基于标记基因重新推理
6. **质量控制**：自动标记线粒体/核糖体主导的低质量聚类为 "Unknown"，生成注释质量报告

#### 模式二：AI 大模型模式

如果启用 AI 注释：
```python
CFG.ai.enabled = True
CFG.ai.ai_annotation = True
```

管线会将每个聚类的 top 标记基因发送给大语言模型，模型基于生物学知识推断细胞类型，并返回结构化的注释结果（细胞类型 + 亚型 + 状态 + 置信度 + 推理过程）。

支持多种 LLM 后端：OpenAI API、DeepSeek、vLLM、Ollama 等（通过 `CFG.ai.api_base` 配置）。

#### 模式三：Score_genes 简单打分（兜底）

如果既没有知识库也没有 AI，管线自动回退到经典的标记基因打分模式——你只需要在配置中提供 `CFG.marker_dict`（手工整理的各细胞类型的标记基因列表）。

> 💡 **注释结果包含**：`cell_type`（主类型）、`cell_subtype`（亚型）、`cell_state`（状态）、`annot_confidence`（置信度）、`annot_reasoning`（推理过程）。

### Step 06：子聚类分析（可选）

**输入**：`05_annotated.h5ad` | **输出**：`05_sub_{细胞类型}.h5ad`

对某个特定细胞类型（如"Müller Glia"）进行精细亚型分析：

```bash
python core/run_pipeline.py --modality rna --step 7 \
    --cell-type "Müller Glia" \
    --config projects/rna/{数据集ID}/config_{数据集ID}.py
```

该步骤将在指定细胞类型的子集上重新运行 PCA → 邻居图 → UMAP → Leiden 聚类，并可选择性地使用 AI 对亚群进行重注释。结果自动写回主 `05_annotated.h5ad` 文件的 `cell_subtype` 列。

### Step 07：差异表达分析（三层 DE）

**输入**：`05_annotated.h5ad` | **输出**：CSV 表格 + 热图/点图

执行三个层次的差异表达分析：

**第一层：标记基因** — 每类细胞 vs 所有其他细胞（Wilcoxon 秩和检验）
- 输出：`marker_genes_per_group.csv`
- 用于鉴定各类细胞的身份标志基因

**第二层：阶段配对比较** — 同一类细胞在相邻发育阶段之间的比较（t 检验）
- 输出：`pairwise_stage_de.csv`
- 追踪细胞在发育过程中的转录变化（需要数据有 `stage` 注释）
- 自动并行处理多组比较

**第三层：时间趋势基因** — 基因表达随发育时间的 Spearman 相关
- 输出：`temporal_trend_genes.csv`
- 每类细胞取前 20 个上调基因和前 20 个下调基因
- 要求至少 3 个发育阶段

> 💡 同时生成：标记基因热图（每类 top 5 基因）、已知标记基因点图。

### Step 08：轨迹分析（PAGA + DPT）

**输入**：`04_clustered.h5ad`（通常使用） | **输出**：`05_final.h5ad`

重建细胞的发育/分化轨迹：

1. **PAGA 图**：构建细胞类型间的拓扑连接图，揭示分化关系
2. **根细胞识别**：自动确定轨迹起点（优先级：配置指定根类型 → 根标记基因最高表达 → 最早发育阶段 → 数据集第一个细胞）
3. **扩散伪时间（DPT）**：沿 PAGA 拓扑计算每个细胞的伪时间位置
4. **分支分析**：鉴定分化分支点的谱系特异性基因
5. **发育基因可视化**：绘制已知发育基因（SOX2、PAX6、NEUROD1 等）沿伪时间的表达趋势

### Step 09：GO/KEGG 通路富集

**输入**：Step 07 的标记基因表 | **输出**：CSV 表格 + 气泡图 + AI 解读

两种互补的富集方法：

| 方法 | 原理 | 输出 |
|------|------|------|
| **ORA**（过表达分析） | 取 top 标记基因，超几何检验看哪些通路被显著富集 | `enrichment_ora.csv` |
| **GSEA**（基因集富集分析） | 所有基因按分数排序，检验通路基因是否在列表顶/底部聚集 | `enrichment_gsea.csv` |

支持 200+ 基因集库，常用包括：
- `GO_Biological_Process` — 基因本体生物学过程
- `KEGG_2021_Human` — KEGG 代谢/信号通路
- `Reactome_2022` — Reactome 通路数据库
- `MSigDB_Hallmark_2020` — 标志性基因集

> 💡 如果启用了 AI，富集结果会自动生成一段生物学解读报告（`ai_interpretation.txt`）。

### Step 10：探索性分析

**输入**：`05_annotated.h5ad` | **输出**：CSV 表格 + 多种 PDF 图表

生成全面的探索性汇总图表，帮助你快速理解数据全貌：

- **细胞组成分析**：各类细胞在不同样本/阶段的占比堆叠柱状图
- **QC 指标可视化**：基因数、UMI 数、线粒体比例在 UMAP 上的分布
- **标记基因表达**：已知标记基因在 UMAP 上的表达热图
- **聚类统计**：各聚类/细胞类型的细胞数及占比

### Step 11：GRN 调控网络分析（decoupler）

**输入**：`05_annotated.h5ad` | **输出**：`11_grn.h5ad` + CSV 表格 + 热图

基于已注释细胞类型的伪细胞聚合进行转录因子（TF）活性推断：

1. **伪细胞聚合**：按 `cell_type` 取平均表达量，平滑单细胞 dropout 噪声
2. **Regulon 网络**：通过 decoupler 获取 CollecTRI 数据库（~1,185 个 TF，含带符号的靶基因调控关系）
3. **TF 活性推断**：运行 ULM（单变量线性模型）——检验每个细胞类型中 TF 的靶基因是否在高表达基因中显著富集
4. **输出**：
   - `11_grn.h5ad` — 伪细胞 AnnData（obs=细胞类型，var=基因），含 `obsm['X_tf_activity']`
   - `tables/11_grn/tf_activity_per_cell_type.csv` — 完整 TF 活性矩阵（细胞类型 × TF）
   - `tables/11_grn/tf_activity_pvals.csv` — 对应的 P 值矩阵
   - `tables/11_grn/tf_target_edges.csv` — top-variance TF 的 TF→靶基因调控边表
   - `tables/11_grn/tf_target_counts.csv` — 各 TF 的靶基因数量汇总
   - `figures/11_grn/tf_activity_heatmap.png` — 按方差选出的 top N 个 TF 的聚类热图

**配置参数：**

```python
CFG.run_grn = True               # 启用/关闭此步骤
CFG.grn_method = "decoupler"     # 方法（目前仅 decoupler；pySCENIC 待定）
CFG.grn_species = "human"        # 'human' | 'mouse'
CFG.grn_n_top_regulons = 50      # 热图中显示的方差最高 TF 数量
CFG.grn_min_regulon_size = 5     # 每个 regulon 的最少靶基因数
CFG.grn_confidence_levels = ["A","B","C"]  # DoRothEA 置信度等级（若使用 DoRothEA）
```

> 💡 **无需下载额外的数据库文件。**decoupler 在首次使用时联网获取 regulon 网络并缓存到本地。CollecTRI（默认）比 DoRothEA 覆盖更多 TF，是推荐的网络。

### Step 12：细胞互作分析（LIANA+）

**输入**：`05_annotated.h5ad` | **输出**：CSV 表格 + 热图/点图

基于配体-受体（LR）数据库的细胞间通讯推断，通过置换检验评估每对细胞类型之间的 LR 互作显著性：

1. **共识打分**：`liana.mt.rank_aggregate` 整合 9 种方法（CellPhoneDB、CellChat、NATMI 等）的评分，取稳健秩聚合（RRA）得到共识互作排序
2. **输出**：
   - `tables/12_cell_interaction/cci_interactions.csv` — 全量 LR 互作评分（含 source、target、ligand、receptor、magnitude_rank 等列）
   - `tables/12_cell_interaction/cci_top_interactions.csv` — 按 magnitude_rank 排序的 top N 显著对
   - `figures/12_cell_interaction/cci_heatmap.png` — source→target 细胞类型互作数量热图
   - `figures/12_cell_interaction/cci_dotplot.png` — top LR 对的配体-受体互作点图

**配置参数：**

```python
CFG.run_cci = True                 # 启用/关闭此步骤
CFG.cci_lr_database = "consensus"   # LR 数据库：'consensus' | 'cellphonedb' | 'cellchat' | ...
CFG.cci_permutations = 1000         # 置换检验迭代次数
CFG.cci_n_top_interactions = 50     # top N 互作对
# CFG.cci_spatial_method = "liana_spatial"  # (spatial) reserves 'commot'
# CFG.cci_multi_condition = False           # future: 多条件差异 CCI
```

> 💡 **首次运行会自动下载 LR 数据库（~200MB）至 ~/.cache/liana 并缓存。**如果 var_names 是 Ensembl ID 格式，管线会自动通过 mygene 转换为 gene symbol。

---

### 空间转录组 CCI（Spatial Step 10）

空间管线在 Step 10 中额外提供了**空间约束的配体-受体共表达分析**（`liana.mt.bivariate`）：

- **局部指标**：空间加权的余弦相似度——检测在相邻 spot 中配体和受体基因是否协同表达
- **全局指标**：Moran's R 及置换 p-value——评估配体-受体对的空间自相关显著性
- **输出**：
  - `tables/10_cell_interaction/cci_spatial_interactions.csv`
  - `tables/10_cell_interaction/cci_spatial_top.csv`
  - `figures/10_cell_interaction/cci_spatial_heatmap.png` — ligand×receptor 热图（Moran's R）
  - `figures/10_cell_interaction/cci_spatial_dotplot.png` — top LR 对条形图

**空间专属配置参数：**

```python
CFG.cci_spatial_method = "liana_spatial"  # 方法（reserves 'commot'）
CFG.cci_spatial_distance = 0.0            # 空间距离阈值（0 = 使用已有 spatial_connectivities）
```

---

## 5. scATAC-seq 管线详解

scATAC-seq 管线包含 10 个步骤（编号 00-09）：

```
原始数据 → 00_load → 01_qc → 02_process → 03_cluster
         → 04_annotate → 05_marker_peaks → 06_motif
         → 07_trajectory → 08_enrichment → 09_integrate
```

### Step 00：数据加载

**输入**：原始 ATAC 数据 | **输出**：`00_raw.h5ad`

支持三种格式：

| 格式 | 典型文件 | 说明 |
|------|---------|------|
| `10x_fragments` | `*fragments.tsv.gz` | ATAC 片段文件，通过 SnapATAC2 流式导入 |
| `10x_peak_h5` | `*filtered_peak_bc_matrix.h5` | 10X 峰-细胞矩阵 HDF5 |
| `h5ad` | `*.h5ad` | 预处理过的 AnnData |

### Step 01：QC + Peak Calling

**输入**：`00_raw.h5ad` | **输出**：`01_filtered.h5ad`

1. **片段数过滤**：去除片段数异常（过少=空液滴，过多=双细胞）的细胞
2. **MACS3 Peak Calling**：在 pseudobulk 上调用 MACS3 识别开放染色质区域（仅保留标准染色体 chr1-22, X, Y）
3. **峰-细胞矩阵构建**：生成二值化的可及性矩阵（开放=1, 关闭=0）
4. **双细胞检测**：在峰矩阵上运行 Scrublet

### Step 02：特征选择与降维

**输入**：`01_filtered.h5ad` | **输出**：`02_processed.h5ad`

1. **去除双细胞**
2. **IDF 特征选择**：选取信息量最大的峰（默认 top 50k），降低噪声
3. **谱嵌入**：基于矩阵的 Lanczos 谱分解（30 维）
4. **KNN 邻居图**：在谱空间中构建细胞邻居关系

### Step 03：聚类与 UMAP

**输入**：`02_processed.h5ad` | **输出**：`03_clustered.h5ad`

与 RNA 管线相同的多参数网格搜索策略（遍历 `n_neighbors` × `resolution`），在谱嵌入空间中以 Pareto 拐点算法自动选择最优聚类。

### Step 04：AI 染色质状态注释

**输入**：`03_clustered.h5ad` | **输出**：`04_annotated.h5ad`

1. 计算每个聚类的差异可及性峰（marker regions）
2. 将 top 峰及其邻近基因发送给大语言模型（LLM）
3. LLM 基于染色质开放区域的基因关联推断细胞类型/染色质状态
4. AI 响应自动磁盘缓存（SHA256 去重），重复运行不额外调用 API

### Step 05-08：下游分析

| 步骤 | 内容 | 输出 |
|------|------|------|
| Step 05 | 差异可及性峰 per 细胞类型 | `marker_peaks.csv` |
| Step 06 | TF 结合基序富集（CIS-BP 数据库） | `motif_enrichment_{细胞类型}.csv` |
| Step 07 | ATAC 伪时间轨迹分析 | `07_trajectory.h5ad` |
| Step 08 | 峰关联基因的 GO/KEGG 富集 | `enrichment_*.csv` |

### Step 09：RNA+ATAC 整合

**输入**：ATAC `04_annotated.h5ad` + RNA h5ad（自动发现） | **输出**：`09_integrated.h5ad`

如果你有配对的多组学数据（同一批细胞的 RNA-seq + ATAC-seq）：
1. 自动发现同数据集下的 RNA 结果文件
2. 通过 barcode 交集找到共有的细胞
3. 构建 MuData 多模态对象（`rna` + `atac` 两个 modality）
4. 运行联合 PCA

> 💡 对于**多组学（multiome）数据集**，预处理脚本会自动生成 RNA 和 ATAC 两份配置。先分别跑完 RNA 和 ATAC 管线，再跑 ATAC Step 09 即可自动整合。

---

## 6. 空间转录组管线详解

空间转录组管线包含 10 个步骤（编号 00-09），针对 10X Visium 数据设计（可扩展至其他平台）：

```
原始数据 → 00_load → 01_qc → 02_image → 03_normalize
         → 04_cluster → 05_annotate → 06_spatial_de
         → 07_trajectory → 08_enrichment → 09_exploratory
```

详细的空间管线运行报告（含遇到的实际问题与修复方法），见 `notes/suggestions/spatial_<GSE_ID>.md`。

### 支持平台

| 平台 | 配置值 | 说明 |
|------|-------|------|
| 10X Visium | `"visium"` | SpaceRanger 输出目录或含空间坐标的 h5ad |
| Slide-seq | `"slideseq"` | 磁珠空间条形码 |
| MERFISH | `"merfish"` | 成像型，基因面板 |
| seqFISH | `"seqfish"` | 成像型，基因面板 |

### 输入格式

支持两种数据格式：

| 格式 | 配置 | 输入 |
|------|------|------|
| 10X Visium 目录 | `CFG.data_format = "visium"` | 包含 `filtered_feature_bc_matrix.h5` + `spatial/` 的目录 |
| 预构建 h5ad | `CFG.data_format = "h5ad"` | 含 `obsm['spatial']` 坐标和 `uns['spatial']` 图像的 `.h5ad` |

### Step 00：数据加载

**输入**：原始 Visium 数据或 h5ad | **输出**：`00_raw.h5ad`

- `visium` 格式：使用 `sq.read.visium()`，自动检测 `library_id`
- `h5ad` 格式：使用 `sc.read()` 读取并验证空间坐标
- 自动转换为稀疏 CSR 格式，确保观察名称唯一
- 为 Visium 数据添加默认 `in_tissue` 标记

### Step 01：质量控制

**输入**：`00_raw.h5ad` | **输出**：`01_qc.h5ad`

应用适配空间数据的 QC 指标：
- 基因数过滤（默认 200-7,500）
- 线粒体百分比过滤（默认 <25%）
- 基因-UMI 复杂度过滤
- 过滤后 spot 数过低时发出警告

### Step 02：图像处理

**输入**：`01_qc.h5ad` | **输出**：`02_image.h5ad`

从组织 H&E/IF 图像中提取图像特征：
- 从 `uns['spatial']` 自动检测 library_id
- 将图像裁剪至组织区域（可配置）
- 通过 `sq.im.process()` 提取基本图像特征（纹理、直方图）
- 若无图像则优雅降级（跳过处理，继续管线）

### Step 03：归一化

**输入**：`02_image.h5ad` | **输出**：`03_processed.h5ad`

- 每 spot 文库大小归一化至 10,000
- Log1p 变换
- 在 `.raw` 中保留原始计数
- PCA 降维

### Step 04：聚类与 UMAP

**输入**：`03_processed.h5ad` | **输出**：`04_clustered.h5ad`

- 多分辨率 Leiden 聚类（分辨率网格搜索）
- UMAP 嵌入（2D）
- 空间邻居图构建（`sq.gr.spatial_neighbors`）
- 网格搜索汇总保存至 `param_grid_summary.csv`

### Step 05：细胞类型注释

**输入**：`04_clustered.h5ad` | **输出**：`05_annotated.h5ad`

核心注释步骤。三种注释模式，按优先级选择：

#### 模式 1：KB 知识库模式（最高精度）

若 `CFG.tissue_kb` 已设置（如 `"retina"`），管线复用 RNA 管线的完整注释引擎：
- 为每个空间 cluster 计算标记基因
- 将 cluster 与组织知识库进行打分比对
- 应用专家确定性规则
- 跨打分层级的证据融合
- 对低置信度 cluster 进行 AI 兜底

#### 模式 2：AI 大模型模式

若 `CFG.ai.enabled` 和 `CFG.ai.ai_annotation` 已设置：
- 将每个 cluster 的标记基因发送给大模型
- 返回结构化注释（cell_type、subtype、state、confidence、reasoning）

#### 模式 3：Score_genes 简单打分（回退）

使用 `CFG.marker_dict` 进行逐 cluster 标记基因打分。此模式可通过以下方式增强：
- 配置文件中**用户自定义的标记基因**
- 通过 Phase 1 标记列表迁移获取的 **scRNA 来源标记基因**（见下文）

#### Phase 1：scRNA 标记列表迁移（新功能）

当已在匹配样本（相同组织、相同时间点）上运行过 scRNA-seq 时，可将每个细胞类型的标记基因自动迁移至空间注释：

```python
# 在空间配置中：
CFG.rna_ref = "<RNA_dataset_id>"    # 例如 "GSE235585"
CFG.rna_marker_top_n = 10           # 每个细胞类型取 top-N 标记基因
CFG.rna_marker_pval_threshold = 0.05
CFG.rna_marker_logfc_min = 0.0      # 0 = 仅正向 LFC
```

工作原理：
1. 从匹配的 RNA 项目中自动发现 scRNA 的 `marker_genes_per_group_cell_type.csv`
2. 提取每个细胞类型的 top-N 显著性标记基因
3. 合并进 `CFG.marker_dict`（用户配置的条目优先）
4. 增强 `score_genes_mode()` 回退，不影响 KB 或 AI 模式

> 💡 这是**标记列表迁移**，而非单细胞粒度的标签迁移。它利用 scRNA 管线已完成计算的差异表达结果来辅助空间 cluster 注释。零新增依赖。

### Step 06：空间 DE + 空间可变基因

**输入**：`05_annotated.h5ad` | **输出**：`marker_genes_per_group.csv`、`svg_rankings.csv`、`06_svg.h5ad`

- 逐 cluster 差异表达（Wilcoxon 秩和检验）
- 空间自相关 Moran's I 检测空间可变基因（SVG）
- Top SVG 空间散点图
- SVG 标记写入 `adata.var['spatially_variable']`

### Step 07：轨迹分析

**输入**：`05_annotated.h5ad` | **输出**：`07_trajectory.h5ad`

- 在空间 cluster 上构建 PAGA 图
- 扩散伪时间（DPT）
- 通过标记基因或细胞类型识别起始根细胞

### Step 08：GO/KEGG 富集

**输入**：Step 05 产出的 `marker_genes_per_group.csv` | **输出**：富集 CSV + 气泡图

复用 RNA 管线的富集引擎：
- ORA（过表达分析）通过 Enrichr API
- 预排序 GSEA（本地计算）
- 支持 200+ 基因集库（GO、KEGG、Reactome、MSigDB Hallmark 等）

### Step 09：探索性分析

**输入**：`05_annotated.h5ad` + `06_svg.h5ad` | **输出**：图表 + CSV 表格

- 细胞类型在组织坐标上的空间散点图
- 基因表达空间图谱（top 标记基因 + SVG）
- Spot 组成统计（cluster / 细胞类型大小）
- 空间邻居图摘要（边数、平均度）
- UMAP 汇总图

### 空间专用配置字段

```python
CFG.spatial_platform = "visium"           # visium | slideseq | merfish | seqfish
CFG.library_id = ""                      # Visium 库 ID（为空时自动检测）
CFG.crop_image = True                    # 将图像裁剪至组织区域
CFG.spatial_neighbors_n = 6              # 空间邻居数
CFG.spatial_neighbors_radius = 0.0       # 半径模式（0 = 使用 n_neighbors）
CFG.run_spatial_autocorr = True          # 运行 Moran's I SVG 检测
CFG.svg_n_top = 2000                     # 下游分析保留的最大 SVG 数

# Phase 1: scRNA 标记列表迁移
CFG.rna_ref = ""                         # scRNA 项目路径或 dataset_id
CFG.rna_marker_top_n = 10                # 每个细胞类型取 top-N 标记基因
CFG.rna_marker_pval_threshold = 0.05     # pvals_adj 阈值
CFG.rna_marker_logfc_min = 0.0           # 最小 logfoldchanges
```

---

## 7. 结果文件说明

管线运行完毕后，所有结果统一存放在数据集的 `results/` 目录下，按类型分为三个子目录：

```
results/
├── h5ad/                          # 中间检查点文件
│   ├── 00_raw.h5ad                # 原始数据
│   ├── 01_doublet.h5ad            # 双细胞检测后
│   ├── 02_qc.h5ad                 # QC 过滤后
│   ├── 03_integrated.h5ad         # 批次整合后
│   ├── 04_clustered.h5ad          # 聚类 + UMAP 后
│   ├── 05_annotated.h5ad          # 细胞注释后 ★
│   ├── 05_final.h5ad              # 轨迹分析后
│   └── 11_grn.h5ad               # 伪细胞 + TF 活性 (GRN) ★
│
├── figures/                       # 可视化图表
│   ├── 02_qc/                     # QC 诊断图
│   │   ├── nFeature_distribution.png
│   │   ├── nCount_vs_nFeature.png
│   │   └── pct_mito_distribution.png
│   ├── 03_integrate/              # 批次整合
│   │   ├── pca_elbow.png
│   │   └── harmony_comparison.png
│   ├── 04_cluster/                # 聚类 + UMAP
│   │   ├── umap_param_grid_summary.png
│   │   ├── umap_min_dist_comparison.png   # min_dist 扫描对比
│   │   ├── umap_grid_n*_r*.png
│   │   └── umap_leiden_n*_all_resolutions.pdf
│   ├── 05_annotation/             # 细胞注释
│   │   └── _05_celltype*.pdf
│   ├── 06_subcluster/             # 亚群分析
│   ├── 07_markers/                # 标记基因
│   │   ├── _07_marker_heatmap.pdf
│   │   └── _07_dotplot.pdf
│   ├── 08_trajectory/             # 轨迹分析
│   │   ├── _08_pseudotime.pdf
│   │   ├── _08_paga_umap.pdf
│   │   └── _08_dev_genes_heatmap.pdf
│   ├── 09_enrichment/             # 富集分析
│   │   ├── ora_*_bubble.pdf
│   │   └── prerank_*_bubble.pdf
│   ├── 10_exploratory/            # 探索性分析
│   │   ├── composition_by_stage_*.png
│   │   └── _06_marker_dotplot.pdf
│   └── 11_grn/                    # GRN 调控网络
│       └── tf_activity_heatmap.png
│
└── tables/                        # 数据表格
    ├── marker_genes_per_group.csv # 标记基因
    ├── param_grid_summary.csv     # 聚类参数网格
    ├── cell_type_annotations.csv  # 细胞类型注释
    ├── cell_metadata.csv          # 细胞元数据
    ├── marker_genes_ai.csv        # AI 模式标记基因
    ├── marker_genes_unified.csv   # 统一 KB 模式标记基因
    ├── 05_annotation_quality.json # 注释质量评估
    ├── pairwise_stage_de.csv      # 阶段配对差异表达
    ├── temporal_trend_genes.csv   # 时间趋势基因
    ├── branch_deg.csv             # 分支 DEG
    ├── cell_type_sizes.csv        # 细胞类型统计
    ├── 09_enrichment/             # 富集分析
    │   ├── ora_*_summary.csv
    │   ├── prerank_*_summary.csv
    │   ├── ai_interpretation.txt
    │   └── ai_interpretation_summary.txt
    ├── 10_exploratory/            # 探索性分析
    │   └── composition_by_stage_*.csv
    └── 11_grn/                    # GRN 调控网络
        ├── tf_activity_per_cell_type.csv
        ├── tf_activity_pvals.csv
        ├── tf_target_edges.csv
        └── tf_target_counts.csv
```

> 💡 加 ★ 的 `05_annotated.h5ad` 是最重要的输出文件——它包含每个细胞的最终注释标签，是绝大多数下游分析（差异表达、轨迹、富集）的起点。

---

## 8. 常用运行技巧

### 7.1 查看管线进度

```bash
# 列出所有步骤及其对应的检查点文件
python core/run_pipeline.py --modality rna --list --config projects/rna/{数据集ID}/config_{数据集ID}.py
```

### 7.2 从断点恢复

```bash
# 自动检测第一个未完成的步骤，从那里继续
python core/run_pipeline.py --modality rna --resume --config projects/rna/{数据集ID}/config_{数据集ID}.py
```

管线会自动扫描检查点文件，跳过已完成的步骤。无论中断原因是网络问题、内存不足还是手动终止，都可以用同一条命令恢复。

### 7.3 只重跑特定步骤

如果你对某一步的结果不满意，调整配置参数后：

```bash
# 只重跑注释步骤（Step 05）
python core/run_pipeline.py --modality rna --step 5 --config projects/rna/{数据集ID}/config_{数据集ID}.py

# 重跑后面的所有步骤
python core/run_pipeline.py --modality rna --steps 6-11 --config projects/rna/{数据集ID}/config_{数据集ID}.py
```

### 7.4 跳过慢步骤

如果你只关注部分分析结果：

```bash
# 只跑加载到注释（前 7 步）
python core/run_pipeline.py --modality rna --steps 0-6 --config projects/rna/{数据集ID}/config_{数据集ID}.py
```

### 7.5 清理中间文件

管线运行过程中会产生多个中间检查点文件（每个 h5ad 可能数百 MB 到数 GB）。如果你磁盘空间有限，可以在每步完成后自动删除上游文件：

```bash
python core/run_pipeline.py --modality rna --cleanup --config projects/rna/{数据集ID}/config_{数据集ID}.py
```

### 7.6 子聚类分析

对某个已注释的细胞类型进行精细亚型分析：

```bash
python core/run_pipeline.py --modality rna --step 7 \
    --cell-type "Müller Glia" \
    --config projects/rna/{数据集ID}/config_{数据集ID}.py
```

你可以对多个细胞类型分别运行，结果会自动合并回主注释文件。

### 7.7 多组学数据处理

```bash
# 第一步：分别跑 RNA 和 ATAC
python core/run_pipeline.py --modality rna  --config projects/rna/{数据集ID}/config_{数据集ID}.py
python core/run_pipeline.py --modality atac --config projects/atac/{数据集ID}/config_{数据集ID}.py

# 第二步：ATAC Step 09 自动整合
# 在 ATAC 全流程的最后一步自动完成，或单独运行：
python core/run_pipeline.py --modality atac --step 9 --config projects/atac/{数据集ID}/config_{数据集ID}.py
```

### 7.8 跨组学 scRNA → spatial 标记迁移

如果你有匹配样本的 scRNA-seq 数据，可将每个细胞类型的标记基因迁移到空间注释中：

```bash
# 第一步：先跑 RNA 管线（产出 marker_genes_per_group_cell_type.csv）
python core/run_pipeline.py --modality rna --config projects/rna/{RNA数据集ID}/config_{RNA数据集ID}.py

# 第二步：在空间配置中设置 rna_ref 指向 RNA 数据
# 在空间配置文件中：
#   CFG.rna_ref = "{RNA数据集ID}"
# 第三步：跑空间管线——注释步骤自动使用 scRNA 标记基因
python core/run_pipeline.py --modality spatial --config projects/spatial/{空间数据集ID}/config_{空间数据集ID}.py
```

`--list` 命令会显示 scRNA 标记文件是否自动发现成功。

---

## 9. 配置文件详解

配置文件（`config_{数据集ID}.py`）是一个 Python 脚本，通过修改全局 `CFG` 对象来控制管线的所有行为。以下是最常需要调整的配置项：

### 8.1 数据输入配置

```python
# 数据格式（必须与你的文件格式匹配）
CFG.data_format = '10X_mtx'     # 选项: 10X_h5, 10X_mtx, csv_matrix, h5ad, 10x_fragments, 10x_peak_h5

# 10X MTX 格式相关
CFG.mtx_dir = ''               # 留空则自动解析为 $FUXI_DATA_ROOT/{数据集ID} (推荐)
CFG.mtx_prefix = 'sample_'     # MTX 文件前缀

# CSV 格式相关
CFG.matrix_file = 'counts.csv.gz'
CFG.barcodes_file = 'barcodes.tsv.gz'
CFG.features_file = 'features.tsv.gz'
```

### 8.2 样本与阶段映射

```python
# barcode 后缀 → 样本名称
CFG.sample_map = {
    1: 'Control',
    2: 'Treatment',
}

# 样本名称 → 发育阶段
CFG.stage_map = {
    'Control':   'E14.5',
    'Treatment': 'P0',
}
CFG.stage_order = ['E14.5', 'P0']  # 阶段的时间顺序（用于时序分析）
```

### 8.3 注释配置（核心）

```python
# 方式一：使用知识库（如果你的组织有 KB 支持）
CFG.tissue_kb = "retina"

# 方式二：手工指定标记基因（KB 不支持时）
CFG.marker_dict = {
    'Rod Photoreceptor': ['RHO', 'PDE6A', 'NRL', 'GNAT1'],
    'Bipolar Cell':      ['VSX2', 'GRIK1', 'TRPM1', 'CABP5'],
    'Müller Glia':       ['GLUL', 'RLBP1', 'CLU', 'VIM'],
    # ... 按你的组织添加
}

# 方式三：启用 AI 注释
CFG.ai.enabled = True
CFG.ai.ai_annotation = True
CFG.ai.model = "deepseek-chat"              # 模型名称
CFG.ai.api_base = "https://api.deepseek.com/v1"
```

### 8.4 质量控制

```python
CFG.expression_type = "raw_counts"     # raw_counts | TPM | FPKM | CPM | log1p_counts
                                       # TPM/FPKM/CPM: 自动跳过 total_counts 和复杂度过滤
                                       # 非 raw_counts 数据自动禁用 Scrublet
CFG.min_genes = 500                    # 细胞最少检测基因数
CFG.max_genes = 7500                   # 细胞最多检测基因数（去除双细胞漏网）
CFG.max_pct_mito = 20.0                # 最大线粒体百分比
CFG.min_genes_per_umi = 0.70           # 复杂度阈值 — 仅在 raw_counts 下生效
CFG.min_cells_per_gene = 3             # 基因在最少几个细胞中表达
CFG.use_adaptive_thresholds = False    # True → 基于 MAD 的阈值（自动适应数据分布）
CFG.mad_n_mads = 3.0                   # MAD 倍数
CFG.qc_ncount_max_mad = 5.0            # nCount 上限的 MAD 倍数（更宽，因为 nCount 有重尾）
```

### 8.5 批次校正

```python
CFG.use_harmony = True              # 启用 Harmony 批次校正
CFG.harmony_batch_key = "sample"    # 按哪个列做批次校正
CFG.use_regress_out = False         # 是否回归 total_counts 和 MT%
```

### 8.6 聚类参数

```python
CFG.n_neighbors_grid = [15, 20, 30]        # UMAP 邻居数候选值
CFG.resolution_grid = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]  # Leiden 分辨率候选值
CFG.best_resolution = None                  # 手动模式下设为具体值
CFG.best_n_neighbors = 0                     # 手动模式下设为具体值（0 = 自动）
CFG.cluster_selection_method = "pareto_elbow"  # "pareto_elbow" (默认) | "silhouette" | None

# UMAP 可视化参数扫描——在最佳聚类参数选定后进行，复用 KNN 图（很快）
CFG.umap_selection_method = "convex_hull"  # "convex_hull" (自动，默认) | None (手动)
CFG.param_grid_min_dist = [0.1, 0.3, 0.5]  # min_dist 扫描值
CFG.param_grid_spread = [1.0]              # spread 扫描值
CFG.umap_min_dist = 0.3                    # 手动模式下的固定值
CFG.umap_spread = 1.0
```

### 8.7 物种与基因组

```python
CFG.species = 'human'    # 物种（影响线粒体基因模式、富集分析数据库等）
CFG.tissue = 'retina'    # 组织名称
CFG.genome = 'hg38'      # 参考基因组（ATAC 管线必需）
```

### 8.8 空间转录组

```python
# 平台与输入
CFG.spatial_platform = "visium"           # visium | slideseq | merfish | seqfish
CFG.library_id = ""                      # Visium 库 ID（为空时自动检测）
CFG.data_format = "visium"               # visium | h5ad

# 图像处理
CFG.crop_image = True                    # 将图像裁剪至组织区域

# 空间邻居图
CFG.spatial_neighbors_n = 6              # 空间邻居数
CFG.spatial_neighbors_radius = 0.0       # 半径模式（0 = 使用 n_neighbors）

# 空间可变基因
CFG.run_spatial_autocorr = True          # 运行 Moran's I
CFG.svg_n_top = 2000                     # 下游分析保留的最大 SVG 数

# Phase 1: scRNA 标记列表迁移
CFG.rna_ref = ""                         # scRNA 项目路径或 dataset_id
CFG.rna_marker_top_n = 10                # 每个细胞类型取 top-N 标记基因
CFG.rna_marker_pval_threshold = 0.05     # pvals_adj 阈值
CFG.rna_marker_logfc_min = 0.0           # 最小 logfoldchanges
```

---

## 9. 常见问题（FAQ）

### Q1：运行报 "HDF5 file locking" 错误

```
OSError: Unable to open file (file locking disabled on this file system)
```

**原因**：WSL 环境下访问 Windows 文件系统（`/mnt/` 路径）时，HDF5 的文件锁机制不兼容。

**解决**：运行管线前先设置环境变量：

```bash
export HDF5_USE_FILE_LOCKING=FALSE
```

> 💡 建议将此设置写入 `~/.bashrc` 或 `~/.zshrc`，一劳永逸。

### Q2：预处理生成的配置文件里有 `# TODO` 标记，怎么处理？

预处理脚本只能自动填充它能确定的信息。标记 `# TODO` 的部分需要你根据实验设计手工填写：

- **`CFG.marker_dict`**：查阅文献，整理你目标组织各细胞类型的已知标记基因。如果你的组织有 KB 支持（如 retina），直接设置 `CFG.tissue_kb` 即可
- **`CFG.sample_map`**：从 GEO 页面的 SRA Run Selector 或论文 Methods 部分整理 barcode → 样本对应关系
- **`CFG.stage_map`**：如果你研究发育过程，建立样本 → 发育阶段映射

### Q3：管线的某一步失败了，怎么恢复？

```bash
# 直接使用 --resume，管线会自动从失败的那一步继续
python core/run_pipeline.py --modality rna --resume --config projects/rna/{数据集ID}/config_{数据集ID}.py
```

`--resume` 会扫描检查点文件，自动找到第一个未完成的步骤。已完成的步骤不会重复运行，所以恢复速度很快。

### Q4：我的细胞注释结果不好，怎么办？

首先确认你选择了合适的注释模式：

1. **如果你的组织有 KB 支持**（如 retina）→ 使用 KB 模式：`CFG.tissue_kb = "retina"`
2. **如果你有 AI API Key** → 开启 AI 兜底：`CFG.ai.ai_annotation = True`（AI 会自动接手 KB 模式下低置信度的聚类）
3. **如果是非标准组织** → 手工整理 `CFG.marker_dict`（参考 CellMarker 2.0、PanglaoDB 等数据库）

> 💡 KB + AI 的组合效果最好：KB 提供专家级规则匹配，AI 负责处理未知和前体细胞类型。

### Q5：我的数据是 TPM 格式，QC 把所有细胞都过滤掉了？

TPM 数据中每个细胞的 `total_counts` 恰为 1,000,000，这会导致复杂度指标 `log10(基因数) / log10(1,000,000)` 约等于固定常数，完全失去区分力。

**解决**：在配置中设置 `expression_type = "TPM"`，管线会自动：

1. 跳过 `total_counts` 过滤（TPM/FPKM/CPM 下 nCount 无解释意义）
2. 跳过 `log_genes_per_umi` 复杂度过滤（同上）
3. 自动禁用 Scrublet（归一化后数据违反负二项假设）

```python
CFG.expression_type = "TPM"  # 自动调整 QC 策略
```

无需再手工调整 `min_genes_per_umi`。

### Q6：ATAC 管线的 Peak 并集内存爆炸？

多样本 ATAC 数据中，每个样本独立做 peak calling 会导致 peak 坐标几乎不重叠，并集可能达到百万级甚至更多。

**建议**：如果可能，将多个样本的 fragment 文件合并，用 SnapATAC2 的 MACS3 统一进行 peak calling，获得一致性 peak set。

### Q7：AI 注释的 API 调用会花多少钱？

管线设计了多层**缓存和降级机制**来控制成本：

- **磁盘缓存**：相同输入（SHA256 匹配）不会重复调用 API，重跑免费
- **优雅降级**：API 失败时自动回退到标记基因打分，不会阻塞管线
- **成本可控**：单次 scRNA-seq 注释通常调用 1 次 API（所有聚类一次发送），成本极低

### Q8：我的数据和参考基因组不是人类，怎么处理？

```python
CFG.species = 'mouse'    # 支持: human, mouse, macaque, zebrafish 等
CFG.genome = 'mm10'      # 相应的参考基因组

# 如果有同源基因映射文件，配置正交映射
CFG.ortholog_map = 'path/to/ortholog_map.csv'
```

> 💡 跨物种分析时，KB 注释模式的准确度会随进化距离递减。建议采取 KB + AI 组合模式，让 AI 弥补 KB 在非人物种上的盲区。

### Q9：gseapy 安装失败？

GSEApy 0.11.0+ 需要 Rust 编译器。如果 `pip install gseapy` 失败：

```bash
# 安装 Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# 然后重新安装
pip install gseapy
```

> 💡 如果不想安装 Rust，可以安装较旧版本的 gseapy：`pip install gseapy==0.10.8`

### Q10：如何知道我的数据适合用什么参数？

管线的网格搜索机制已经为你解决了大部分参数选择问题（如 Leiden 分辨率、UMAP 邻居数）。以下参数可能需要根据你的数据特性手工调整：

| 参数 | 调整依据 |
|------|---------|
| `expression_type` | 设为 `"TPM"`/`"FPKM"`/`"CPM"` 可自动调整 QC 策略 |
| `min_genes` / `max_genes` | 查看 QC 步骤生成的分布图，取 1%-99% 分位数附近 |
| `max_pct_mito` | 通常 20% 是安全默认值；代谢活跃组织（如心肌）可适当放宽 |
| `use_adaptive_thresholds` | 设为 `True` 可启用自适应 QC（MAD 法）；`mad_n_mads=3.0` 是较好的起点 |
| `n_pcs` | 查看 `pca_elbow.png`，取肘部位置 |
| `harmony_batch_key` | 通常是 "sample"；如果你的实验包含不同测序平台，也可是 "platform" |

---

## 附录：快速命令参考

### scRNA-seq

```bash
# 全流程
python core/run_pipeline.py --modality rna --config projects/rna/{数据集ID}/config_{数据集ID}.py

# 断点续跑
python core/run_pipeline.py --modality rna --resume --config projects/rna/{数据集ID}/config_{数据集ID}.py

# 单步
python core/run_pipeline.py --modality rna --step 5 --config ...

# 范围
python core/run_pipeline.py --modality rna --steps 4-8 --config ...

# 子聚类
python core/run_pipeline.py --modality rna --step 7 --cell-type "细胞类型名" --config ...

# 列出步骤
python core/run_pipeline.py --modality rna --list --config ...
```

### scATAC-seq

```bash
# 全流程
python core/run_pipeline.py --modality atac --config projects/atac/{数据集ID}/config_{数据集ID}.py

# 断点续跑
python core/run_pipeline.py --modality atac --resume --config projects/atac/{数据集ID}/config_{数据集ID}.py

# 单步
python core/run_pipeline.py --modality atac --step 4 --config ...

# 整合 RNA+ATAC
python core/run_pipeline.py --modality atac --step 9 --config ...
```

### 空间转录组

```bash
# 全流程
python core/run_pipeline.py --modality spatial --config projects/spatial/{数据集ID}/config_{数据集ID}.py

# 断点续跑
python core/run_pipeline.py --modality spatial --resume --config projects/spatial/{数据集ID}/config_{数据集ID}.py

# 单步
python core/run_pipeline.py --modality spatial --step 5 --config ...

# 带 scRNA 标记迁移（Phase 1）
python core/run_pipeline.py --modality spatial --config projects/spatial/{数据集ID}/config_{数据集ID}.py
# 配置文件需设置: CFG.rna_ref = "{RNA数据集ID}"
```
