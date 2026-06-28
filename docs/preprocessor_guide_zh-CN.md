# Fuxi 预处理脚本使用指南

> 适用于：**单细胞组学研究人员** | 无需编程背景即可上手

---

## 目录

1. [这个脚本做什么？](#1-这个脚本做什么)
2. [前置准备](#2-前置准备)
3. [四种使用场景](#3-四种使用场景)
4. [运行结果：你会得到什么？](#4-运行结果你会得到什么)
5. [生成文件后：如何运行完整 Pipeline](#5-生成文件后如何运行完整-pipeline)
6. [高级选项](#6-高级选项)
7. [常见问题（FAQ）](#7-常见问题faq)

---

## 1. 这个脚本做什么？

当你从 GEO / ArrayExpress 或其他来源下载好单细胞数据后，**预处理脚本会自动完成以下工作**：

| 步骤 | 做了什么 | 你原来需要手工做 |
|------|---------|-----------------|
| 🔍 检测文件格式 | 自动识别 10X MTX、HDF5、CSV 矩阵、ATAC fragments…… | 打开文件夹逐个看文件名猜 |
| 📦 解压归档 | 自动解压 `.tar.gz`、`.zip`、`.gz` 等 | 手动 `tar -xzf` 或右键解压 |
| 🧬 检测模态 | 判断是 scRNA-seq、scATAC-seq、还是多模态 (multiome) | 阅读论文 Methods 部分推断 |
| 📋 生成 `dataset.yaml` | 创建数据集元信息清单（样本列表、文件路径、格式） | 手工编辑 YAML |
| ⚙️ 生成 `config_*.py` | 创建可直接运行的 Pipeline 配置文件（模板自动匹配格式） | 从零开始写 ~80 行 Python config |
| 🌐 NCBI 查询（可选） | 获取数据集标题、物种、是否为 SuperSeries | 打开浏览器查 GEO 网页 |

**简单来说：下载完数据 → 运行一个命令 → 配置文件自动生成 → 下一步直接跑 Pipeline。**

---

## 2. 前置准备

### 2.1 安装环境

```bash
# WSL / Linux
cd /path/to/Fuxi
python -m venv .venv
source .venv/bin/activate
pip install -r requirements/rna.txt  # 或 requirements.txt（全部模态）
```

### 2.2 设置数据根目录

```bash
export FUXI_DATA_ROOT=/data/geo_datasets    # Linux 用户
export FUXI_DATA_ROOT=/mnt/c/geo_datasets  # WSL 用户（根据你的实际挂载路径修改）
# 或
export FUXI_DATA_ROOT=/data/geo_datasets    # Linux 用户
```

> 💡 这就是你存放所有 GEO 下载文件的顶层目录。预处理脚本会从 `$FUXI_DATA_ROOT/GSE12345/` 找到对应数据集。

### 2.3 准备 API Key（如果用 AI 注释）

```bash
export LLM_API_KEY=sk-your-api-key-here
```

> 💡 生成 config 文件时 API Key 不会被填入。这是安全考量——你需要手动在 config 中取消注释 AI 配置段。


## 3. 四种使用场景

根据你的实际情况，选择对应的命令：

### 场景 1：标准 GEO 数据集（有 GSE 编号 + 能上网）

```bash
python core/preprocess/preprocessor.py --gse GSE00001 --query-ncbi
```

| 能自动做的事 | 说明 |
|------------|------|
| 文件格式 | ✅ 自动识别 |
| 归档解压 | ✅ tar.gz / zip 等 |
| 物种推断 | ✅ 从 NCBI API 获取（最准确） |
| 数据集标题 | ✅ 从 NCBI API 获取 |
| SuperSeries 检测 | ✅ NCBI + 目录结构 + Series Matrix |
| `dataset.yaml` | ✅ 完整生成 |
| `config_*.py` | ✅ 完整生成 |

**运行后的输出示例：**
```
============================================================
Fuxi Preprocessing: GSE00001
Data root: /data/geo_datasets
============================================================

[Phase 1] Scanning for archives...
  No archives found.
  Total files: 5

[Phase 2] Checking for SuperSeries structure...
  Not a SuperSeries (single accession).

[Phase 3] Detecting file formats...
  Inferred modality: rna

[Phase 4] Generating dataset.yaml...
  Written: projects/rna/GSE00001/dataset.yaml

[Phase 5] Generating config file...
  Written: projects/rna/GSE00001/config_GSE00001.py

============================================================
[Summary] GSE00001
============================================================
  Type:         SingleAccession
  Modality:     rna
  Data format:  10X_mtx
  Species:      homo_sapiens
  Tissue:       retina
  Elapsed:      0.0s

  Generated:
    projects/rna/GSE00001/dataset.yaml
    projects/rna/GSE00001/config_GSE00001.py

  Next steps:
    1. Review and edit the generated files
    2. Run the pipeline:
       python core/run_pipeline.py --modality rna --config projects/rna/GSE00001/config_GSE00001.py
```

---

### 场景 2：有 GSE 编号 + 不能上网（内网 / 离线环境）

```bash
python core/preprocess/preprocessor.py --gse GSE00001
# 不传 --query-ncbi
```

| 相比场景 1 的差异 | 说明 |
|------------------|------|
| 物种推断 | ⚠️ 仅从文件名推断（如 `_human_`、`_mouse_`）。若不匹配，设为 `unknown` |
| 数据集标题 | ❌ 无法获取，`dataset.yaml` 中 title 为空 |
| SuperSeries 子数据集 | ⚠️ 仅靠目录结构检测，可能不完整 |

> ⚠️ **强烈建议**：离线环境下，手动确认 `tissue` 和 `species` 字段。如果预处理脚本设为 `unknown`，请编辑 `dataset.yaml` 和 `config_*.py`。

---

### 场景 3：自己实验室的数据 + 无 GSE 编号 + 能上网

```bash
python core/preprocess/preprocessor.py \
    --input-dir /data/my_retina_organoid_experiment \
    --name retina_organoid_batch1 \
    --query-ncbi
# --query-ncbi 对无 GSE 编号的数据无效（NCBI 查不到），但不会报错
```

| 参数 | 说明 |
|------|------|
| `--input-dir` | 直接指定包含实验数据的目录路径 |
| `--name` | 给数据集起个名字（用于输出目录和文件名）。省略则使用目录名 |

> 💡 `--query-ncbi` 在没有 GSE 编号时不会产生任何作用，可以省略。

---

### 场景 4：纯内网 + 自己实验室的新数据

```bash
python core/preprocess/preprocessor.py \
    --input-dir /data/experiment_20260625 \
    --name my_new_data \
    --no-extract
```

| 相比场景 1 的差异 | 说明 |
|------------------|------|
| 物种推断 | ⚠️ 只能从文件名/元数据文件中推断；若无线索则为 `unknown` |
| 数据集标题 | ❌ 空 |
| GSE 编号 | ❌ 没有（用 `--name` 代替） |
| 归档解压 | 若文件已是标准格式用 `--no-extract` 跳过解压提速 |

> ⚠️ **最重要的差异**：纯内网时，物种可能完全检测不到。请在生成文件后手工填写 `config_*.py` 中的 `CFG.species` 和 `CFG.tissue`。


## 4. 运行结果：你会得到什么？

预处理成功后，`projects/{modality}/{dataset_id}/` 下会生成两个文件：

### 4.1 `dataset.yaml` — 数据集清单

```yaml
id: GSE00001
type: SingleAccession
title: ''                           # ← 需人工填写
species: homo_sapiens
tissue: retina
modalities:
  - name: scRNA-seq
    status: downloaded
    format: 10X_mtx
    file_count: 5
    total_size_gb: 0.0
samples:
  - id: all
    label: ''
    rna:
      - file: GSE00001_Sample1_filtered_feature_bc_matrix.h5
        format: auto
      - file: GSE00001_Sample1_barcodes.tsv.gz
        format: auto
      ...
```

**你需要检查/修改的地方：**

| 字段 | 怎么做 |
|------|--------|
| `title` | 从 GEO 页面粘贴论文标题 |
| `species` | 若为 `unknown`，手工填写（如 `mus_musculus`） |
| `tissue` | 若为 `unknown`，手工填写（如 `retina`、`brain`） |

### 4.2 `config_GSE00001.py` — Pipeline 配置文件

```python
from core.config import CFG

CFG.data_format = '10X_mtx'
CFG.mtx_prefix = 'GSE00001_Sample1_'
CFG.mtx_dir = ''               # 留空则自动解析

CFG.tissue = 'retina'        # ← 人工确认
CFG.species = 'human'

# CFG.sample_map = {        # ← 若有多样本，需人工填写
#     1: 'sample1',
# }

# CFG.marker_dict = {        # ← 需人工填写
#     'CellTypeA': ['GENE1', 'GENE2'],
# }

# CFG.ai.enabled = True      # ← 取消注释以启用 AI 注释
# CFG.ai.api_base = 'https://api.deepseek.com/v1'
```

**你必须做的修改（标有 `# TODO`）：**

| 需要修改 | 重要性 | 说明 |
|---------|--------|------|
| `CFG.marker_dict` | 🔴 必须 | 填写你的组织已知标记基因。若组织有 KB 支持，改为设置 `CFG.tissue_kb` |
| `CFG.sample_map` | 🟡 多样本时 | 映射 10X barcode 后缀 → 样本名 |
| `CFG.stage_map` | 🟡 发育数据时 | 映射样本 → 发育阶段 |
| `CFG.tissue_kb` | 🟢 推荐 | 如为 `retina`/`hypothalamus`，设为对应的 KB 名称即可跳过 `marker_dict` |
| AI 设置 | 🟢 推荐 | 取消注释 AI 段落，填入 API Key |

> 💡 **KB 模式优先**：如果你的组织在 `rna/tissue_ontologies/` 下有对应知识库，只需设置 `CFG.tissue_kb = "retina"` 即可跳过 `marker_dict`。KB 模式比简单打分准确度更高。


## 5. 生成文件后：如何运行完整 Pipeline

### 5.1 scRNA-seq 全流程

```bash
# Step 1: 跑完所有 12 步骤
python core/run_pipeline.py --modality rna --config projects/rna/GSE00001/config_GSE00001.py

# Step 2: 对某个特定细胞类型做子聚类（可选）
python core/run_pipeline.py --modality rna --step 7 --cell-type "Müller Glia" \
    --config projects/rna/GSE00001/config_GSE00001.py
```

### 5.2 scATAC-seq 全流程

```bash
python core/run_pipeline.py --modality atac --config projects/atac/GSE00001/config_GSE00001.py
```

### 5.3 多模态（multiome）数据集

如果数据集同时包含 RNA 和 ATAC（如某个多模态数据集），预处理脚本会自动生成 **两份** config：

```
projects/rna/GSE00001/config_GSE00001.py    ← RNA 配置
projects/atac/GSE00001/config_GSE00001.py   ← ATAC 配置
```

先分别跑 RNA 和 ATAC Pipeline：
```bash
python core/run_pipeline.py --modality rna  --config projects/rna/GSE00001/config_GSE00001.py
python core/run_pipeline.py --modality atac --config projects/atac/GSE00001/config_GSE00001.py
```

然后 ATAC Step 09 会自动发现 RNA 的结果并进行整合。

### 5.4 断点续跑

```bash
python core/run_pipeline.py --modality rna --resume --config projects/rna/GSE00001/config_GSE00001.py
```

### 5.5 只跑某一步

```bash
# 列出所有步骤
python core/run_pipeline.py --modality rna --list

# 只跑注释步骤
python core/run_pipeline.py --modality rna --step 6 --config projects/rna/GSE00001/config_GSE00001.py
```

### 5.6 完整工作流总结

```
1. 从 GEO 下载数据
       ↓
2. 预处理（本文档）
   python core/preprocess/preprocessor.py --gse GSE12345
       ↓
3. 检查生成的 dataset.yaml + config_*.py
   编辑 marker_dict / tissue_kb / AI 设置
       ↓
4. 运行 Pipeline
   python core/run_pipeline.py --modality rna --config ...
       ↓
5. 分析结果
   results/figures/   → UMAP、热图、轨迹图
   results/tables/    → 注释表、差异基因、富集通路
```


## 6. 高级选项

### 完整命令行参数

```
python core/preprocess/preprocessor.py --help
```

| 参数 | 说明 | 适用场景 |
|------|------|---------|
| `--gse GSE12345` | GEO 数据集编号 | 场景 1、2 |
| `--input-dir /path/` | 直接指定输入目录 | 场景 3、4 |
| `--name my_label` | 自定义数据集名称 | 场景 3、4 |
| `--data-root /path/` | 覆盖 `FUXI_DATA_ROOT` | 所有场景 |
| `--query-ncbi` | 查询 NCBI API 获取元数据 | 场景 1 |
| `--dry-run` | 仅检测和报告，不写文件 | 检查效果时 |
| `--force` | 覆盖已有文件 | 重新生成时 |
| `--no-extract` | 跳过归档解压 | 文件已解压时 |
| `--modality rna\|atac` | 强制指定模态 | 多模态数据集分开处理时 |
| `--output-dir /path/` | 指定输出目录 | 不想污染 `projects/` 目录时 |
| `--verbose` / `-v` | 显示详细检测信息 | 排查问题时 |
| `--quiet` / `-q` | 最简输出 | 批量处理时 |

### 查看效果但不写文件

```bash
python core/preprocess/preprocessor.py --gse GSE12345 --dry-run --verbose
```

> 💡 **建议**：第一次使用时先用 `--dry-run` 确认检测结果符合预期，再用 `--force` 正式生成。


## 7. 常见问题（FAQ）

### Q1: 预处理脚本报 "Directory not found"

```
[ERROR] Directory not found: /data/geo_datasets/GSE12345
```

**原因**：
1. `FUXI_DATA_ROOT` 环境变量未设置
2. 数据集尚未下载到 `$FUXI_DATA_ROOT/GSE12345/`

**解决**：
```bash
echo $FUXI_DATA_ROOT                    # 确认已设置
ls $FUXI_DATA_ROOT/GSE12345/            # 确认数据已存在
```

或者使用 `--input-dir` 直接指定：
```bash
python core/preprocess/preprocessor.py --input-dir /path/to/data --name my_data
```

### Q2: 生成的文件里有 `# TODO` 标记，怎么处理？

这是**正常的**。预处理脚本只能自动填充它能确定的信息。标记 `# TODO` 的部分需要你根据实际数据手工填写：

- **`CFG.marker_dict`**：查找文献中你目标组织的已知标记基因
- **`CFG.sample_map`**：从 GEO Metadata 中整理 barcode → 样本的对应关系
- **`CFG.stage_map`**：如果你的实验有时间序列/发育阶段，定义阶段映射

### Q3: 预处理把我的文件识别成错误的格式怎么办？

```bash
# 先查看详细检测结果
python core/preprocess/preprocessor.py --gse GSE12345 --dry-run -v

# 如果检测错误，手动指定模态
python core/preprocess/preprocessor.py --gse GSE12345 --modality atac
```

然后检查生成的 `config_*.py`，手工修改 `CFG.data_format` 和文件路径。

### Q4: 我的数据是 SuperSeries（包含多个子数据集），预处理能处理吗？

能。预处理脚本会自动检测 SuperSeries（通过目录结构 / Series Matrix / NCBI API）。

对于 SuperSeries（如某个 SuperSeries 数据集），生成的 `dataset.yaml` 会列出所有子数据集，但 **config 文件不会为子数据集自动生成**——你需要对每个子数据集分别运行一次预处理。

### Q5: 我的数据格式很特别（非标准命名），预处理能识别吗？

部分能。支持的自定义格式包括：
- 带前缀的 10X 文件（如 `GSE12345_my_sample_filtered_feature_bc_matrix.h5`）
- 自定义的 counts matrix（如 `*_counts.mtx.gz`）
- 每个样本单独 CSV/TSV 文件（如 GSM 文件）

如果仍然无法识别，你可以：
1. 参考 `templates/config_templates/` 中的模板手工写 config
2. 在 `projects/{modality}/{GSE_ID}/` 下创建对应的 `config_*.py`

### Q6: 预处理会覆盖我已有的 config 文件吗？

**不会。** 除非你传了 `--force`。

```bash
# 生成到临时目录查看效果
python core/preprocess/preprocessor.py --gse GSE12345 \
    --output-dir /tmp/preprocess_output

# 确认无误后再正式输出
python core/preprocess/preprocessor.py --gse GSE12345 --force
```

### Q7: 环境变量 `FUXI_DATA_ROOT` 没设置怎么办？

有两种方法：

```bash
# 方法 1: 设置环境变量（推荐）
export FUXI_DATA_ROOT=/data/geo_datasets   # 替换为你的实际数据目录

# 方法 2: 每次运行时指定
python core/preprocess/preprocessor.py --gse GSE12345 --data-root /data/geo_datasets

# 方法 3: 使用 --input-dir（无需 FUXI_DATA_ROOT）
python core/preprocess/preprocessor.py --input-dir /data/geo_datasets/GSE12345
```

### Q8: 批量处理多个数据集怎么写？

```bash
# 在 bash 中循环
for gse in GSE00001 GSE00002 GSE00003; do
    python core/preprocess/preprocessor.py --gse "$gse" --force
done
```

---

