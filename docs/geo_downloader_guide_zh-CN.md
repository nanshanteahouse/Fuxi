# Fuxi GEO 下载器 — 使用指南

> 适用于：**单细胞组学研究人员** | 无需编程背景即可上手

---

## 目录

1. [这个工具做什么？](#1-这个工具做什么)
2. [前置准备](#2-前置准备)
   - [2.1 安装下载工具（wget 或 curl）](#21-安装下载工具wget-或-curl)
   - [2.2 设置数据根目录（FUXI_DATA_ROOT）](#22-设置数据根目录fuxi_data_root)
   - [2.3 设置 NCBI API Key（可选）](#23-设置-ncbi-api-key可选)
3. [三种使用方式](#3-三种使用方式)
   - [方式 1：独立 CLI](#方式-1独立-cli)
   - [方式 2：预处理脚本带 --download](#方式-2预处理脚本带---download)
   - [方式 3：Registry register --pmid 带 --download](#方式-3registry-register---pmid-带---download)
   - [方式 3.1：选择性数据集注册](#方式-31选择性数据集注册)
   - [方式 3.2：注销（deregister）](#方式-32注销deregister)
4. [独立命令行详解](#4-独立命令行详解)
   - [4.1 完整参数](#41-完整参数)
   - [4.2 用法示例](#42-用法示例)
5. [元数据解析](#5-元数据解析)
   - [5.1 SOFT 文件提供了哪些信息](#51-soft-文件提供了哪些信息)
   - [5.2 物种归一化](#52-物种归一化)
   - [5.3 元数据缓存](#53-元数据缓存)
6. [补充文件列表](#6-补充文件列表)
   - [6.1 文件发现机制](#61-文件发现机制)
   - [6.2 RAW.tar 标记](#62-rawtar-标记)
   - [6.3 多文件系列示例](#63-多文件系列示例)
   - [6.4 没有补充文件时怎么办？](#64-没有补充文件时怎么办)
7. [断点续传](#7-断点续传)
   - [7.1 续传原理](#71-续传原理)
   - [7.2 安全中断（Ctrl+C）](#72-安全中断ctrlc)
   - [7.3 超时设置](#73-超时设置)
   - [7.4 大文件下载建议](#74-大文件下载建议)
8. [Registry 联动](#8-registry-联动)
   - [8.1 自动联动路径](#81-自动联动路径)
   - [8.2 空字段优先策略](#82-空字段优先策略)
   - [8.3 CLI 触发方式](#83-cli-触发方式)
   - [8.4 查看结果](#84-查看结果)
9. [特殊场景](#9-特殊场景)
   - [9.1 SuperSeries 检测](#91-superseries-检测)
   - [9.2 5 位 GSE 编号兼容性](#92-5-位-gse-编号兼容性)
   - [9.3 非单细胞数据（bulk RNA-seq、STARR-seq 等）](#93-非单细胞数据bulk-rna-seqstarr-seq-等)
   - [9.4 Visium 空间转录组（GSM 级文件限制）](#94-visium-空间转录组gsm-级文件限制)
   - [9.5 多模态（Multiome）数据](#95-多模态multiome数据)
   - [9.6 多数据集论文过滤](#96-多数据集论文过滤)
10. [常见问题（FAQ）](#10-常见问题faq)
    - [Q1: 下载中断了怎么办？](#q1-下载中断了怎么办)
    - [Q2: 提示 "Neither wget nor curl found" 怎么处理？](#q2-提示-neither-wget-nor-curl-found-怎么处理)
    - [Q3: 能否只下载部分文件？](#q3-能否只下载部分文件)
    - [Q4: 我的数据集是 SuperSeries，应该下载哪个 GSE？](#q4-我的数据集是-superseries应该下载哪个-gse)
    - [Q5: 下载完了下一步是什么？](#q5-下载完了下一步是什么)
    - [Q6: 下载时提示 "No supplementary files found" 怎么办？](#q6-下载时提示-no-supplementary-files-found-怎么办)
    - [Q7: FUXI_DATA_ROOT 没设置怎么办？](#q7-fuxi_data_root-没设置怎么办)
    - [Q8: 能否批量下载多个数据集？](#q8-能否批量下载多个数据集)
    - [Q9: 工具标记了 SuperSeries，但它有自己的文件，还需要下载子系列吗？](#q9-工具标记了-superseries但它有自己的文件还需要下载子系列吗)

---

## 1. 这个工具做什么？

GEO 下载器负责从 NCBI GEO（Gene Expression Omnibus）自动下载单细胞数据集到本地。它替你完成了以下工作：

| 步骤 | 做了什么 | 你原来需要手工做 |
|------|---------|-----------------|
| 自动解析 GSE 编号 | 从 5 位到 6 位编号全适配 | 查 NCBI FTP 目录结构 |
| 下载 SOFT 元数据 | 获取标题、物种、PMID、样本数…… | 打开浏览器查 GEO 页面 |
| 列出补充文件 | 显示所有 `suppl/` 下的文件及大小 | `wget --spider` 试错 |
| 批量下载 | 自动调用 wget/curl 下载所有文件 | 手动逐个下载 |
| 写入 Registry | 更新 Master Registry 的字段和关联 | 手工编辑 registry YAML |
| 与预处理流程联动 | 下载后自动衔接 `preprocessor.py` | 先手动下载再跑预处理 |

**简单来说：给一个 GSE 编号 → 运行一个命令 → 数据自动下载到 `$FUXI_DATA_ROOT/GSEXXXXXX/` → 下一步直接跑预处理。**

---

## 2. 前置准备

### 2.1 安装下载工具（wget 或 curl）

下载器本身是纯 Python 脚本，但它调用系统的 `wget`（推荐）或 `curl` 完成实际文件传输。

**两者任选其一即可：**

```bash
# Debian / Ubuntu / WSL
sudo apt install wget        # 推荐 — 断点续传体验更好
# 或
sudo apt install curl

# macOS (Homebrew)
brew install wget
# 或
brew install curl
```

> 运行时会自动检测：优先用 `wget`，找不到则用 `curl`。如果两个都没有，会报错提示安装。

### 2.2 设置数据根目录（FUXI_DATA_ROOT）

下载的数据会存放到 `$FUXI_DATA_ROOT/GSEXXXXXX/` 下。

```bash
export FUXI_DATA_ROOT=/data/geo_datasets    # Linux 用户
export FUXI_DATA_ROOT=/mnt/c/geo_datasets  # WSL 用户（根据你的实际挂载路径修改）
```

> 这就是你存放所有 GEO 下载文件的顶层目录。每个 GSE 会自动创建子目录。

你也可以每次运行通过 `--data-root` 覆盖。

### 2.3 设置 NCBI API Key（可选）

不设置也能下载，但 NCBI 默认限制 3 次请求/秒。设置 API Key 后提速到 10 次请求/秒，元数据解析更快。

```bash
export NCBI_API_KEY=your-ncbi-api-key-here
```

[申请 NCBI API Key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-ncbi-e-utilities/) 免费。

---

## 3. 三种使用方式

### 方式 1：独立 CLI

最直接的使用方式，给 GSE 编号就下载：

```bash
python core/geo_downloader.py --gse GSE107618
```

### 方式 2：预处理脚本带 `--download`

在预处理前自动检查并下载数据：

```bash
python core/preprocess/preprocessor.py --gse GSE107618 --modality rna --download
```

预处理脚本的 Phase 0a 会自动调用下载器，数据就绪后进入格式检测和 config 生成。如果数据已存在则跳过下载。

### 方式 3：Registry `register --pmid` 带 `--download`

通过 Registry 注册论文时自动下载关联的 GSE 数据集。`register --pmid` 会先调用 `paper_insights` 解析论文，然后进入交互式数据集选择：

```bash
python -m core.registry register --pmid 31269016 --download
```

流程：论文解读 → 交互式数据集选择 → 注册到 Registry → 自动下载所选数据集。

关键参数：
- `--datasets GSE1,GSE2` — 非交互式，只注册并下载指定数据集
- `--all` — 注册所有数据集（旧版行为，跳过选择界面）
- `--download` — 注册后自动下载数据

### 方式 3.1：选择性数据集注册

当一篇论文关联多个 GSE 数据集时（例如一篇单细胞论文同时提交了 scRNA-seq 和 Hi-C 数据），你可以只选择需要的部分进行注册。

**交互式模式（默认）**：

```bash
python -m core.registry register --pmid 00000000
```

运行后会解析论文，列出所有关联的 GSE 数据集及 SOFT 元数据摘要：

```
Found 2 datasets:
  [1] GSE123456  | single-cell RNA-seq of human retina
  [2] GSE35156   | Hi-C of retinal cell lines

Select datasets (comma-separated numbers, or 'all'): 1
```

输入 `1` 后只注册 GSE123456，GSE35156 被跳过。输入 `all` 或直接回车（如果只有一个数据集）则注册全部。

**非交互式模式（脚本用）**：

```bash
# 只注册指定数据集
python -m core.registry register --pmid 00000000 --datasets GSE123456

# 注册全部（跳过选择）
python -m core.registry register --pmid 00000000 --all
```

`--datasets` 接受逗号分隔的 GSE 编号列表，适合 shell 脚本中批量调用。`--all` 保持旧版行为，注册所有数据集。

**配合 `--download`**：

```bash
# 交互选择后下载所选数据集
python -m core.registry register --pmid 00000000 --download

# 非交互：只下载 GSE123456
python -m core.registry register --pmid 00000000 --datasets GSE123456 --download
```

### 方式 3.2：注销（deregister）

注册错误或不再需要某个数据集/论文时，可以使用 `deregister` 命令从 Registry 中移除。

```bash
# 移除单个数据集
python -m core.registry deregister --gse GSE35156

# 无确认（脚本用）
python -m core.registry deregister --gse GSE35156 --force

# 移除论文 + 级联删除未被其他论文引用的数据集
python -m core.registry deregister --pmid 00000000 --cascade

# 预览（不执行）
python -m core.registry deregister --pmid 00000000 --cascade --dry-run
```

`--cascade` 只删除孤立数据集（不被其他任何论文引用的数据集）。如果某个数据集同时被多篇论文引用，级联删除只会移除数据集与目标论文的关联，但保留数据集条目本身。
---

## 4. 独立命令行详解

### 4.1 完整参数

```bash
python core/geo_downloader.py --help
```

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--gse GSE107618` | GEO 数据集编号（必填） | — |
| `--dry-run` | 只报告，不下载 | `false` |
| `--skip-soft` | 跳过 SOFT 元数据获取 | `false` |
| `--data-root /path/` | 覆盖 `FUXI_DATA_ROOT` | 环境变量 |
| `--force` | 重新下载已存在的文件 | `false` |
| `--quiet` / `-q` | 最小输出 | `false` |

### 4.2 用法示例

**1. 标准下载 — 获取 GSE107618 所有数据：**

```bash
python core/geo_downloader.py --gse GSE107618
```

运行过程会依次显示 SOFT 元数据、文件列表、下载进度和最终汇总。

**2. Dry-run 预览 — 先看再下：**

```bash
python core/geo_downloader.py --gse GSE107618 --dry-run
```

输出示例：

```
============================================================
  GEO Download: GSE107618
============================================================

  [METADATA] GSE107618_family.soft.gz
    Title:    Dissecting the transcriptome landscape of human neural retina...
    Organism: Homo sapiens
    Platform: Illumina NovaSeq 6000
    PMID:     31269016
    Samples:  64

  [SUPPL] Listing files...
    2 file(s) found (271.5 MB total)
    ★ RAW (tar)                        264.8 MB
      GSE107618_Supplementary_Data.xlsx    6.7 MB

  [DRY-RUN] Would download to: /data/geo_datasets/GSE107618/
    → GSE107618_RAW.tar  (264.8 MB)
    → GSE107618_Supplementary_Data.xlsx  (6.7 MB)
```

> 先用 `--dry-run` 确认你要的文件都在，再正式下载。

**3. 跳过元数据下载 — 仅下载数据文件：**

```bash
python core/geo_downloader.py --gse GSE107618 --skip-soft
```

如果之前已经下载过 SOFT 文件（缓存在 `.geo_meta.json`），跳过可节省一次 NCBI 请求。

**4. 指定数据根目录：**

```bash
python core/geo_downloader.py --gse GSE107618 --data-root /mnt/e/data
```

**5. 强制重新下载：**

```bash
python core/geo_downloader.py --gse GSE107618 --force
```

默认情况下已存在的文件会跳过（大小匹配时）。`--force` 让所有文件重新下载。

**6. 多模态数据集下载（如 GSE310245 同时含 RNA 和 ATAC）：**

```bash
python core/geo_downloader.py --gse GSE310245
```

下载器不区分模态，它会下载 FTP 上所有补充文件。模态识别由后续的 `preprocessor.py` 自动完成。

> **提示**：在大文件下载前始终先用 `--dry-run` 预览，确认文件大小和总量再开始下载。

---

## 5. 元数据解析

下载器会从 NCBI FTP 获取 `{GSE}_family.soft.gz` 文件并解析为结构化信息。

### 5.1 SOFT 文件提供了哪些信息

| 字段 | 说明 | 示例 |
|------|------|------|
| `title` | 数据集标题 | "Dissecting the transcriptome landscape of human neural retina" |
| `gse_id` | GEO 编号 | GSE107618 |
| `organism` | 物种（从样本块汇总） | Homo sapiens |
| `n_samples` | 样本数 | 64 |
| `pmid` | 关联 PMID（列表） | ["31269016"] |
| `summary` | 数据集简介 | （SOFT Summary 全文） |
| `series_type` | 数据类型 | "Expression profiling by high throughput sequencing" |
| `platform_title` | 测序平台 | "Illumina NovaSeq 6000" |
| `contributors` | 提交者列表 | ["Doe, John", ...] |
| `submission_date` | 提交日期 | "May 15 2019" |
| `is_superseries` | 是否为 SuperSeries | false |
| `sample_list` | 每个样本的 accession/title/organism | [{"accession": "GSM...", "title": "..."}, ...] |

### 5.2 物种归一化

下载器从每个样本的 `!Sample_organism_ch1` 字段汇总物种信息。如果所有样本属于同一物种，则直接使用该值。如果样本跨物种（例如跨物种比较研究），则用逗号拼接多个物种名称。

归一化后的物种值会写入 Registry 的 `species` 字段，并转换为短标识 slug（`Homo sapiens` → `human`，`Mus musculus` → `mouse`）。

### 5.3 元数据缓存

下载的 SOFT 数据会保存为 `$FUXI_DATA_ROOT/GSEXXXXXX/.geo_meta.json`。后续运行 `--skip-soft` 时将直接从缓存读取，不再请求 NCBI。Registry 的字段填充步骤也会读取该缓存文件。

---

## 6. 补充文件列表

在下载之前，下载器会先读取 NCBI FTP 上的 `suppl/` 目录，列出所有可供下载的文件。

### 6.1 文件发现机制

下载器通过以下步骤发现补充文件：

1. 连接 NCBI 的 FTP 服务器（`ftp.ncbi.nlm.nih.gov/geo/series/.../suppl/`），获取目录列表（HTML 格式）。
2. 解析每个 `<a href="...">` 标签，提取文件名、日期和大小。
3. 过滤掉目录导航链接（`..`、`Parent Directory`）、`.html` 文件和符号链接。
4. 将 NCBI 格式的文件大小（如 `373M`、`1.1G`）转换为字节数和人类可读字符串。

### 6.2 RAW.tar 标记

文件名匹配 `RAW*.tar` 模式（不区分大小写）的文件会用 ★ 标记：

```
    ★ RAW GSE81905_RAW.tar                                           61.0 GB
      filelist.txt                                                 482 B
```

- **★ 标记**：RAW.tar 文件通常包含原始数据（Cell Ranger 输出、FASTQ、BAM 或未处理的 count 矩阵），通常是你最需要下载的文件。
- RAW.tar 文件在输出中始终排在前面，其他文件按文件名排序。
- 无标记的文件是补充表格、元数据文件等辅助文件。

**文件筛选逻辑**：下载器只下载数据文件，跳过目录导航链接、`.html` 文件和符号链接。

**文件大小显示示例**：

| 实际大小 | 显示 |
|----------|------|
| 647 字节 | 647 B |
| 264,799 字节 | 258.6 KB |
| 1,572,864 字节 | 1.5 MB |
| 65,534,590,976 字节 | 61.0 GB |

### 6.3 多文件系列示例

以 GSE235583（24 samples，Visium 空间转录组）为例：

```bash
python core/geo_downloader.py --gse GSE235583 --dry-run
```

输出：

```
============================================================
  GEO Download: GSE235583
============================================================

  [METADATA] GSE235583_family.soft.gz
    Title:    Deciphering the spatio-temporal transcriptional and chromatin
              accessibility of human retinal organoid development at the
    Organism: Homo sapiens
    Platform: Illumina NovaSeq 6000 (Homo sapiens)
    PMID:     none
    Samples:  24

  [SUPPL] Listing files...
    14 file(s) found (163.7 MB total)
    ★ RAW GSE235583_RAW.tar                                         117.0 MB
       filelist.txt                                                8.0 KB
       GSE235583_AD3_D10_counts.csv.gz                             3.5 MB
       GSE235583_AD3_D10_md.csv.gz                               317.0 KB
       GSE235583_AD3_D150_counts.csv.gz                           11.0 MB
       ...
```

### 6.4 没有补充文件时怎么办？

部分 GEO 数据集在 NCBI FTP 上没有补充文件。常见原因包括：

- 数据存放于受控访问仓库（dbGaP、EGA）。
- 数据托管在外部平台（ArrayExpress、Zenodo、实验室自己的服务器）。
- 仅有处理后的数据表格，嵌入在 SOFT 文件中而非独立文件。

此时下载器会提示并正常退出。你可以查看该 GSE 的 GEO 页面（https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXXX），在 "Supplementary data" 或 "Supplementary file" 下寻找数据链接。

---

## 7. 断点续传

下载器原生支持断点续传。

### 7.1 续传原理

- **wget 模式**：使用 `--continue` 参数，自动从上次中断处继续。
- **curl 模式**：使用 `-C -` 参数，行为同 wget 的断点续传。

下载前会检查本地已有文件的大小是否与 NCBI FTP 上的一致：

- **大小匹配** → 跳过（`[SKIP]`）
- **大小不匹配** → 自动续传（`[RESUME]`），显示已下载部分和总大小

续传示例：

```
    [RESUME] GSE81905_RAW.tar (partial: 4.2 GB / 61.0 GB)
```

如果文件存在但小于预期，则从断点继续。如果传入 `--force`，则无论文件是否存在都重新下载。

### 7.2 安全中断（Ctrl+C）

下载过程中随时可以按 `Ctrl+C` 中断，之后重新运行相同的命令即可从断点继续。每个文件的进度独立保存，已完整下载的文件不会被重复下载。

### 7.3 超时设置

下载器不设文件下载超时——wget 和 curl 自带网络超时机制（连接超时、读取超时），由它们根据网络状况自行处理。如果网络长时间无响应，wget/curl 会自动重试或报错退出，此时重新运行命令即可续传，已下载的部分不会浪费。

### 7.4 大文件下载建议

对于 RAW.tar 超过 10GB 的数据集（如 GSE81905 的 61GB RAW.tar）：

```bash
# 使用 screen / tmux 让下载在后台持续运行
tmux new -s geo_download
python core/geo_downloader.py --gse GSE81905
# Ctrl+B, D 分离，随时回来检查
tmux attach -t geo_download
```

---

## 8. Registry 联动

下载完成后，下载器会自动更新 **Master Registry**（`core/paper/registry.py` 管理的统一登记表）。

### 8.1 自动联动路径

下载成功后，更新流程如下：

```
下载完成（零失败）
       ↓
update_registry_after_download()
       ↓
1. 数据集状态 → "data_downloaded"
2. enrich_dataset_from_soft()
   ├─ species       （如果为空，从 SOFT organism 填充）
   ├─ n_samples     （如果为 None，从 SOFT 样本数填充）
   ├─ paper_pmids   （追加不重复的 PubMed ID）
   ├─ links         （自动创建 PaperDatasetLink 关联已有论文）
   └─ notes         （如果为空，填充 SOFT summary 前 200 字）
```

### 8.2 空字段优先策略

字段填充遵循"仅填充空字段"原则，不会覆盖手工维护的数据。例如，如果你已经在 Registry 中设置了 `species: human`，SOFT 中的物种值会被忽略。

| Registry 字段 | SOFT 来源 | 填充逻辑 |
|--------------|-----------|----------|
| `species` | `organism` | 标准化为短标识 slug（`Homo sapiens` → `human`，`Mus musculus` → `mouse`） |
| `n_samples` | 样本块计数 | 仅当原来为 `None` 时填充 |
| `paper_pmids` | `Series_pubmed_id` | 追加不重复的 PMID |
| `links` | PMID 匹配已有论文 | 自动创建 `PaperDatasetLink`（role=PRIMARY） |
| `notes` | `summary` 前 200 字 | 仅当 notes 为空时填充 |

此外，`update_registry_after_download()` 还会将数据集的状态更新为 `DATA_DOWNLOADED`。

### 8.3 CLI 触发方式

**独立 CLI**：运行完成后（零失败）自动调用 `update_registry_after_download()`，无需额外参数。

**预处理脚本**：使用 `--download` 参数时，Registry 更新作为流程的一部分自动执行。数据集先标记为已下载，然后预处理脚本进入自己的阶段。

**Registry register --pmid**：使用 `--download` 参数时，新增论文并下载数据后，关联关系自动写入 Registry。

### 8.4 查看结果

```bash
# 查看注册表摘要
python -m core.registry status --gse GSE123456      # 查看 GSE 状态（注册、数据、配置）
python -m core.registry report

# 查看特定数据集的注册信息
python -c "
from core.registry import load_master_registry
reg = load_master_registry()
ds = reg.datasets.get('GSE107618')
print(ds.model_dump_json(indent=2))
"
```

---

## 9. 特殊场景

### 9.1 SuperSeries 检测

如果数据集是一个 SuperSeries（包含多个子数据集），下载器会在元数据阶段自动检测并发出警告：

```
    ⚠  SuperSeries detected — individual sub-series should be downloaded separately
```

例如 GSE81905（683 samples，61GB RAW.tar，4 个子系列），下载器会提示它包含了 GSE81906、GSE81907、GSE81908、GSE81909 四个子数据集。

**处理建议**：先下载父系列了解整体结构（父系列通常含一个合并的 RAW.tar），再对感兴趣的子系列分别运行下载器以获取各自的 SOFT 元数据。

### 9.2 5 位 GSE 编号兼容性

Fuxi 的 GEO 下载器已适配所有 GSE 编号格式，包括早期 5 位编号（如 GSE81905）和现代 6 位编号（如 GSE107618）。NCBI 的 FTP 目录按数字前缀分组，5 位编号取前 2 位，6 位编号取前 3 位。

| Accession | 数字部分 | 前缀 (nnn) |
|-----------|----------|------------|
| GSE81905 | 81905 | GSE81nnn |
| GSE107618 | 107618 | GSE107nnn |
| GSE310245 | 310245 | GSE310nnn |

```python
# 内部 URL 构建逻辑
# GSE107618 → 取前 3 位 "107" → GSE107nnn/GSE107618/
# GSE81905  → 取前 2 位 "81"  → GSE81nnn/GSE81905/
# GSE235583 → 取前 3 位 "235" → GSE235nnn/GSE235583/
```

### 9.3 非单细胞数据（bulk RNA-seq、STARR-seq 等）

下载器本身对数据类型不做过滤，它会下载 NCBI FTP 上的所有补充文件。对于非单细胞数据（如 bulk RNA-seq、STARR-seq、ChIP-seq），数据同样会正常下载到 `$FUXI_DATA_ROOT/GSEXXXXXX/`。后续的 `preprocessor.py` 在格式检测时可能无法匹配单细胞模态，此时需要手动指定模态或使用 `--input-dir` 自定义流程。

Registry 中的 `DatasetEntry.non_pipeline` 字段标记了这类数据集（bulk / STARR / SuperSeries 容器不会进入管线流程）。

### 9.4 Visium 空间转录组（GSM 级文件限制）

对于 10X Visium 空间转录组数据（如 GSE235583，24 samples），NCBI FTP 上常见的文件结构是每个 GSM 一个目录：

```
GSE235583/
  suppl/
    GSE235583_RAW.tar           ← 包含所有 SpaceRanger 输出
    GSM7453000_Sample1_image.tif   ← H&E 染色图
    GSM7453001_Sample2_image.tif   ← H&E 染色图
    ...
```

下载器会下载 `RAW.tar` 以及所有 GSM 级的补充文件。如果下载器只发现独立 GSM 文件（如 `GSMXXXXXXX_*.h5`）而找不到 `RAW.tar`，你需要下载每个样本的文件后手动重建 SpaceRanger 输出目录结构，再运行预处理。Visium 数据在预处理阶段需要使用 `--modality spatial` 参数。

### 9.5 多模态（Multiome）数据

对于同时包含 RNA 和 ATAC 的多模态数据集（如 GSE310245，6 samples），下载器会一并下载所有文件。后续 `preprocessor.py` 会自动检测多模态结构并生成 RNA 和 ATAC 两份 config。

### 9.6 多数据集论文过滤

一篇论文可能同时提交 scRNA-seq、ATAC-seq、Hi-C、ChIP-seq 等多种数据。你可以使用 `register --pmid` 的交互模式或 `--datasets` 参数只注册单细胞相关数据集，跳过不适合管线流程的数据。例如：

```bash
# 交互模式：运行后输入序号选择需要的 GSE
python -m core.registry register --pmid 00000000

# 非交互模式：直接指定
python -m core.registry register --pmid 00000000 --datasets GSE123456
```

---

## 10. 常见问题（FAQ）

### Q1: 下载中断了怎么办？

直接重新运行相同的命令即可。

```bash
# 中断后重新运行
python core/geo_downloader.py --gse GSE107618
```

下载器会自动检测已存在的文件：
- 已完整下载且大小匹配 → 跳过
- 部分下载 → 从断点续传
- 大小异常（如截断） → 重新下载

### Q2: 提示 "Neither wget nor curl found" 怎么处理？

安装其中一个：

```bash
# Debian / Ubuntu / WSL
sudo apt install wget        # 推荐（断点续传体验更好）
# 或
sudo apt install curl

# macOS
brew install wget
# 或
brew install curl
```

下载器调用系统工具作为子进程完成实际传输，不提供纯 Python 模式（原生断点续传和进度显示更可靠）。装完后重新运行即可。如果系统限制不能安装，你也可以手动下载文件放到 `$FUXI_DATA_ROOT/GSEXXXXXX/` 下，然后跳过下载步骤直接跑预处理。

### Q3: 能否只下载部分文件？

目前下载器是全量下载（下载 FTP 上的所有补充文件）。如果你只想下载特定文件：

**方案 1**：先用 `--dry-run` 查看文件列表，然后手动用 wget 下载你需要的文件：

```bash
# 查看有哪些文件
python core/geo_downloader.py --gse GSE107618 --dry-run

# 手动下载某个文件
wget --continue \
  -O /data/geo_datasets/GSE107618/GSE107618_RAW.tar \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107618/suppl/GSE107618_RAW.tar
```

GEO FTP 的 URL 模式为 `https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{GSE_ID}/suppl/{filename}`。

**方案 2**：下载完成后删除不需要的文件。

### Q4: 我的数据集是 SuperSeries，应该下载哪个 GSE？

对于 SuperSeries（如 GSE81905），你有两种选择：

**选项 A：下载父系列（元数据 + 文件索引）**

```bash
python core/geo_downloader.py --gse GSE81905
```

这会下载父系列的 SOFT 元数据和补充文件。父系列通常有一个合并的 RAW.tar 包含所有子系列的原始数据（GSE81905 的 RAW.tar 达 61GB）。

**选项 B：分别下载各子系列**

查看 `--dry-run` 输出的元数据中列出的子系列编号，然后分别下载：

```bash
python core/geo_downloader.py --gse GSE81906 --dry-run
python core/geo_downloader.py --gse GSE81907 --dry-run
# 选择合适的子系列下载
python core/geo_downloader.py --gse GSE81906
```

**推荐**：先下载父系列获取合并数据，再根据实验设计下载需要的子系列以获取各自的 SOFT 元数据。

### Q5: 下载完了下一步是什么？

下载完成后，数据存放在 `$FUXI_DATA_ROOT/GSEXXXXXX/`。目录结构如下：

```
$FUXI_DATA_ROOT/
└── GSE107618/
    ├── .geo_meta.json           # SOFT 元数据缓存（隐藏文件）
    ├── GSE107618_RAW.tar        # 原始数据（可能需要解压）
    └── GSE107618_Supplementary_Data.xlsx  # 补充表格
```

如果 RAW.tar 包含多个样本的数据，预处理脚本的 Phase 1（归档解压）会将其解压到子目录中：

```
$FUXI_DATA_ROOT/GSE107618/
    ├── .geo_meta.json
    ├── GSE107618_RAW.tar
    ├── sample1_filtered_feature_bc_matrix.h5
    ├── sample2_filtered_feature_bc_matrix.h5
    └── ...
```

接下来运行预处理脚本生成 pipeline 配置文件：

```bash
# 自动检测模态并生成 config
python core/preprocess/preprocessor.py --gse GSE107618

# 如果是 Visium 空间数据，指定空间模态
python core/preprocess/preprocessor.py --gse GSE235583 --modality spatial
```

预处理会完成：
1. 解压归档文件
2. 检测文件格式（10X MTX / HDF5 / Visium / 等）
3. 推断模态（scRNA-seq / scATAC-seq / Spatial）
4. 生成 `dataset.yaml` 和 `config_GSEXXXXXX.yaml`

然后你就可以运行完整 Pipeline 了：

```bash
python core/run_pipeline.py --modality rna --config projects/rna/GSE107618/config_GSE107618.yaml
```

### Q6: 下载时提示 "No supplementary files found" 怎么办？

部分 GEO 数据集在 NCBI FTP 上没有补充文件。常见原因：

- 数据存放于受控访问仓库（dbGaP、EGA）。
- 数据托管在外部平台（ArrayExpress、Zenodo、实验室自己的服务器）。
- 仅有处理后的数据表格，嵌入在 SOFT 文件中而非独立文件。

请查看该 GSE 的 GEO 页面（https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXXX），在 "Supplementary data" 或 "Supplementary file" 下寻找数据链接。

### Q7: FUXI_DATA_ROOT 没设置怎么办？

```bash
# 方法 1：设置环境变量（推荐）
export FUXI_DATA_ROOT=/data/geo_datasets

# 方法 2：每次运行时指定
python core/geo_downloader.py --gse GSE107618 --data-root /data/geo_datasets

# 方法 3：写到 .env 文件（如果项目使用 .env）
echo "export FUXI_DATA_ROOT=/data/geo_datasets" >> .env
source .env
```

### Q8: 能否批量下载多个数据集？

```bash
# 在 bash 中循环
for gse in GSE107618 GSE118614 GSE235583; do
    python core/geo_downloader.py --gse "$gse"
done
```

或者使用 dry-run 先查看所有数据集的信息：

```bash
for gse in GSE107618 GSE118614 GSE235583; do
    echo "----- $gse -----"
    python core/geo_downloader.py --gse "$gse" --dry-run --quiet
done
```

### Q9: 工具标记了 SuperSeries，但它有自己的文件，还需要下载子系列吗？

是的。父 SuperSeries 可能有一个合并的 RAW.tar 包含所有原始数据，但每个子系列有自己独立的 SOFT 元数据，包含样本级别的详细信息（组织、处理条件、批次等）。预处理脚本和下游分析工具依赖子系列的元数据才能正确工作。

建议的做法：下载父系列获取合并的 RAW.tar，同时分别下载各子系列获取各自的 SOFT 元数据。
