# Bulk RNA-seq 入门与流程概览

> 单细胞告诉我们"谁在唱"，Bulk 告诉我们"合唱声有多大"。

## 什么是 Bulk RNA-seq？

Bulk RNA-seq 是最经典、最广泛使用的转录组学方法。它测量的是**整个组织样本或细胞群体**的基因表达总和，而不是单个细胞。

相比单细胞 RNA-seq：
- **成本低得多** —— 每个样本几十到几百美元，而非几千美元
- **通量高** —— 可以轻松比较几十个条件下的数百个样本
- **统计功效强** —— 样本间变异（biological replicates）提供了真正的生物学结论基础
- **但丢失了细胞类型信息** —— 你看到的表达值是所有细胞的平均值

## 什么时候用 Bulk？

- **药物筛选** —— 比较药物处理组 vs 对照组
- **基因敲除/过表达** —— 比较 KO vs WT 的全局转录变化
- **疾病 vs 对照** —— 寻找疾病相关的差异表达基因
- **时间序列** —— 连续时间点的基因表达变化（虽然后续分析需要额外配置）

## Fuxi Bulk 管线概览

Fuxi 的 Bulk 管线只有 **5 个核心步骤 + 1 个可选步骤**，比单细胞管线短很多：

```
计数矩阵 → 00_load → 01_qc → 02_de → 03_enrichment → 04_exploratory
                                            (可选 05_batch)
```

### Step 00 — 加载数据

读入计数矩阵（CSV/TSV）或已有的 h5ad 文件。支持三种格式：
- `count_matrix`：原始整数计数矩阵
- `tpm_matrix`：TPM/FPKM 表达矩阵
- `h5ad`：已有的 AnnData 格式

### Step 01 — 质量控制

**样本级别**的 QC（记住，Bulk 没有"细胞"这个概念）：
- 文库大小（每个样本的总 counts）
- 基因检出率
- 样本间相关性热图
- 标记不合格样本

### Step 02 — 差异表达（核心步骤）

使用 **PyDESeq2**（纯 Python 实现的 DESeq2）：
1. **标准化** —— DESeq2 的"中位数比率"法（median of ratios）
2. **分散度估计** —— 负二项模型的基因级分散度
3. **Wald 检验** —— 比较处理组 vs 对照组
4. **LFC 收缩** —— 减少低表达基因的假阳性（apeGLM 方法）

输出：
- `02_de_results.csv`：所有基因的 DE 结果
- `02_de_significant.csv`：显著差异表达基因（padj < 0.05）
- `02_volcano.png`：火山图
- `02_ma_plot.png`：MA 图

### Step 03 — 通路富集

使用 **GSEApy** 进行两种分析：
- **Over-Representation Analysis (ORA)**：对上下调基因分别做超几何检验
- **Preranked GSEA**：对所有基因按表达变化排序，检测整体趋势

### Step 04 — 探索性可视化

生成标准的 Bulk 分析图表：
- PCA 图（按实验条件着色）
- 样本间距离热图
- 差异表达基因热图（top 50）
- 表达箱线图（top 10 DEGs）

### Step 05 — 批次校正（可选）

如果数据来自多个批次（例如不同日期测序的样本），可用 **ComBat**（pycombat）去除批次效应。

> **注意**：Step 05 是独立的可选步骤。如果只有一个批次，直接跳过不会影响下游分析。

## 为什么用 PyDESeq2 而不是 R 的 DESeq2？

- **纯 Python**：不需要装 R、Bioconductor、rpy2，部署简单
- **API 友好**：Pythonic 的接口，与 Fuxi 的其余部分无缝集成
- **性能可靠**：SciPy 社区的维护项目，已有 1500+ GitHub stars
- **结果等效**：统计学上与原版 DESeq2 结果一致

## 与 scRNA-seq 的 Fuxi 配置区别

如果你熟悉 Fuxi 的 RNA 管线，Bulk 管线有以下简化：

| 不需要的配置段 | 原因 |
|----------------|------|
| `scrublet` | 无双细胞检测 |
| `hvg` | 所有基因都用于 DE |
| `clustering` | 无聚类 |
| `harmony` | 无批次整合（用 ComBat 替代） |
| `marker` | 无细胞注释 |
| `trajectory` | 无轨迹推断 |
| `cci` | 无细胞互作 |
| `grn` | 无 GRN 分析 |

Bulk 管线只需要：`data_input` + `bulk` + `enrichment` + `execution`

## 快速开始

```bash
# 1. 安装依赖
pip install -r requirements/bulk.txt

# 2. 查看步骤
python core/run_pipeline.py --modality bulk --list

# 3. 运行管线
python core/run_pipeline.py --modality bulk --config projects/bulk/<GSE_ID>/config_<GSE_ID>.yaml

# 4. 断点重跑
python core/run_pipeline.py --modality bulk --resume --config ...
```

> 完整使用指南参见 [Pipeline 用户指南](../pipeline_guide_zh-CN.md)
