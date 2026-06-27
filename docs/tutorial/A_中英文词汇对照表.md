# 中英文词汇对照表

> 本表收录 Fuxi 管道和单细胞组学领域的核心技术、方法和概念，按主题分组编排。

## 分子生物学基础

| 中文 | English | 简要说明 |
|------|---------|---------|
| 中心法则 | central dogma | DNA→RNA→蛋白质的遗传信息流 |
| 转录组 | transcriptome | 一个细胞或组织中所有 mRNA 分子的集合 |
| 基因组 | genome | 一个细胞中所有 DNA 序列的集合 |
| 表观基因组 | epigenome | DNA 甲基化、组蛋白修饰等可遗传的染色质修饰全貌 |
| 基因表达 | gene expression | 基因被转录为 mRNA 的过程（通常用 mRNA 丰度量化） |
| 表达矩阵 | expression matrix | 行×列（细胞×基因）的数值表格，记录每个基因在每个细胞中的表达量 |
| 计数矩阵 | count matrix | 原始整数形式的表达矩阵，每个值是 mRNA 分子计数 |
| 信使RNA | messenger RNA (mRNA) | 携带蛋白质编码信息的 RNA，从细胞核运送到核糖体翻译 |
| 互补DNA | complementary DNA (cDNA) | 由 mRNA 逆转录而来的人工 DNA，用于建库测序 |
| 开放染色质 | open chromatin | DNA 松散缠绕、转录因子可结合的染色质状态 |
| 异染色质 | heterochromatin | DNA 高度压缩、转录机器无法进入的染色质状态 |
| 调控元件 | regulatory element | 控制基因表达的基因组区域（启动子、增强子、沉默子等） |
| 启动子 | promoter | 基因转录起始点上游的调控序列，RNA 聚合酶结合位置 |
| 增强子 | enhancer | 可远距离增强靶基因转录的调控序列 |
| 转录因子 | transcription factor (TF) | 结合特定 DNA 序列、调控基因转录的蛋白质 |
| 基序 | motif | 转录因子偏好的短的、保守的 DNA 结合序列模式 |
| 调控网络 | regulatory network / regulon | 转录因子及其靶基因组成的调控关系图 |

## 单细胞测序技术

| 中文 | English | 简要说明 |
|------|---------|---------|
| 单细胞RNA测序 | single-cell RNA sequencing (scRNA-seq) | 在单个细胞分辨率下测量基因表达 |
| 单细胞ATAC测序 | single-cell ATAC sequencing (scATAC-seq) | 在单个细胞分辨率下测量染色质可及性 |
| 空间转录组学 | spatial transcriptomics | 保留组织空间位置信息的基因表达测量 |
| 转座酶可直接染色的染色质分析 | Assay for Transposase-Accessible Chromatin using sequencing (ATAC-seq) | 用 Tn5 转座酶标记并识别开放染色质区域 |
| Tn5 转座酶 | Tn5 transposase | 一种在开放染色质区域切割 DNA 并连接接头的酶 |
| 峰 | peak | ATAC-seq 中信号显著高于背景的基因组区域 |
| 峰调用 | peak calling | 从测序信号中识别峰的统计方法 |
| 峰×细胞矩阵 | peak-by-cell matrix | ATAC 数据的基本形式：每行是一个峰，每列是一个细胞，值是切割计数 |
| 基序富集 | motif enrichment | 在一组峰区域中搜索特定转录因子结合基序是否出现过频 |
| 10X Genomics | 10X Genomics | 最主流的单细胞测序平台提供商 |
| 液滴 | droplet | 微流控系统中包裹单个细胞和试剂的油包水小滴 |
| 条形码 | barcode | 一段已知的短 DNA 序列，用于标记分子或细胞的身份 |
| 唯一分子标识符 | unique molecular identifier (UMI) | 一段随机序列，用于标记并校正 PCR 重复的单个 mRNA 分子 |
| 空间条形码 | spatial barcode | 编码空间位置信息的 DNA 序列，用于 Visium 等空间技术 |
| 捕获点 | spot | Visium 载玻片上的固定位置点，直径 55μm，捕获组织 mRNA |
| 组织切片 | tissue section | 切成薄片（通常 5-10μm）的组织样品 |
| 双细胞 | doublet | 一个液滴中捕获了两个细胞，是单细胞数据中的常见污染物 |
| H&E 染色 | Hematoxylin and Eosin (H&E) stain | 组织学标准染色，细胞核呈紫色，细胞质呈粉红色 |
| 扩散 | diffusion | mRNA 从组织扩散到捕获阵列上的过程 |
| 测序深度 | sequencing depth | 一个文库或样本的总测序读段数，影响检测灵敏度 |
| 高通量 | high-throughput | 一次实验可同时测量大量分子或细胞的特性 |
| 混样测序 | bulk sequencing | 将一群细胞一起裂解测序，丢失单细胞分辨率 |
| 10X 兼容 HDF5 | 10X HDF5 (`.h5`) | 10X Cell Ranger 输出的压缩格式，包含表达矩阵和元数据 |
| Cell Ranger | Cell Ranger | 10X Genomics 的官方数据预处理软件 |


## 生物信息学基础

| 中文 | English | 简要说明 |
|------|---------|---------|
| 归一化 | normalization | 消除技术差异（如测序深度）使细胞间可比较 |
| 标准化 | standardization / scaling | 使每个基因的表达量在细胞间均值为 0、方差为 1；用于降维前 |
| 批次效应 | batch effect | 由测序日期、操作人员、试剂批次等技术因素引入的系统性差异 |
| 整合 | integration | 消除批次效应并将不同批次的数据对齐到共同空间 |
| 高变基因 | highly variable genes (HVG) | 在细胞间表达变异最大的基因，常用于降维前的特征筛选 |
| 降维 | dimensionality reduction | 将高维数据（数万个基因）压缩到便于可视化和计算的低维空间 |
| 稀疏矩阵 | sparse matrix | 只存储非零元素及其位置的矩阵格式，节省内存 |
| 检查点 | checkpoint | 管线的中间输出文件，断点续跑的依据 |
| 复现性 | reproducibility | 在不同时间/条件下运行相同分析得到相同结果的能力 |
| 伪复制 | pseudo-replication | 因不了解数据内在结构（如同一样本多次测量）而产生的统计假象 |


## 核心技术方法

| 中文 | English | 简要说明 |
|------|---------|---------|
| 主成分分析 | Principal Component Analysis (PCA) | 寻找数据方差最大的方向，是最经典的线性降维方法 |
| 均匀流形近似与投影 | Uniform Manifold Approximation and Projection (UMAP) | 非线性降维可视化方法，保留数据的全局和局部结构 |
| 莱顿聚类 | Leiden clustering | 社区检测算法，通过优化模块度来寻找图中的细胞群体 |
| Lumapal 整合方法 | Harmony | 用模糊聚类在 PCA 空间中校正批次效应的整合算法 |
| 逆文档频率 | Inverse Document Frequency (IDF) | 给稀有峰更高权重的特征选择方法，类似 TF-IDF 中的 IDF 部分 |
| 谱嵌入 | spectral embedding | 基于图拉普拉斯算子的降维方法，适合二值化/稀疏数据的结构发现 |
| 兰佐斯算法 | Lanczos algorithm | 大规模稀疏矩阵特征值分解的迭代近似算法，谱嵌入的底层方法 |
| 帕蒂·杰罗因·算法 | Partition-based Graph Abstraction (PAGA) | 在聚类基础上构建全局拓扑图，同时呈现连续分化和离散结构 |
| 扩散伪时间 | Diffusion Pseudotime (DPT) | 基于随机游走扩散距离的伪时间推断，量化细胞"进展程度" |
| 伪时间 | pseudotime | 将横截面单细胞数据沿连续变化过程排列为"时间"轴的计算推断 |
| 蒙特卡洛随机排列 | Moran's I | 衡量基因表达是否具有空间自相关（邻近 spot 是否更相似）的统计量 |
| 空间可变基因 | spatially variable gene (SVG) | 呈现出非随机空间表达模式的基因 |
| 轮廓系数 | silhouette score | 评估聚类质量的指标：衡量一个细胞和自己聚类的匹配度 vs 与其他聚类的分离度 |
| 威尔科克森秩和检验 | Wilcoxon rank-sum test | 非参数统计检验，用于比较两组间的分布差异（RNA 差异表达的标准方法） |
| 斯特林指数 | F1 score / adjusted Rand index (ARI) | 衡量聚类结果与真实标签匹配程度的指标 |
| K 最近邻图 | K-Nearest Neighbors graph (KNN) | 连接每个细胞与其最近的 K 个邻居的图结构，聚类和降维的基础 |
| Scrublet | Scrublet | 基于模拟 doublet 的机器学习方法，用于单细胞数据的双细胞检测 |
| 过表达分析 | Over-Representation Analysis (ORA) | 给定基因列表中某个通路/功能的基因是否显著过多的统计检验 |
| 基因集富集分析 | Gene Set Enrichment Analysis (GSEA) | 基于基因排序而非 p 值阈值的通路富集方法，敏感度更高 |
| 基因本体论 | Gene Ontology (GO) | 标准化的基因功能分类体系（生物过程、分子功能、细胞组分） |
| 京都基因与基因组百科全书 | Kyoto Encyclopedia of Genes and Genomes (KEGG) | 手工整理的生物学通路数据库 |
| CollecTRI | CollecTRI | 转录因子—靶基因调控关系的精选数据库 |
| 转录因子活性 | TF activity | 从转录组数据推断的转录因子的调控活跃程度（而非其 mRNA 表达量） |
| 峰到基因连接 | peak-to-gene linkage | 将可及峰映射到其调控的基因（通常是近邻基因）的方法 |
| 负二项分布 | negative binomial distribution | 单细胞数据建模中最常用的统计分布，刻画计数数据的过离散特性 |
| dropout 现象 | dropout | 低表达基因在单细胞测序中未被检测到的技术性缺失 |


## 数据处理与存储

| 中文 | English | 简要说明 |
|------|---------|---------|
| AnnData | AnnData | 专为单细胞组学设计的 Python 数据结构，将表达矩阵、细胞元数据、基因元数据整合在一起 |
| h5ad | h5ad | AnnData 的序列化文件格式（基于 HDF5） |
| MuData | MuData | 多组学 AnnData 的扩展，用于同时存储 RNA+ATAC 等多模态数据 |
| 压缩稀疏行 | Compressed Sparse Row (CSR) | 常用的稀疏矩阵存储格式，按行压缩，行操作（如按细胞查询）高效 |
| 压缩稀疏列 | Compressed Sparse Column (CSC) | 按列压缩的稀疏矩阵存储格式，列操作（如按基因查询）高效 |
| HDF5 | Hierarchical Data Format version 5 (HDF5) | 科学计算中常用的层级数据存储格式，支持压缩和快速随机访问 |
| H5AD 文件 | HDF5-based AnnData file | Fuxi 管线中，后缀为 `.h5ad` 的文件是贯穿全流程的标准数据格式 |
| TSR 谱系 | fragments.tsv | scATAC-seq 的标准输出文件，记录了每个切割事件（染色体、位置、细胞、副本数） |
| metadata | metadata | 描述数据的数据，如 barcode 列表、基因名称、实验条件等 |
| SpaceRanger | SpaceRanger | 10X Visium 的官方数据预处理软件，输出 spot × 基因矩阵和组织图像 |

## 生物学分析

| 中文 | English | 简要说明 |
|------|---------|---------|
| 细胞类型 | cell type | 具有稳定基因表达程序和功能的细胞类别（如 T 细胞、神经元） |
| 细胞亚型 | cell subtype | 大类别下的细分（如 CD8+ 细胞毒性 T 细胞） |
| 细胞状态 | cell state | 细胞在特定条件下的瞬时功能状态（如活化、增殖、应激） |
| 标记基因 | marker gene | 在特定细胞类型中特异高表达的基因，用于鉴别细胞身份 |
| 标记峰 | marker peak | 在特定细胞类型/聚类中可及性显著更高的峰区域 |
| 差异表达 | differential expression (DE) | 比较两组细胞间的基因表达差异，寻找统计上显著的差异基因 |
| 差异可及性 | differential accessibility | 比较两组细胞间的峰可及性差异 |
| 富集分析 | enrichment analysis | 从基因列表中识别在特定通路或功能中"富集"的生物学主题 |
| 通路 | pathway | 一组共同执行特定生物功能的基因和蛋白质的集合 |
| 基因调控网络 | gene regulatory network (GRN) | 描述转录因子、靶基因及其调控关系的网络模型 |
| 轨迹推断 | trajectory inference | 从横截面单细胞数据重建细胞分化或状态转变的连续路径 |
| 降采样 | downsampling | 从数据中随机抽取部分细胞以减小计算负担 |
| 知识库 | Knowledge Base (KB) | Fuxi 中用于细胞类型注释的手工整理的专家标记基因集合 |
| 伪细胞聚合 | pseudobulk aggregation | 将同一细胞类型/聚类的细胞计数合并成一个"虚拟 bulk 样本" |
| 细胞类型注释 | cell type annotation | 为每个聚类分配有生物学意义的细胞类型标签 |
| 染色质状态注释 | chromatin state annotation | 在 ATAC 数据中确定每个聚类对应的染色质可及性特征和细胞类型 |
| 数据驱动聚类 | leiden clustering | 在 KNN 图上优化模块度的图聚类算法，将细胞划分为不同的群体 |
| 邻域图（空间） | spatial neighborhood graph | 基于空间坐标连接相邻 spot 的图，空间分析的基础 |

## Fuxi 管线相关

| 中文 | English | 简要说明 |
|------|---------|---------|
| 管线 | pipeline | 从原始数据到生物学结论的自动化分析流程 |
| 配置文件 | config file / CFG | 用于设置分析参数和数据集元信息的 Python 文件 |
| 模态 | modality | 测序数据的类型（RNA、ATAC、空间） |
| 步骤注册表 | step registry | `run_pipeline.py` 中定义各模态所有步骤名称、脚本和输出文件的结构 |
| 证据融合 | evidence fusion | Fuxi 注释系统中的五层决策引擎，综合 KB 打分、规则和 AI 输出 |
| 专家规则 | expert rules | 知识库中定义的、按优先级排序的确定性匹配规则 |
| 组织知识库 | tissue Knowledge Base (tissue_kb) | 特定组织的专家标记基因集和规则集合 |
| 标准本体 | StandardOntology | Fuxi 用于标准化细胞类型名称的本体映射系统 |
| AI 缓存 | AI cache | 对 LLM 响应按输入摘要进行 SHA-256 哈希摘要并缓存，避免重复调用计费和等待 |
| 提示词模板 | prompt template | 结构化的 AI 交互模板（RNA 注释、 ATAC 注释、结果解读等） |
| 数据根目录 $(FUXI\_DATA\_ROOT) | data root | 通过环境变量指定的目录，存放所有下载的原始数据集 |
| dataset.yaml | dataset.yaml | 记录数据集元信息的 YAML 文件（格式、物种、组织、实验设计等） |
| 预处理脚本 | preprocessor | 自动下载、检测格式、提取存档并生成配置文件的工具 |
| 模型驱动分析 | expression_type | 配置项，指定输入数据的表达量类型（raw_counts, TPM, FPKM, CPM 等），影响后续归一化和 QC 策略 |
| Stage Map | stage_map | 配置项：将 barcode 后缀映射到实验阶段/时间点（如 D0, D7），用于发育轨迹分析 |
| Sample Map | sample_map | 配置项：将 barcode 后缀映射到样本名称，用于批次校正 |
