# Fuxi 生产部署指南

> 适用于：**只需运行管线、不需要开发代码**的使用者
>
> 如果你需要修改代码或运行测试，请参阅[环境配置指南](environment_setup_zh-CN.md)。

---

## 目录

1. [这是什么？](#1-这是什么)
2. [前置条件](#2-前置条件)
3. [首次部署（一键）](#3-首次部署一键)
4. [配置环境](#4-配置环境)
5. [验证安装](#5-验证安装)
6. [运行管线](#6-运行管线)
7. [日常更新（从开发机同步）](#7-日常更新从开发机同步)
8. [目录结构说明](#8-目录结构说明)
9. [常见问题（FAQ）](#9-常见问题faq)

---

## 1. 这是什么？

Fuxi 支持两种使用方式：

| 方式 | 适用人群 | 说明 |
|------|----------|------|
| **开发环境** | 代码贡献者 | 完整 git 仓库 + 测试 + linter，用 `pip install -e .` |
| **生产环境** ← 本指南 | 管线使用者 | 只含运行时代码 + 依赖，用 `bin/bootstrap-prod.sh` 一键搭建 |

本指南面向**生产环境**使用者：你拿到了一份 Fuxi 源码（通过 rsync、scp 或 U 盘），想用最少的步骤把它跑起来。

---

## 2. 前置条件

| 需求 | 说明 |
|------|------|
| **操作系统** | Linux 或 WSL2（Windows 用户） |
| **Python** | 3.14+ |
| **磁盘空间** | ≥ 50 GB（原始数据 + 中间结果） |
| **内存** | ≥ 16 GB（ATAC-seq 建议 32 GB+） |
| **Fuxi 源码** | 已通过 rsync / scp / U 盘获取 |

> ⚠️ **Windows 用户**：scATAC-seq 依赖 Snapatac2，仅支持 Linux。请使用 WSL2。

---

## 3. 首次部署（一键）

### 3.1 默认路径

如果你的布局是标准的（源码在 `/mnt/e/fuxi-prod`，venv 在 `~/.local/venvs/fuxi-prod`）：

```bash
# 进入源码目录（或开发机同步过来的目录）
cd /path/to/fuxi-source

# 一键部署
bin/bootstrap-prod.sh
```

脚本会自动完成 4 步：

1. **同步源码** → 到生产目录（rsync，跳过已存在文件）
2. **创建虚拟环境** → 在 Linux 原生文件系统上（`~/.local/venvs/fuxi-prod`）
3. **安装运行时依赖** → `pip install .[all]`（约 3-10 分钟，复用 pip 缓存）
4. **生成配置模板** → `.env`、`global.yaml`、`projects/` 目录骨架

### 3.2 自定义路径

```bash
FUXI_PROD_DIR=/opt/fuxi \
FUXI_PROD_VENV=/opt/fuxi-venv \
bin/bootstrap-prod.sh
```

### 3.3 完成后你会看到

```
[bootstrap] => Production environment ready.

  Source:   /mnt/e/fuxi-prod
  Venv:     /home/user/.local/venvs/fuxi-prod

  Next steps:
    1. Edit  /mnt/e/fuxi-prod/.env          (FUXI_DATA_ROOT, LLM_API_KEY, ...)
    2. Edit  /mnt/e/fuxi-prod/global.yaml   (machine-specific settings)
    ...
```

---

## 4. 配置环境

### 4.1 必须编辑：`.env`

```bash
vim /mnt/e/fuxi-prod/.env
```

关键字段：

```ini
# 数据根目录（必填）—— 存放 GEO 数据集的目录
FUXI_DATA_ROOT=/mnt/e/data

# AI 注释 API 密钥（可选）
LLM_API_KEY=sk-your-key-here

# WSL2 用户必需
HDF5_USE_FILE_LOCKING=FALSE
```

### 4.2 按需编辑：`global.yaml`

```bash
vim /mnt/e/fuxi-prod/global.yaml
```

包含执行参数（CPU 核数、内存策略）、聚类参数、可视化设置等。大多数情况下默认值即可使用，无需修改。

### 4.3 添加数据集配置

将你的数据集配置放到对应模态目录：

```bash
# RNA 数据集
projects/rna/GSE123456/config_GSE123456.yaml

# ATAC 数据集
projects/atac/GSE123456/config_GSE123456.yaml
```

> 📖 配置文件格式请参阅[预处理脚本使用指南](preprocessor_guide_zh-CN.md)和[配置模板](../templates/config_templates/)。

---

## 5. 验证安装

```bash
VENV=~/.local/venvs/fuxi-prod
TARGET=/mnt/e/fuxi-prod

# 检查管线是否可启动
$VENV/bin/python $TARGET/core/run_pipeline.py --modality rna --list
```

如果看到步骤列表（00-12），说明环境就绪。

---

## 6. 运行管线

### 6.1 基本命令

```bash
# 定义便捷变量（可选，加到 ~/.bashrc）
export FUXI_PYTHON=~/.local/venvs/fuxi-prod/bin/python
export FUXI_HOME=/mnt/e/fuxi-prod

# 列出步骤
$FUXI_PYTHON $FUXI_HOME/core/run_pipeline.py --modality rna --list

# 运行完整管线
$FUXI_PYTHON $FUXI_HOME/core/run_pipeline.py \
    --modality rna \
    --config $FUXI_HOME/projects/rna/GSE123456/config_GSE123456.yaml

# 运行单步
$FUXI_PYTHON $FUXI_HOME/core/run_pipeline.py \
    --modality rna --step 3 \
    --config $FUXI_HOME/projects/rna/GSE123456/config_GSE123456.yaml

# 从断点恢复
$FUXI_PYTHON $FUXI_HOME/core/run_pipeline.py \
    --modality rna --resume \
    --config $FUXI_HOME/projects/rna/GSE123456/config_GSE123456.yaml
```

### 6.2 设置环境变量（持久化）

将以下内容追加到 `~/.bashrc`：

```bash
# ── Fuxi 生产环境 ──
export FUXI_DATA_ROOT=/mnt/e/data
export HDF5_USE_FILE_LOCKING=FALSE
# export LLM_API_KEY=sk-...

# 便捷别名
alias fuxi='~/.local/venvs/fuxi-prod/bin/python /mnt/e/fuxi-prod/core/run_pipeline.py'
```

之后可以直接：

```bash
fuxi --modality rna --list
```

---

## 7. 日常更新（从开发机同步）

当开发机上的代码有更新时，**在开发机上**运行：

```bash
# 在开发机上
cd /mnt/d/Projects/Fuxi
bin/deploy.sh
```

`deploy.sh` 会：

1. 用 `git describe` 生成版本标记（`version.txt`）
2. rsync 同步变更的源码文件到生产目录（增量同步，秒级完成）
3. 在生产 venv 中重新 `pip install`（更新 fuxi 包，复用缓存）

> 💡 只想同步源码、不重装依赖？用 `bin/deploy.sh --no-reinstall`
>
> 💡 想预览会同步哪些文件？用 `bin/deploy.sh --dry-run`

### 更新后验证版本

```bash
cat /mnt/e/fuxi-prod/version.txt
# 输出示例: v0.2.0-5-g3a8f2b1
```

---

## 8. 目录结构说明

```
开发机                                    生产机
/mnt/d/Projects/Fuxi/                    /mnt/e/fuxi-prod/
├── .git/                    不同步 →    ├── core/, rna/, atac/...  ← rsync 同步
├── .venv/                   不同步      ├── bin/                   ← rsync 同步
├── tests/                   不同步      ├── pyproject.toml         ← rsync 同步
├── core/, rna/, atac/...   ──rsync──→  ├── .env                   ← 生产专属
├── bin/deploy.sh            ──rsync──→  ├── global.yaml            ← 生产专属
├── bin/bootstrap-prod.sh   ──rsync──→  ├── projects/              ← 生产专属
└── .rsync-filter                        ├── results/, logs/         ← 运行时生成
                                         └── version.txt            ← 部署时生成

                                         ~/.local/venvs/fuxi-prod/  ← Python venv
                                         (Linux 原生文件系统)
```

### 为什么 venv 不放在生产目录里？

WSL2 通过 9p 协议访问 `/mnt/*`（Windows 挂载盘）。Python 有 10 万+ 小文件
（`site-packages/`），9p 对大量小文件的 I/O 极慢——pip install 可能耗时
30 分钟以上。放在 Linux 原生文件系统（ext4）上只需 3-5 分钟。

| venv 位置 | pip install 耗时 | 日常运行 |
|-----------|:----------------:|:--------:|
| `/mnt/e/fuxi-prod/.venv/` (Windows 盘) | 30+ 分钟 | 慢 |
| `~/.local/venvs/fuxi-prod/` (Linux fs) | 3-5 分钟 | 正常 |

---

## 9. 常见问题（FAQ）

### Q1: `pip install` 报 "No matching distribution found" 或版本冲突

**原因**：Python 3.14 很新，部分科学计算包（如 liana、numba）在 PyPI 上
声明的兼容版本范围是 `<3.14`，但实际可以正常运行。

**解决**：`bin/bootstrap-prod.sh` 和 `bin/deploy.sh` 已自动添加
`--ignore-requires-python` 参数。如果你手动执行 pip install，请加上此参数：

```bash
pip install --ignore-requires-python .[all]
```

### Q2: rsync 或 pip install 非常慢（WSL2）

**原因**：venv 或生产目录在 `/mnt/*`（Windows 挂载盘）上。

**解决**：确保 venv 在 Linux 原生文件系统上。检查：

```bash
echo $FUXI_PROD_VENV
# 应输出 ~/.local/venvs/fuxi-prod 或其他 Linux 路径，而非 /mnt/...
```

如果已经在 /mnt 上创建了 venv，删除后重新运行 bootstrap：

```bash
rm -rf /mnt/e/fuxi-prod/.venv
bin/bootstrap-prod.sh
```

### Q3: 如何查看当前生产环境跑的是什么版本？

```bash
cat /mnt/e/fuxi-prod/version.txt
```

输出格式为 git describe（如 `v0.2.0-5-g3a8f2b1-dirty`），可与开发机的
`git log --oneline -1` 对照。

### Q4: 生产环境可以自己修改代码吗？

可以，但不推荐。生产环境没有 `.git/`，修改不会被版本管理追踪。
最佳做法是：在开发机上修改 → `bin/deploy.sh` 同步 → 在生产环境运行。

### Q5: 如何更新所有依赖包（不只是 fuxi 本身）？

```bash
# 在开发机上
bin/deploy.sh
# ↑ deploy.sh 默认会重新 pip install，自动拉取新增依赖
```

如果只改了代码、没加新依赖，用 `--no-reinstall` 跳过安装步骤。

### Q6: 可以不用 deploy.sh，手动 rsync 吗？

可以。deploy.sh 本质上就是 rsync + pip install：

```bash
rsync -av --delete \
    --filter='merge .rsync-filter' \
    /mnt/d/Projects/Fuxi/ /mnt/e/fuxi-prod/

~/.local/venvs/fuxi-prod/bin/pip install --ignore-requires-python /mnt/e/fuxi-prod[all]
```

但 deploy.sh 还帮你处理了版本标记和错误检查，推荐使用脚本。

### Q7: 生产环境需要哪些额外 extras？

`bin/bootstrap-prod.sh` 默认安装 `.[all]`（全部模态，不含方法学）。
如需方法学包（CellTypist、scVI 等），在生产环境手动追加：

```bash
~/.local/venvs/fuxi-prod/bin/pip install --ignore-requires-python /mnt/e/fuxi-prod[methods]
```

---

> 📖 **更多文档**：
> - [环境配置指南（开发用）](environment_setup_zh-CN.md)
> - [管线使用指南](pipeline_guide_zh-CN.md)
> - [预处理脚本指南](preprocessor_guide_zh-CN.md)
