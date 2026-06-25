# Fuxi 环境配置指南

> 适用于：**单细胞组学研究人员** | 无需编程背景即可上手

---

## 目录

1. [你需要什么？](#1-你需要什么)
2. [第一步：安装 Python](#2-第一步安装-python)
3. [第二步：创建虚拟环境](#3-第二步创建虚拟环境)
4. [第三步：安装依赖包](#4-第三步安装依赖包)
5. [第四步：配置环境变量](#5-第四步配置环境变量)
6. [验证安装是否成功](#6-验证安装是否成功)
7. [平台特别说明](#7-平台特别说明)
8. [常见问题（FAQ）](#8-常见问题faq)

---

## 1. 你需要什么？

在开始之前，请确认你具备以下条件：

| 需求 | 说明 | 如何获取 |
|------|------|---------|
| **操作系统** | Linux（推荐）、WSL2（Windows 用户）、或 macOS | Windows 10/11 自带 WSL2 支持 |
| **Python** | 3.14 或更高版本 | [python.org](https://www.python.org/downloads/) 或系统包管理器 |
| **pip** | Python 包管理器（随 Python 一起安装） | 通常无需单独安装 |
| **Git** | 用于克隆项目代码 | [git-scm.com](https://git-scm.com/) |
| **磁盘空间** | ≥ 50 GB（用于存储原始数据和中间结果） | — |
| **内存** | ≥ 16 GB（推荐 32 GB 以上用于 ATAC-seq 分析） | — |

> 💡 **Windows 用户请注意**：scATAC-seq 分析依赖 Snapatac2，该库需要 Linux 环境。如果你需要跑 ATAC 数据，**强烈建议使用 WSL2**。纯 scRNA-seq 分析在 Windows 原生 Python 下也可以运行。

---

## 2. 第一步：安装 Python

### 2.1 检查是否已安装

打开终端，输入：

```bash
python3 --version
```

如果输出类似 `Python 3.14.x`，说明已安装且版本满足要求，可以跳到[第二步](#3-第二步创建虚拟环境)。

如果显示 `command not found` 或版本低于 3.14，请参照下方说明安装。

### 2.2 Linux（Ubuntu / Debian）

```bash
# 添加 deadsnakes PPA（提供最新 Python 版本）
sudo apt update
sudo apt install software-properties-common -y
sudo add-apt-repository ppa:deadsnakes/ppa -y
sudo apt update

# 安装 Python 3.14
sudo apt install python3.14 python3.14-venv python3.14-dev -y
```

### 2.3 WSL2（Windows）

在 WSL2 终端中，按上述 Linux 方式安装。如果你还没有安装 WSL2：

```powershell
# 在 Windows PowerShell（管理员）中运行
wsl --install -d Ubuntu-24.04
```

安装完成后重启电脑，打开 "Ubuntu" 应用即可进入 Linux 终端。

### 2.4 macOS

```bash
# 使用 Homebrew
brew install python@3.14
```

---

## 3. 第二步：创建虚拟环境

> 💡 **为什么要用虚拟环境？** 虚拟环境将本项目的依赖包与系统其他 Python 项目隔离开来，避免包版本冲突。这是 Python 项目的最佳实践。

### 3.1 获取项目代码

```bash
git clone <项目仓库地址>
cd <项目目录>
```

### 3.2 创建并激活虚拟环境

```bash
# 在项目根目录下创建虚拟环境（命名为 .venv）
python3.14 -m venv .venv
```

> ⚠️ 请使用 `python3.14` 而非 `python3`，确保使用正确的 Python 版本。

**激活虚拟环境：**

```bash
# Linux / WSL2 / macOS
source .venv/bin/activate
```

激活成功后，终端提示符前会出现 `(.venv)` 标识：

```
(.venv) user@host:~/project$
```

### 3.3 确认环境正确

```bash
# 确认 Python 版本
python --version
# 应输出: Python 3.14.x

# 确认 pip 位于虚拟环境中
which pip
# 应输出: .../项目目录/.venv/bin/pip
```

### 3.4 退出与重新进入

```bash
# 退出虚拟环境
deactivate

# 下次使用时重新激活
cd <项目目录>
source .venv/bin/activate
```

> 💡 每次打开新终端窗口都需要重新 `source .venv/bin/activate`。这是正常的，不是环境出了问题。

---

## 4. 第三步：安装依赖包

### 4.1 基础安装

确保虚拟环境已激活（提示符前有 `(.venv)`），然后：

```bash
pip install -r requirements.txt
```

此命令会安装所有必要的 Python 包，包括：

| 类别 | 主要包 | 用途 |
|------|--------|------|
| **数据处理** | `numpy`, `scipy`, `pandas` | 数值计算、稀疏矩阵、数据表格 |
| **可视化** | `matplotlib` | 绘制 UMAP、热图、小提琴图等 |
| **并行计算** | `joblib` | 加速差异基因分析等步骤 |
| **scRNA-seq** | `scanpy`, `scrublet`, `leidenalg` | 单细胞 RNA 分析核心 |
| **批次校正** | `harmony-pytorch` | 去除样本间批次效应 |
| **功能富集** | `gseapy` | GO / KEGG 通路富集分析 |
| **scATAC-seq** | `snapatac2`, `macs3` | 单细胞 ATAC 分析核心 |
| **多组学整合** | `muon`, `mudata` | RNA + ATAC 联合分析 |
| **AI 注释** | `openai` | 调用大语言模型辅助细胞类型注释 |

### 4.2 安装时间参考

首次安装通常需要 **5-15 分钟**，具体取决于网络速度和机器性能。

### 4.3 验证安装

```bash
# 确认核心包已正确安装
python -c "import scanpy; print('scanpy', scanpy.__version__)"
python -c "import snapatac2; print('snapatac2', snapatac2.__version__)"
python -c "import numpy; print('numpy', numpy.__version__)"
```

如果上述命令均无报错，说明依赖安装成功。

### 4.4 后续更新

当项目代码更新后，可能需要安装新的依赖：

```bash
git pull                    # 拉取最新代码
pip install -r requirements.txt  # 安装可能新增的依赖
```

---

## 5. 第四步：配置环境变量

### 5.1 必需环境变量

| 环境变量 | 是否必需 | 说明 |
|----------|:----:|------|
| `FUXI_DATA_ROOT` | ✅ 必需 | 存放原始下载数据的根目录 |
| `LLM_API_KEY` | ❌ 可选 | AI 注释功能的 API 密钥 |
| `HDF5_USE_FILE_LOCKING` | ⚠️ 建议 | WSL2 用户需设为 `FALSE` |

### 5.2 设置数据根目录

```bash
# 设定为你的数据存放目录（替换为实际路径）
export FUXI_DATA_ROOT=/data/geo_datasets

# 示例：如果你把数据放在 /home/user/geo_data 下
export FUXI_DATA_ROOT=/home/user/geo_data
```

`FUXI_DATA_ROOT` 的目录结构应该是这样的：

```
$FUXI_DATA_ROOT/
├── GSE12345/          # 每个数据集一个文件夹（文件夹名即为数据集 ID）
│   ├── dataset.yaml   # 数据集描述文件
│   ├── *.h5           # 原始 HDF5 文件
│   └── *.csv.gz       # 原始 CSV 文件
├── GSE23456/
│   └── ...
└── ...
```

> 💡 管道运行时会从 `$FUXI_DATA_ROOT/<数据集ID>/` 读取原始数据。

### 5.3 设置 AI API 密钥（可选）

如果你希望使用 AI 辅助细胞类型注释功能：

```bash
export LLM_API_KEY=sk-your-api-key-here
```

> 💡 不设置此变量也可以运行大部分分析步骤。AI 注释是一个可选增强功能。

### 5.4 WSL2 用户额外设置

WSL2 访问 Windows 文件系统时，HDF5 文件锁可能导致错误。请添加：

```bash
export HDF5_USE_FILE_LOCKING=FALSE
```

### 5.5 持久化环境变量

上述 `export` 命令只在当前终端窗口有效。要让它们每次打开终端时自动生效，将以下内容追加到 `~/.bashrc` 文件末尾：

```bash
# ── Fuxi 管道环境变量 ──
export FUXI_DATA_ROOT=/data/geo_datasets    # 替换为你的实际路径
export HDF5_USE_FILE_LOCKING=FALSE           # WSL2 用户需要
# export LLM_API_KEY=sk-...                  # 可选：取消注释并填入你的 API Key
```

保存后运行 `source ~/.bashrc` 使其立即生效。

---

## 6. 验证安装是否成功

完成以上所有步骤后，运行以下命令全面检查环境是否就绪：

```bash
# 1. 确认在项目目录中
cd <项目目录>

# 2. 确认虚拟环境已激活
echo $VIRTUAL_ENV
# 应输出项目目录下的 .venv 路径

# 3. 确认环境变量已设置
echo $FUXI_DATA_ROOT
# 应输出你设置的数据目录路径

# 4. 确认核心功能可用
python -c "
import scanpy as sc
import snapatac2 as snap
import anndata
import numpy as np
import pandas as pd
import matplotlib
print('✅ 所有核心包加载成功')
print(f'  Scanpy:    {sc.__version__}')
print(f'  Snapatac2: {snap.__version__}')
print(f'  AnnData:   {anndata.__version__}')
print(f'  NumPy:     {np.__version__}')
print(f'  Pandas:    {pd.__version__}')
print(f'  Matplotlib:{matplotlib.__version__}')
"
```

如果看到绿色的 `✅ 所有核心包加载成功` 和各包的版本号，说明环境配置完毕。

---

## 7. 平台特别说明

### 7.1 Linux（原生）

✅ **最佳体验**。所有功能（scRNA-seq + scATAC-seq + 多组学整合）完整可用。推荐 Ubuntu 22.04 或更高版本。

### 7.2 WSL2（Windows 用户）

✅ **功能完整**。所有分析与原生 Linux 一致。额外注意事项：

- 数据建议存放在 WSL2 内部（如 `/home/user/data/`），而非 `/mnt/c/` 或 `/mnt/e/` 下的 Windows 文件系统。跨文件系统访问 HDF5 文件会显著变慢。
- 务必设置 `HDF5_USE_FILE_LOCKING=FALSE`。
- 如果数据已在 Windows 磁盘上，可以用 `cp -r /mnt/e/geo_data ~/data/` 复制到 WSL 内部。

### 7.3 macOS

⚠️ **部分功能受限**。Snapatac2 依赖 Linux 底层特性，在 macOS 上可能无法安装或运行。scRNA-seq 分析功能完整可用。macOS 用户如需跑 ATAC 分析，建议使用 Linux 虚拟机或云服务器。

### 7.4 Windows（原生，非 WSL）

⚠️ **仅 scRNA-seq**。Snapatac2 不支持 Windows。如果你只需要跑 scRNA-seq 分析，可以在 Windows 原生 Python 下使用。安装步骤与 Linux 相同，只需将 `source .venv/bin/activate` 替换为：

```cmd
.venv\Scripts\activate
```

---

## 8. 常见问题（FAQ）

### Q1: `pip install -r requirements.txt` 报错，某个包安装失败

常见的包安装问题及解决方案：

**`snapatac2` 安装失败：**

```bash
# 确认使用 Python 3.14+
python --version

# 如果在 macOS 上
# snapatac2 仅支持 Linux，可跳过 ATAC 相关功能：
pip install scanpy anndata numpy scipy pandas matplotlib joblib scrublet leidenalg gseapy openai python-dotenv harmony-pytorch
```

**`harmony-pytorch` 安装失败：**

```bash
# PyTorch 可能需要单独安装
pip install torch
pip install harmony-pytorch
```

**`gseapy` 安装失败：**

```bash
# gseapy 依赖 Rust 编译器
# Ubuntu / Debian:
sudo apt install cargo rustc -y
pip install gseapy
```

**其他包安装失败：**

请检查：
1. 是否已安装 Python 开发头文件：`sudo apt install python3.14-dev`
2. 是否有编译工具链：`sudo apt install build-essential`
3. 网络连接是否正常（部分包需要从 PyPI 下载）

### Q2: 激活虚拟环境后 `python` 命令找不到

```bash
# 确认 .venv 目录存在
ls -la .venv/bin/python*

# 如果不存在，重新创建
python3.14 -m venv .venv --clear
```

### Q3: WSL2 中 "Unable to open file" 或 HDF5 相关错误

```bash
# 原因：HDF5 文件锁在跨文件系统访问时出问题
# 解决：
export HDF5_USE_FILE_LOCKING=FALSE

# 或把数据移到 WSL2 内部
cp -r /mnt/e/geo_data ~/data/
export FUXI_DATA_ROOT=~/data
```

### Q4: 虚拟环境激活后运行 `python` 仍然使用系统 Python

```bash
# 确认 which python 指向 .venv
which python
# 应输出: /path/to/project/.venv/bin/python

# 如果仍然指向系统 Python:
hash -r                  # 清除 bash 命令缓存
source .venv/bin/activate
```

### Q5: 如何彻底删除并重建虚拟环境？

```bash
deactivate                # 先退出
rm -rf .venv              # 删除旧的虚拟环境
python3.14 -m venv .venv  # 重新创建
source .venv/bin/activate
pip install -r requirements.txt
```

### Q6: 内存不足导致安装或运行失败

建议：
- 确保至少有 **16 GB 物理内存**；ATAC-seq 分析建议 **32 GB+**
- 关闭其他占用内存的应用
- 对于超大 scRNA-seq 数据集（>50K 细胞），在配置文件中启用降采样
- WSL2 用户可以在 `%UserProfile%\.wslconfig` 中调整分配给 WSL 的内存上限

### Q7: 需要 GPU 吗？

**不需要。** 管道的所有计算都在 CPU 上运行。有 GPU 可以加速 Harmony（PyTorch 后端）的批次校正步骤，但不是必需的。

### Q8: 可以用 conda / Anaconda 代替 venv 吗？

可以。如果你习惯使用 conda，用以下命令替代 virtualenv 步骤即可，后续 `pip install -r requirements.txt` 步骤不变：

```bash
conda create -n fuxi python=3.14 -y
conda activate fuxi
pip install -r requirements.txt
```

> 💡 本文档以 venv 为主线，因为 conda 渠道对 Python 3.14 的支持可能滞后，且 Snapatac2 推荐 pip 安装。两种方式本质上都是 `pip install`，选你熟悉的即可。

---

## 环境配置检查清单

在开始运行管道之前，请逐项确认：

- [ ] Python 3.14+ 已安装
- [ ] 虚拟环境 `.venv` 已创建并激活
- [ ] `pip install -r requirements.txt` 执行成功
- [ ] `FUXI_DATA_ROOT` 环境变量已设置且目录存在
- [ ] （WSL2 用户）`HDF5_USE_FILE_LOCKING=FALSE` 已设置
- [ ] （可选）`LLM_API_KEY` 已设置（如需 AI 注释功能）
- [ ] 验证脚本运行无误（见[第6节](#6-验证安装是否成功)）

全部完成后，即可参照《预处理脚本使用指南》开始准备数据和运行管道。

---

> 📖 **下一步**：环境准备好后，请阅读《Fuxi 预处理脚本使用指南》了解如何从原始下载数据生成管道配置文件。
