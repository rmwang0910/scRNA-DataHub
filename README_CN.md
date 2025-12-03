# scRNA-DataHub

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-linux%20%7C%20macOS-lightgrey)](https://github.com)

[English](README.md) | 简体中文

**scRNA-DataHub** 是一个通用的单细胞RNA测序数据处理工具，支持 **17种** 数据格式的读取、转换和标准化，统一输出为 H5AD 格式，可直接用于 Scanpy、Seurat 等主流分析工具。

---

## ✨ 主要特性

- 🔄 **格式统一** - 支持17种单细胞数据格式，统一转换为 H5AD
- 🤖 **自动检测** - 智能识别数据格式，无需手动指定
- ⚡ **高性能** - 支持稀疏矩阵、缓存机制、backed 模式
- 🔧 **灵活使用** - 命令行工具 + Python API
- 📦 **开箱即用** - 完整的测试数据和示例
- 📚 **文档完善** - 详细的中文使用文档

---

## 🎯 支持的数据格式（17种）

| 平台/工具 | 格式 | 数量 | 支持状态 |
|----------|------|------|---------|
| **10X Genomics** | MTX v2/v3, H5 v2/v3 | 4 | ✅ 完全支持 |
| **标准格式** | H5AD, Loom, Zarr | 3 | ✅ 完全支持 |
| **文本格式** | CSV, TSV, TXT, Excel | 6 | ✅ 完全支持 |
| **其他格式** | MTX单文件, UMI-tools, 自定义10X | 3 | ✅ 完全支持 |
| **MGI DNBelab** | dnbc4tools 输出（10X兼容） | 1 | ✅ 完全支持 |

> **测试覆盖**: 17/17 全部通过 ✅

详细说明请查看 [docs/data_formats.md](docs/data_formats.md)

---

## 🚀 快速开始

### 安装

#### 方式1：使用 Conda（推荐）

```bash
# 克隆仓库
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub

# 创建独立的 conda 环境
conda env create -f environment.yml

# 激活环境
conda activate scrna-datahub

# 验证安装
python src/universal_reader.py --help
```

#### 方式2：使用 pip + venv

```bash
# 克隆仓库
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub

# 创建 Python 虚拟环境
python -m venv venv

# 激活虚拟环境（Linux/macOS）
source venv/bin/activate

# 安装依赖
pip install -r requirements.txt
```

### 基础使用

#### 命令行

```bash
# 读取 10X Genomics 数据
python src/universal_reader.py \
  filtered_feature_bc_matrix/ \
  -o output.h5ad

# 读取 DNB C4 数据
python src/universal_reader.py \
  02.count/filter_matrix/ \
  -o output.h5ad

# 读取 STARsolo 输出（未压缩）
python src/universal_reader.py \
  Solo.out/Gene/filtered/ \
  -o output.h5ad \
  --no-compressed
```

#### Python API

```python
from src.universal_reader import UniversalScRNAReader

# 创建读取器
reader = UniversalScRNAReader(verbose=True)

# 自动读取任意格式
adata = reader.read_auto('filtered_feature_bc_matrix/')

# 保存为 H5AD
reader.save_h5ad(adata, 'output.h5ad')
```

---

## 📖 文档

- [快速开始](docs/quickstart.md) - 5分钟入门教程
- [安装指南](docs/installation.md) - 详细的安装说明
- [使用教程](docs/user_guide.md) - 完整的使用文档
- [API文档](docs/api_reference.md) - API 参考
- [数据格式](docs/data_formats.md) - 所有支持的格式详解
- [常见问题](docs/faq.md) - FAQ

---

## 🧪 测试

```bash
# 快速测试（5分钟，测试 6 种格式）
cd scripts/test_data_all_formats
bash QUICK_START.sh

# 完整测试（测试 17 种格式）
cd scripts/test_data_all_formats
bash run_all_format_tests_simple.sh
```

---

## 💡 使用示例

更多示例请查看 [examples/](examples/) 目录：

- `basic_usage.py` - 基础使用方法
- `batch_processing.py` - 批量处理多个样本
- `multi_sample.py` - 多样本合并和批次校正

---

## 📊 项目结构

```
scRNA-DataHub/
├── src/                    # 源代码
│   └── universal_reader.py # 核心读取器
├── docs/                   # 文档
├── examples/               # 示例代码
├── scripts/                # 测试脚本
│   └── test_data_all_formats/  # 格式测试
├── tests/                  # 单元测试
├── README_CN.md            # 中文说明
├── README.md               # 英文说明
├── requirements.txt        # Python依赖
└── environment.yml         # Conda环境配置
```

查看 [目录结构说明.md](目录结构说明.md) 了解详细结构。

---

## 🤝 贡献

欢迎贡献！请查看 [CONTRIBUTING.md](CONTRIBUTING.md)。

---

## 📜 开源协议

MIT 协议 - 查看 [LICENSE](LICENSE) 文件。

---

## 🙏 致谢

感谢以下项目：

- [Scanpy](https://github.com/scverse/scanpy) - 单细胞分析工具
- [AnnData](https://github.com/scverse/anndata) - 数据结构
- [10X Genomics](https://www.10xgenomics.com/) - 数据格式标准
- [MGI DNBelab](https://github.com/MGI-tech-bioinformatics) - DNB C4 支持

---

## 📮 联系方式

- GitHub Issues: [提交问题](https://github.com/yourusername/scRNA-DataHub/issues)
- Email: your.email@example.com

---

**Made with ❤️ for single-cell community**
