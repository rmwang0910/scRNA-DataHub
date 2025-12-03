# scRNA-DataHub

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-linux%20%7C%20macOS-lightgrey)](https://github.com)

[English](README.md) | 简体中文

**scRNA-DataHub** 是一个通用的单细胞RNA测序数据处理工具，支持12+种数据格式的读取、转换和标准化，统一输出为H5AD格式，可直接用于Scanpy、Seurat等主流分析工具。

## ✨ 主要特性

- 🔄 **格式统一**：支持12+种单细胞数据格式，统一转换为H5AD
- 🤖 **自动检测**：智能识别数据格式，无需手动指定
- ⚡ **高性能**：支持稀疏矩阵、缓存机制、backed模式
- 🔧 **灵活使用**：命令行工具 + Python API
- 📦 **开箱即用**：完整的测试数据和示例
- 📚 **文档完善**：详细的中文使用文档和教程

## 🎯 支持的数据格式

| 平台/工具 | 格式 | 支持状态 |
|----------|------|---------|
| **10X Genomics** | Cell Ranger MTX (v2/v3) | ✅ 完全支持 |
| **10X Genomics** | Cell Ranger H5 | ✅ 完全支持 |
| **MGI DNBelab** | dnbc4tools输出 (10X兼容) | ✅ 完全支持 |
| **STARsolo** | MTX输出 | ✅ 完全支持 |
| **Scanpy** | H5AD | ✅ 原生支持 |
| **Seurat** | RDS/H5Seurat | ⚠️ 需要R转换 |
| **Velocyto** | Loom | ✅ 完全支持 |
| **通用** | CSV/TSV/TXT/Excel | ✅ 完全支持 |
| **云存储** | Zarr | ✅ 完全支持 |
| **GEO数据库** | SOFT.GZ | ✅ 完全支持 |
| **其他** | MTX/UMI-tools | ✅ 完全支持 |

## 🚀 快速开始

### 安装

#### 方式1：使用Conda（推荐）

```bash
# 克隆仓库
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub

# 创建独立的conda环境
conda env create -f environment.yml

# 激活环境
conda activate scrna-datahub

# 验证安装
python src/universal_reader.py --help
```

#### 方式2：使用pip + venv

```bash
# 克隆仓库
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub

# 创建Python虚拟环境
python -m venv venv

# 激活虚拟环境
source venv/bin/activate  # Linux/macOS
# 或
venv\Scripts\activate     # Windows

# 安装依赖
pip install -r requirements.txt
```

#### 方式3：直接安装（不推荐）

```bash
# 不创建隔离环境（可能与其他包冲突）
pip install -r requirements.txt
```

### 基础使用

```bash
# 读取10X Genomics数据
python src/universal_reader.py \
  filtered_feature_bc_matrix/ \
  -o output.h5ad

# 读取DNB C4数据
python src/universal_reader.py \
  02.count/filter_matrix/ \
  -o output.h5ad

# 读取STARsolo输出
python src/universal_reader.py \
  Solo.out/Gene/filtered/ \
  -o output.h5ad \
  --no-compressed
```

### Python API

```python
from src.universal_reader import UniversalScRNAReader

# 创建读取器
reader = UniversalScRNAReader(verbose=True)

# 自动读取任意格式
adata = reader.read_auto('filtered_feature_bc_matrix/')

# 保存为H5AD
reader.save_h5ad(adata, 'output.h5ad')
```

## 📖 文档

- [安装指南](docs/installation.md) - 详细的安装说明
- [快速开始](docs/quickstart.md) - 5分钟入门教程
- [使用教程](docs/user_guide.md) - 完整的使用文档
- [API文档](docs/api_reference.md) - API参考
- [数据格式](docs/data_formats.md) - 所有支持的格式详解
- [常见问题](docs/faq.md) - FAQ

## 🧪 测试

```bash
# 快速测试（使用scanpy内置数据）
bash scripts/quick_test.sh

# 完整测试
python scripts/create_test_data.py
cd test_data_all_formats
bash run_all_format_tests.sh
```

## 💡 使用示例

### 示例1：Cell Ranger输出

```bash
python src/universal_reader.py \
  sample1/outs/filtered_feature_bc_matrix/ \
  -o sample1.h5ad \
  --sample-id sample1
```

### 示例2：dnbc4tools输出

```bash
python src/universal_reader.py \
  CNS1063416_brain/02.count/filter_matrix/ \
  -o CNS1063416_brain.h5ad \
  --sample-id CNS1063416_brain
```

### 示例3：批量处理

```python
from src.universal_reader import UniversalScRNAReader

reader = UniversalScRNAReader()

for sample in ['sample1', 'sample2', 'sample3']:
    reader.process_single_sample(
        input_path=f'{sample}/filtered_feature_bc_matrix/',
        output_path=f'{sample}.h5ad',
        sample_id=sample
    )
```

更多示例请查看 [examples/](examples/) 目录。

## 📊 项目结构

```
scRNA-DataHub/
├── src/                    # 源代码
├── tests/                  # 测试
├── scripts/                # 脚本
├── examples/               # 示例
├── docs/                   # 文档
├── README.md               # 项目说明
└── requirements.txt        # 依赖
```

查看 [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) 了解详细结构。

## 🤝 贡献

欢迎贡献！请查看 [CONTRIBUTING.md](CONTRIBUTING.md)。

## 📜 开源协议

MIT协议 - 查看 [LICENSE](LICENSE) 文件。

## 🙏 致谢

感谢以下项目：

- [Scanpy](https://github.com/scverse/scanpy)
- [AnnData](https://github.com/scverse/anndata)
- [10X Genomics](https://www.10xgenomics.com/)
- [MGI DNBelab](https://github.com/MGI-tech-bioinformatics)

## 📮 联系方式

- GitHub Issues: [提交问题](https://github.com/yourusername/scRNA-DataHub/issues)
- Email: your.email@example.com

---

**Made with ❤️ for single-cell community**

