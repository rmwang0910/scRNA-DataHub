# scRNA-DataHub

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-linux%20%7C%20macOS-lightgrey)](https://github.com)

**scRNA-DataHub** 是一个通用的单细胞RNA测序数据处理工具，支持12+种数据格式的读取、转换和标准化，统一输出为H5AD格式，可直接用于Scanpy、Seurat等主流分析工具。

## ✨ 主要特性

- 🔄 **格式统一**：支持12+种单细胞数据格式，统一转换为H5AD
- 🤖 **自动检测**：智能识别数据格式，无需手动指定
- ⚡ **高性能**：支持稀疏矩阵、缓存机制、backed模式
- 🔧 **灵活使用**：命令行工具 + Python API
- 📦 **开箱即用**：完整的测试数据和示例
- 📚 **文档完善**：详细的使用文档和教程

## 🎯 支持的数据格式

| 格式类型 | 格式名称 | 来源 | 自动检测 |
|---------|---------|------|---------|
| **10X Genomics** | MTX格式 (v2/v3) | Cell Ranger | ✅ |
| **10X Genomics** | H5格式 | Cell Ranger | ✅ |
| **AnnData** | H5AD | Scanpy/Seurat | ✅ |
| **Loom** | Loom | Velocyto/scVelo | ✅ |
| **Zarr** | Zarr | 云存储 | ✅ |
| **文本格式** | CSV/TSV/TXT | 通用 | ✅ |
| **Excel** | XLSX/XLS | 手动整理 | ✅ |
| **稀疏矩阵** | MTX | Matrix Market | ✅ |
| **HDF5** | H5/HDF5 | 通用 | ✅ |
| **GEO** | SOFT.GZ | NCBI GEO | ✅ |
| **UMI-tools** | TSV.GZ | UMI-tools | ✅ |
| **DNB C4** | MTX (10X兼容) | dnbc4tools | ✅ |

## 🚀 快速开始

### Installation

#### Method 1: Using Conda (Recommended)

```bash
# Clone repository
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub

# Create isolated conda environment
conda env create -f environment.yml

# Activate environment
conda activate scrna-datahub

# Verify installation
python src/universal_reader.py --help
```

#### Method 2: Using pip + venv

```bash
# Clone repository
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub

# Create Python virtual environment
python -m venv venv

# Activate virtual environment
source venv/bin/activate  # Linux/macOS
# or
venv\Scripts\activate     # Windows

# Install dependencies
pip install -r requirements.txt
```

#### Method 3: Direct Installation (Not Recommended)

```bash
# Install without isolation (may conflict with other packages)
pip install -r requirements.txt
```

### 基础使用

```bash
# 命令行使用
python src/universal_reader.py input_data/ -o output.h5ad

# 示例：读取10X Genomics数据
python src/universal_reader.py \
  filtered_feature_bc_matrix/ \
  -o sample1.h5ad \
  --sample-id sample1
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

## 📖 详细文档

- [安装指南](docs/installation.md)
- [快速开始](docs/quickstart.md)
- [完整使用教程](docs/user_guide.md)
- [API文档](docs/api_reference.md)
- [数据格式详解](docs/data_formats.md)
- [常见问题](docs/faq.md)

## 💡 使用示例

### 示例1：读取10X Genomics数据

```bash
# Cell Ranger v3+ 输出
python src/universal_reader.py \
  sample1/outs/filtered_feature_bc_matrix/ \
  -o sample1.h5ad
```

### 示例2：读取DNB C4数据

```bash
# dnbc4tools输出
python src/universal_reader.py \
  CNS1063416_brain/02.count/filter_matrix/ \
  -o CNS1063416_brain.h5ad
```

### 示例3：读取STARsolo输出

```bash
# STARsolo输出（未压缩）
python src/universal_reader.py \
  Solo.out/Gene/filtered/ \
  -o sample1.h5ad \
  --no-compressed
```

### 示例4：批量处理

```python
from src.universal_reader import UniversalScRNAReader

reader = UniversalScRNAReader()

samples = ['sample1', 'sample2', 'sample3']
for sample in samples:
    adata = reader.process_single_sample(
        input_path=f'{sample}/filtered_feature_bc_matrix/',
        output_path=f'{sample}.h5ad',
        sample_id=sample
    )
```

更多示例请查看 [examples/](examples/) 目录。

## 🧪 测试

### 快速测试

```bash
# 使用scanpy内置数据快速测试
bash scripts/quick_test.sh
```

### 完整测试

```bash
# 创建所有格式的测试数据
python scripts/create_test_data.py

# 运行完整测试套件
cd tests
bash run_all_tests.sh
```

## 📊 项目结构

```
scRNA-DataHub/
├── README.md                   # 项目说明
├── LICENSE                     # 开源协议
├── requirements.txt            # Python依赖
├── environment.yml             # Conda环境
├── setup.py                    # 安装脚本
├── pyproject.toml             # 项目配置
├── .gitignore                 # Git忽略文件
├── src/                       # 源代码
│   ├── __init__.py
│   ├── universal_reader.py    # 核心读取器
│   └── utils.py               # 工具函数
├── tests/                     # 测试
│   ├── test_reader.py
│   └── test_formats.py
├── scripts/                   # 脚本
│   ├── quick_test.sh
│   └── create_test_data.py
├── examples/                  # 示例
│   ├── basic_usage.py
│   ├── batch_processing.py
│   └── multi_sample.py
├── docs/                      # 文档
│   ├── installation.md
│   ├── quickstart.md
│   ├── user_guide.md
│   ├── api_reference.md
│   ├── data_formats.md
│   └── faq.md
└── .github/                   # GitHub配置
    └── workflows/
        └── tests.yml          # CI/CD配置
```

## 🤝 贡献

欢迎贡献！请查看 [CONTRIBUTING.md](CONTRIBUTING.md) 了解如何参与项目。

### 贡献方式

1. Fork 本仓库
2. 创建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 开启Pull Request

## 📜 开源协议

本项目采用 MIT 协议 - 查看 [LICENSE](LICENSE) 文件了解详情。

## 👥 作者

- **Wang Ruiming** - *Initial work*

## 🙏 致谢

- [Scanpy](https://github.com/scverse/scanpy) - 单细胞分析工具
- [AnnData](https://github.com/scverse/anndata) - 数据结构
- [10X Genomics](https://www.10xgenomics.com/) - 数据格式标准

## 📮 联系方式

- Issues: [GitHub Issues](https://github.com/yourusername/scRNA-DataHub/issues)
- Email: your.email@example.com

## 🔗 相关项目

- [Scanpy](https://github.com/scverse/scanpy) - Python单细胞分析
- [Seurat](https://github.com/satijalab/seurat) - R单细胞分析
- [OmicVerse](https://github.com/Starlitnightly/omicverse) - 多组学分析
- [dnbc4tools](https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software) - DNB C4分析

## 📈 更新日志

查看 [CHANGELOG.md](CHANGELOG.md) 了解版本更新历史。

## ⭐ Star History

如果这个项目对您有帮助，请给我们一个Star！

---

**Made with ❤️ for single-cell community**

