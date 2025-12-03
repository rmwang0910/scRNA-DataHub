# scRNA-DataHub

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-linux%20%7C%20macOS-lightgrey)](https://github.com)

English | [简体中文](README_CN.md)

**scRNA-DataHub** is a universal single-cell RNA-seq data processing tool that supports **17 data formats** for reading, converting, and standardizing, with unified output to H5AD format for direct use with mainstream analysis tools like Scanpy and Seurat.

---

## ✨ Key Features

- 🔄 **Format Unification** - Supports 17 single-cell data formats, unified conversion to H5AD
- 🤖 **Auto Detection** - Intelligent format recognition, no manual specification needed
- ⚡ **High Performance** - Supports sparse matrices, caching, backed mode
- 🔧 **Flexible Usage** - Command-line tool + Python API
- 📦 **Ready to Use** - Complete test data and examples
- 📚 **Well Documented** - Detailed documentation and tutorials

---

## 🎯 Supported Data Formats (17 types)

| Platform/Tool | Format | Count | Support Status |
|--------------|--------|-------|----------------|
| **10X Genomics** | MTX v2/v3, H5 v2/v3 | 4 | ✅ Full support |
| **Standard Formats** | H5AD, Loom, Zarr | 3 | ✅ Full support |
| **Text Formats** | CSV, TSV, TXT, Excel | 6 | ✅ Full support |
| **Other Formats** | MTX single file, UMI-tools, Custom 10X MTX | 3 | ✅ Full support |
| **Special Variants** | DNB C4 (subset of Custom 10X MTX) | - | ✅ Auto-detected |

> **Test Coverage**: 17/17 all passed ✅  
> **Note**: DNB C4 format is a variant of "Custom 10X MTX", automatically detected and handled

See [docs/data_formats.md](docs/data_formats.md) for details.

---

## 🚀 Quick Start

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

# Activate virtual environment (Linux/macOS)
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### Basic Usage

#### Command Line

```bash
# Read 10X Genomics data
python src/universal_reader.py \
  filtered_feature_bc_matrix/ \
  -o output.h5ad

# Read DNB C4 data
python src/universal_reader.py \
  02.count/filter_matrix/ \
  -o output.h5ad

# Read STARsolo output (uncompressed)
python src/universal_reader.py \
  Solo.out/Gene/filtered/ \
  -o output.h5ad \
  --no-compressed
```

#### Python API

```python
from src.universal_reader import UniversalScRNAReader

# Create reader
reader = UniversalScRNAReader(verbose=True)

# Auto-read any format
adata = reader.read_auto('filtered_feature_bc_matrix/')

# Save as H5AD
reader.save_h5ad(adata, 'output.h5ad')
```

---

## 📖 Documentation

- [Quick Start](docs/quickstart.md) - 5-minute tutorial
- [Installation Guide](docs/installation.md) - Detailed installation instructions
- [User Guide](docs/user_guide.md) - Complete usage documentation
- [API Reference](docs/api_reference.md) - API documentation
- [Data Formats](docs/data_formats.md) - All supported formats
- [FAQ](docs/faq.md) - Frequently asked questions

---

## 🧪 Testing

```bash
# Quick test (5 minutes, tests 6 formats)
cd scripts/test_data_all_formats
bash QUICK_START.sh

# Complete test (tests 17 formats)
cd scripts/test_data_all_formats
bash run_all_format_tests_simple.sh
```

---

## 💡 Examples

See [examples/](examples/) directory for more examples:

- `basic_usage.py` - Basic usage methods
- `batch_processing.py` - Batch processing multiple samples
- `multi_sample.py` - Multi-sample integration and batch correction

---

## 📊 Project Structure

```
scRNA-DataHub/
├── src/                    # Source code
│   └── universal_reader.py # Core reader
├── docs/                   # Documentation
├── examples/               # Example code
├── scripts/                # Test scripts
│   └── test_data_all_formats/  # Format tests
├── tests/                  # Unit tests
├── README.md               # English README
├── README_CN.md            # Chinese README
├── requirements.txt        # Python dependencies
└── environment.yml         # Conda environment config
```

See [目录结构说明.md](目录结构说明.md) for detailed structure.

---

## 📜 License

MIT License - see [LICENSE](LICENSE) file.

---

## 🙏 Acknowledgments

Thanks to the following projects:

- [Scanpy](https://github.com/scverse/scanpy) - Single-cell analysis tool
- [AnnData](https://github.com/scverse/anndata) - Data structure
- [10X Genomics](https://www.10xgenomics.com/) - Data format standards
- [MGI DNBelab](https://github.com/MGI-tech-bioinformatics) - DNB C4 support

---

## 📮 Contact

- GitHub Issues: [Submit issues](https://github.com/yourusername/scRNA-DataHub/issues)
- Email: rmwang0910@gmail.com

---

**Made with ❤️ for single-cell community**
