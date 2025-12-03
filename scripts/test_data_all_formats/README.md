# 全格式测试数据使用指南

本目录提供一套测试数据，其中 **自动化测试脚本当前覆盖 17 种主流单细胞/表达矩阵格式**，用于验证 `scRNA-DataHub` 的 `universal_reader` 是否能正确读取、统一转换为 `AnnData (h5ad)`。

---

## 目录结构

```text
test_data_all_formats/
├── README.md                          # 本文件
├── QUICK_START.sh                     # 一键快速测试脚本（推荐）
├── run_all_format_tests_simple.sh     # 全格式自动测试脚本（稳定版）
├── test_data_manifest.txt             # 数据清单
│
├── 10x_pbmc68k_reduced.h5ad           # H5AD格式
├── filtered_feature_bc_matrix/        # 10X MTX v3 (压缩)
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── hg19_chr21/                        # 10X MTX v2 (未压缩)
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── filtered_feature_bc_matrix.h5      # 10X H5 v3
├── filtered_gene_bc_matrices_h5.h5    # 10X H5 v2
├── test_data.loom                     # Loom格式
├── 10x-10k-subset.zarr/               # Zarr格式
├── test_expression.csv                # CSV格式
├── test_expression.csv.gz             # CSV (压缩)
├── test_expression.tsv                # TSV格式
├── test_expression.tsv.gz             # TSV (压缩)
├── krumsiek11.txt                     # TXT格式
├── test_expression.xlsx               # Excel格式
├── test_matrix.mtx                    # MTX单文件
├── test_matrix.mtx.gz                 # MTX单文件 (压缩)
├── umi_tools_counts.tsv.gz            # UMI-tools格式
└── custom_10x_mtx/                    # 自定义10X MTX
    ├── barcodes.tsv.gz
    ├── features.tsv.gz
    └── matrix.mtx.gz
```

项目根目录记为：

```text
PROJECT_ROOT = /path/to/scRNA-DataHub
```

本目录路径为：

```text
${PROJECT_ROOT}/scripts/test_data_all_formats
```

---

## 一键快速体验（推荐）

如果你只是想快速确认环境和测试数据是否正常，直接执行：

```bash
cd ${PROJECT_ROOT}/scripts/test_data_all_formats
bash QUICK_START.sh
```

**特点：**
- ✅ 自动检测并提示 Python/依赖是否就绪  
- ✅ 自动调用 `run_all_format_tests_simple.sh` 运行全量测试  
- ✅ 在当前目录创建 `test_outputs/` 保存所有转换得到的 `h5ad`  

不需要手动输入路径，对新用户最友好。

---

## 全格式稳定测试脚本

如果你需要 **完整、可重复的格式兼容性测试**，建议直接使用稳定版脚本：

```bash
cd ${PROJECT_ROOT}/scripts/test_data_all_formats
bash run_all_format_tests_simple.sh
```

执行过程中脚本会 **提示输入测试数据目录的绝对路径**，例如：

```text
请输入测试数据目录的绝对路径:
示例: /abs/path/to/scRNA-DataHub/scripts/test_data_all_formats
路径:
```

你只需复制本目录的绝对路径粘贴进去即可。

**脚本特性：**
- ✅ 覆盖所有 17 种格式测试  
- ✅ 单个格式失败不会中断整个流程，方便排查  
- ✅ 输出每种格式的细胞数、基因数、文件大小等统计信息  
- ✅ 在 `test_outputs/` 目录下生成结果 `h5ad` 文件  
- ✅ 生成日志/表格化结果（便于记录和比对）  

---

## 查看测试输出

所有转换后的 `h5ad` 文件统一保存在本目录下的 `test_outputs/` 中：

```bash
cd ${PROJECT_ROOT}/scripts/test_data_all_formats
ls -lh test_outputs/
```

你可以任选一个文件用 `scanpy` 打开：

```python
import scanpy as sc

adata = sc.read_h5ad("test_outputs/10x_mtx_v3.h5ad")
print(adata)
```

---

## 测试覆盖的 17 种格式（自动测试）

| #  | 格式        | 测试文件                           | 说明             |
|----|-------------|------------------------------------|------------------|
| 1  | 10X MTX v3  | `filtered_feature_bc_matrix/`      | 压缩版本         |
| 2  | 10X MTX v2  | `hg19_chr21/`                      | 未压缩版本       |
| 3  | 10X H5 v3   | `filtered_feature_bc_matrix.h5`    | Cell Ranger v3+  |
| 4  | 10X H5 v2   | `filtered_gene_bc_matrices_h5.h5`  | Cell Ranger v2   |
| 5  | H5AD        | `10x_pbmc68k_reduced.h5ad`         | Scanpy 标准      |
| 6  | Loom        | `test_data.loom`                   | Loom 单细胞格式  |
| 7  | Zarr        | `10x-10k-subset.zarr/`             | 云原生格式       |
| 8  | CSV         | `test_expression.csv`              | 逗号分隔         |
| 9  | CSV (GZ)    | `test_expression.csv.gz`           | 压缩             |
| 10 | TSV         | `test_expression.tsv`              | 制表符分隔       |
| 11 | TSV (GZ)    | `test_expression.tsv.gz`           | 压缩             |
| 12 | TXT         | `krumsiek11.txt`                   | 普通文本         |
| 13 | Excel       | `test_expression.xlsx`             | Excel 表格       |
| 14 | MTX         | `test_matrix.mtx`                  | 稀疏矩阵         |
| 15 | MTX (GZ)    | `test_matrix.mtx.gz`               | 压缩             |
| 16 | UMI-tools   | `umi_tools_counts.tsv.gz`          | UMI-tools 输出   |
| 17 | Custom 10X  | `custom_10x_mtx/`                  | 自定义 10X 格式  |

---

## 使用命令行单独测试某种格式

以下示例均假设当前目录为项目根目录：

```bash
cd ${PROJECT_ROOT}
```

### 示例 1：测试 10X MTX v3 目录

```bash
python src/universal_reader.py \
  scripts/test_data_all_formats/filtered_feature_bc_matrix \
  -o test_10x_mtx.h5ad \
  --sample-id test_sample

# 静默模式（减少日志）
python src/universal_reader.py \
  scripts/test_data_all_formats/filtered_feature_bc_matrix \
  -o test_10x_mtx.h5ad \
  --sample-id test_sample \
  --quiet
```

### 示例 2：测试 Loom 文件

```bash
python src/universal_reader.py \
  scripts/test_data_all_formats/test_data.loom \
  -o test_loom.h5ad \
  --sample-id loom_sample
```

### 示例 3：测试 CSV 文件

```bash
python src/universal_reader.py \
  scripts/test_data_all_formats/test_expression.csv \
  -o test_csv.h5ad \
  --sample-id csv_sample \
  --delimiter ","
```

---

## Python API 示例

```python
from pathlib import Path
import scanpy as sc
from universal_reader import UniversalScRNAReader

# 项目根目录，根据实际情况修改
PROJECT_ROOT = Path("/path/to/scRNA-DataHub")
test_dir = PROJECT_ROOT / "scripts" / "test_data_all_formats"

reader = UniversalScRNAReader(verbose=True)

# 1. 10X MTX
adata1 = reader.read_auto(test_dir / "filtered_feature_bc_matrix")
print(f"10X MTX: {adata1.n_obs} cells × {adata1.n_vars} genes")

# 2. H5AD
adata2 = reader.read_auto(test_dir / "10x_pbmc68k_reduced.h5ad")
print(f"H5AD: {adata2.n_obs} cells × {adata2.n_vars} genes")

# 3. Loom
adata3 = reader.read_auto(test_dir / "test_data.loom")
print(f"Loom: {adata3.n_obs} cells × {adata3.n_vars} genes")

# 4. CSV
adata4 = reader.read_auto(test_dir / "test_expression.csv", delimiter=",")
print(f"CSV: {adata4.n_obs} cells × {adata4.n_vars} genes")

# 5. Zarr
adata5 = reader.read_auto(test_dir / "10x-10k-subset.zarr")
print(f"Zarr: {adata5.n_obs} cells × {adata5.n_vars} genes")
```

---

## 常见问题（FAQ & 故障排除）

### 问题 1：找不到 `universal_reader.py`

**症状：**

```text
找不到 reader 脚本: src/universal_reader.py
```

**排查：**

- 确认当前目录在项目根目录：`cd ${PROJECT_ROOT}`  
- 确认 `src/universal_reader.py` 文件存在  
- 若通过 `pip` 安装为包使用，请参考项目根目录的 `README_CN.md` 或 `docs/installation.md`  

---

### 问题 2：缺少依赖包（例如 `loompy`、`zarr` 等）

**错误示例：**

```text
ModuleNotFoundError: No module named 'loompy'
```

**解决方法：**

```bash
# 推荐使用项目提供的环境
cd ${PROJECT_ROOT}
conda env create -f environment.yml    # 首次创建
conda activate scrna-datahub           # 环境名称以实际配置为准

# 如仍缺少个别包，可单独安装
pip install loompy zarr
```

---

### 问题 3：脚本没有执行权限

**错误示例：**

```text
Permission denied: run_all_format_tests_simple.sh
```

**解决方法：**

```bash
cd ${PROJECT_ROOT}/scripts/test_data_all_formats
chmod +x QUICK_START.sh run_all_format_tests_simple.sh
```

---

## 数据来源说明

- **10X 数据**：部分来自 `scanpy` 内置公开数据集（`sc.datasets.*`）  
- **PBMC 数据**：来源于 10x Genomics 官方公开数据  
- **合成/裁剪数据**：由脚本 `scripts/create_test_data.py` 生成或下采样获得，仅用于示例和兼容性测试  

数据仅用于演示与测试，不推荐直接用于真实科研分析。

---

## 扩展与自定义测试

### 添加新的测试数据格式

1. 将你的数据文件放入 `test_data_all_formats/` 目录，或建立新的子目录  
2. 修改 `run_all_format_tests_simple.sh`，增加一条测试条目（伪代码示例）：

```bash
test_format \
  "我的新格式" \
  "$TEST_DATA_DIR/my_data.h5ad" \
  "$OUTPUT_DIR/my_data_converted.h5ad" \
  "--sample-id my_sample"
```

### 使用生成脚本创建测试数据

```bash
cd ${PROJECT_ROOT}/scripts
python create_test_data.py
```

脚本会在合适的位置生成/更新测试数据（详见脚本内说明）。

---

## 更多文档

- 项目中文总览：`../../README_CN.md`  
- 详细用户指南：`../../docs/user_guide.md`  
- 安装与环境管理：`../../docs/installation.md`、`../../docs/environment_management.md`  
- 故障排除汇总：`../../docs/troubleshooting.md`  

如需补充或发现文档错误，可在仓库中提交 Issue 或 PR。  

**当前自动测试覆盖度：17/17（100%）✅**

