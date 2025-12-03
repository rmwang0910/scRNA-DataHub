# 全格式测试数据使用指南

## 目录结构

```
test_data_all_formats/
├── README.md                              # 本文件
├── run_all_format_tests.sh               # 全格式自动测试脚本 ⭐
├── test_data_manifest.txt                # 数据清单
│
├── 10x_pbmc68k_reduced.h5ad             # H5AD格式
├── filtered_feature_bc_matrix/          # 10X MTX v3 (压缩)
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── hg19_chr21/                          # 10X MTX v2 (未压缩)
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── filtered_feature_bc_matrix.h5        # 10X H5 v3
├── filtered_gene_bc_matrices_h5.h5      # 10X H5 v2
├── test_data.loom                       # Loom格式
├── 10x-10k-subset.zarr/                 # Zarr格式
├── test_expression.csv                  # CSV格式
├── test_expression.csv.gz               # CSV (压缩)
├── test_expression.tsv                  # TSV格式
├── test_expression.tsv.gz               # TSV (压缩)
├── krumsiek11.txt                       # TXT格式
├── test_expression.xlsx                 # Excel格式
├── test_matrix.mtx                      # MTX单文件
├── test_matrix.mtx.gz                   # MTX单文件 (压缩)
├── test_data.hdf5                       # HDF5格式
├── umi_tools_counts.tsv.gz              # UMI-tools格式
└── custom_10x_mtx/                      # 自定义10X MTX
    ├── barcodes.tsv.gz
    ├── features.tsv.gz
    └── matrix.mtx.gz
```

---

## 快速开始

### 1. 运行全格式测试

有两个版本可选：

#### 版本1: 完整版（带详细验证）

```bash
cd /storeData/ztron/wangrm/tools/scRNA-DataHub/scripts/test_data_all_formats

# 运行完整版测试脚本
bash run_all_format_tests.sh
```

**特点：**
- ✅ 显示每个文件的细胞数、基因数等详细信息
- ⚠️ 可能因为验证步骤提前退出（如果遇到问题）

#### 版本2: 简化版（更稳定）⭐ 推荐

```bash
cd /storeData/ztron/wangrm/tools/scRNA-DataHub/scripts/test_data_all_formats

# 运行简化版测试脚本（更稳定）
bash run_all_format_tests_simple.sh
```

**特点：**
- ✅ 保证运行所有18个测试
- ✅ 不会因为单个验证失败而退出
- ✅ 更容易调试
- ✅ **显示详细统计信息**（细胞数、基因数、文件大小）⭐ 新功能
- ✅ **生成测试报告**（3个日志文件）⭐ 新功能
- ✅ **表格化输出**（清晰的测试结果展示）⭐ 新功能

**脚本会提示输入绝对路径：**
```
请输入测试数据目录的绝对路径:
示例: /storeData/ztron/wangrm/tools/scRNA-DataHub/scripts/test_data_all_formats
路径: 
```

**输入你的实际路径（复制粘贴）：**
```
/storeData/ztron/wangrm/tools/scRNA-DataHub/scripts/test_data_all_formats
```

### 2. 查看测试结果

测试完成后，所有输出文件会保存在 `test_outputs/` 目录：

```bash
ls -lh test_outputs/

# 输出示例:
# 10x_mtx_v3.h5ad
# 10x_mtx_v2.h5ad
# 10x_h5_v3.h5ad
# h5ad.h5ad
# loom.h5ad
# zarr.h5ad
# csv.h5ad
# ...
```

### 3. 验证单个文件

```python
import scanpy as sc

# 读取任意输出文件
adata = sc.read_h5ad('test_outputs/10x_mtx_v3.h5ad')
print(adata)
```

---

## 测试覆盖的18种格式

| # | 格式 | 测试文件 | 说明 |
|---|------|---------|------|
| 1 | 10X MTX v3 | `filtered_feature_bc_matrix/` | 压缩版本 |
| 2 | 10X MTX v2 | `hg19_chr21/` | 未压缩版本 |
| 3 | 10X H5 v3 | `filtered_feature_bc_matrix.h5` | Cell Ranger v3+ |
| 4 | 10X H5 v2 | `filtered_gene_bc_matrices_h5.h5` | Cell Ranger v2 |
| 5 | H5AD | `10x_pbmc68k_reduced.h5ad` | Scanpy标准 |
| 6 | Loom | `test_data.loom` | 单细胞专用 |
| 7 | Zarr | `10x-10k-subset.zarr/` | 云原生格式 |
| 8 | CSV | `test_expression.csv` | 逗号分隔 |
| 9 | CSV (GZ) | `test_expression.csv.gz` | 压缩版本 |
| 10 | TSV | `test_expression.tsv` | 制表符分隔 |
| 11 | TSV (GZ) | `test_expression.tsv.gz` | 压缩版本 |
| 12 | TXT | `krumsiek11.txt` | 文本格式 |
| 13 | Excel | `test_expression.xlsx` | Excel表格 |
| 14 | MTX | `test_matrix.mtx` | 稀疏矩阵 |
| 15 | MTX (GZ) | `test_matrix.mtx.gz` | 压缩版本 |
| 16 | HDF5 | `test_data.hdf5` | 通用HDF5 |
| 17 | UMI-tools | `umi_tools_counts.tsv.gz` | UMI-tools输出 |
| 18 | Custom 10X | `custom_10x_mtx/` | 自定义10X格式 |

---

## 单独测试某个格式

### 示例1: 测试10X MTX格式

```bash
cd /storeData/ztron/wangrm/tools/scRNA-DataHub

# 默认是详细模式（verbose）
python src/universal_reader.py \
  scripts/test_data_all_formats/filtered_feature_bc_matrix \
  -o test_10x_mtx.h5ad \
  --sample-id test_sample

# 或使用静默模式
python src/universal_reader.py \
  scripts/test_data_all_formats/filtered_feature_bc_matrix \
  -o test_10x_mtx.h5ad \
  --sample-id test_sample \
  --quiet
```

### 示例2: 测试Loom格式

```bash
python src/universal_reader.py \
  scripts/test_data_all_formats/test_data.loom \
  -o test_loom.h5ad \
  --sample-id loom_sample
```

### 示例3: 测试CSV格式

```bash
python src/universal_reader.py \
  scripts/test_data_all_formats/test_expression.csv \
  -o test_csv.h5ad \
  --sample-id csv_sample \
  --delimiter ","
```

---

## Python API测试

```python
import sys
sys.path.insert(0, '/storeData/ztron/wangrm/tools/scRNA-DataHub/src')

from universal_reader import UniversalScRNAReader
import scanpy as sc

# 创建reader
reader = UniversalScRNAReader(verbose=True)

# 测试数据目录
test_dir = '/storeData/ztron/wangrm/tools/scRNA-DataHub/scripts/test_data_all_formats'

# 1. 测试10X MTX
adata1 = reader.read_auto(f'{test_dir}/filtered_feature_bc_matrix')
print(f"10X MTX: {adata1.n_obs} cells × {adata1.n_vars} genes")

# 2. 测试H5AD
adata2 = reader.read_auto(f'{test_dir}/10x_pbmc68k_reduced.h5ad')
print(f"H5AD: {adata2.n_obs} cells × {adata2.n_vars} genes")

# 3. 测试Loom
adata3 = reader.read_auto(f'{test_dir}/test_data.loom')
print(f"Loom: {adata3.n_obs} cells × {adata3.n_vars} genes")

# 4. 测试CSV
adata4 = reader.read_auto(f'{test_dir}/test_expression.csv', delimiter=',')
print(f"CSV: {adata4.n_obs} cells × {adata4.n_vars} genes")

# 5. 测试Zarr
adata5 = reader.read_auto(f'{test_dir}/10x-10k-subset.zarr')
print(f"Zarr: {adata5.n_obs} cells × {adata5.n_vars} genes")
```

---

## 测试脚本功能

`run_all_format_tests.sh` 提供：

1. ✅ **交互式路径输入** - 支持任意安装位置
2. ✅ **自动格式检测** - 测试18种数据格式
3. ✅ **环境检查** - 验证Python和依赖包
4. ✅ **详细日志** - 彩色输出，易于调试
5. ✅ **结果验证** - 自动检查输出文件
6. ✅ **统计报告** - 显示通过/失败/跳过数量
7. ✅ **错误处理** - 单个失败不影响其他测试

---

## 测试输出示例

```
========================================================================
scRNA-DataHub 全格式测试
========================================================================

请输入测试数据目录的绝对路径:
示例: /storeData/ztron/wangrm/tools/scRNA-DataHub/scripts/test_data_all_formats
路径: /storeData/ztron/wangrm/tools/scRNA-DataHub/scripts/test_data_all_formats

✓ 测试数据目录: /storeData/ztron/wangrm/tools/scRNA-DataHub/scripts/test_data_all_formats
✓ Reader脚本: /storeData/ztron/wangrm/tools/scRNA-DataHub/src/universal_reader.py
✓ 输出目录: /storeData/ztron/wangrm/tools/scRNA-DataHub/scripts/test_data_all_formats/test_outputs

========================================================================
环境检查
========================================================================

Python版本:
Python 3.10.14

关键依赖包:
  ✓ scanpy: 1.10.0
  ✓ anndata: 0.10.5
  ✓ pandas: 2.2.1
  ✓ numpy: 1.26.4
  ✓ scipy: 1.13.0
  ✓ h5py: 3.10.0
  ✓ loompy: 3.0.7
  ✓ zarr: 2.17.1

========================================================================
开始格式测试
========================================================================

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
测试格式: 10X MTX v3 (压缩)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

ℹ 输入: /path/to/filtered_feature_bc_matrix
ℹ 输出: /path/to/test_outputs/10x_mtx_v3.h5ad

✓ 10X MTX v3 (压缩) 测试通过 (文件大小: 1.2M)
  - 细胞数: 2700
  - 基因数: 32738
  - obs列: ['sample_id']
  - var列: ['gene_ids', 'feature_types']

[... 更多测试输出 ...]

========================================================================
测试结果汇总
========================================================================

总测试数: 18
通过: 18
失败: 0
跳过: 0

✓ 所有测试通过！ 🎉

输出文件列表:
total 24M
-rw-r--r-- 1 user user 1.2M Dec  3 11:00 10x_mtx_v3.h5ad
-rw-r--r-- 1 user user 890K Dec  3 11:00 10x_mtx_v2.h5ad
-rw-r--r-- 1 user user 1.5M Dec  3 11:00 10x_h5_v3.h5ad
[... 更多文件 ...]
```

---

## 故障排除

### 问题1: 找不到reader脚本

**错误：**
```
找不到reader脚本: /path/to/src/universal_reader.py
```

**解决：**
```bash
# 确保从正确的目录运行
cd /storeData/ztron/wangrm/tools/scRNA-DataHub/scripts/test_data_all_formats
bash run_all_format_tests.sh
```

### 问题2: 缺少依赖包

**错误：**
```
ModuleNotFoundError: No module named 'loompy'
```

**解决：**
```bash
# 激活conda环境
conda activate omicverse

# 或安装缺失的包
pip install loompy
```

### 问题3: 权限错误

**错误：**
```
Permission denied: run_all_format_tests.sh
```

**解决：**
```bash
chmod +x run_all_format_tests.sh
```

---

## 数据来源

所有测试数据均来自公开数据集：

- **10X数据**: scanpy内置数据集 (`sc.datasets.*`)
- **PBMC数据**: 10X Genomics公开数据
- **合成数据**: 使用`create_test_data.py`生成

---

## 自定义测试

### 添加新的测试数据

1. 将数据文件放入 `test_data_all_formats/` 目录
2. 编辑 `run_all_format_tests.sh`，添加测试：

```bash
test_format \
    "我的格式" \
    "$TEST_DATA_DIR/my_data.h5ad" \
    "$OUTPUT_DIR/my_test.h5ad" \
    "--sample-id my_sample"
```

### 创建新的测试数据

```bash
cd /storeData/ztron/wangrm/tools/scRNA-DataHub/scripts

# 运行数据生成脚本
python create_test_data.py
```

---

## 联系与支持

如有问题，请查看：
- 主README: `../../README_CN.md`
- 文档目录: `../../docs/`
- 故障排除: `../../docs/troubleshooting.md`

---

**测试数据完整度: 18/18 (100%)** ✅

