# 单细胞数据格式测试数据获取方式汇总

## 📋 12种格式完整清单

根据您的需求，以下是除了已标注✅的4种格式外，其余8种格式的数据获取方式：

---

## ✅ 已有数据（无需额外操作）

| 序号 | 格式 | 状态 | 来源 |
|------|------|------|------|
| 1 | 10X MTX | ✅ | Scanpy本地 |
| 2 | 10X H5 | ✅ | Scanpy本地 |
| 3 | H5AD | ✅ | Scanpy本地 |
| 4 | Loom | ✅ | 需创建 |
| 6 | CSV | ✅ | 需创建 |

---

## 📦 需要获取的格式（7种）

### 5. Zarr格式 - 云原生大数据格式

**✅ Scanpy已有测试数据：**

**路径：**
```
/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x-10k-subset.zarr/
```

**测试命令：**
```bash
python universal_scrna_reader.py \
  /Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x-10k-subset.zarr/ \
  -o test_zarr.h5ad
```

**获取方式：** ✅ 已存在，直接使用

---

### 7. TSV/TXT格式 - 制表符分隔表格

**✅ Scanpy已有测试数据：**

**路径：**
```
/Users/warm/华大智造/TCGA/gdc/scanpy/src/scanpy/datasets/krumsiek11.txt
/Users/warm/华大智造/TCGA/gdc/scanpy/src/scanpy/datasets/toggleswitch.txt
```

**测试命令：**
```bash
python universal_scrna_reader.py \
  /Users/warm/华大智造/TCGA/gdc/scanpy/src/scanpy/datasets/krumsiek11.txt \
  -o test_txt.h5ad
```

**获取方式：** ✅ 已存在，直接使用

**额外创建自定义TSV：**
```bash
# 运行创建脚本会生成更多TSV格式
python create_test_data.py
# 输出: test_data_all_formats/test_expression.tsv
```

---

### 8. Excel格式 - Excel表格

**📦 需要创建：**

**方式1：运行create_test_data.py（推荐）**
```bash
# 安装依赖
pip install openpyxl

# 创建测试数据
python create_test_data.py

# 生成路径
test_data_all_formats/test_expression.xlsx
```

**方式2：手动创建**
```python
import scanpy as sc
import pandas as pd

# 读取数据
adata = sc.datasets.pbmc68k_reduced()
adata_small = adata[:50, :100].copy()  # Excel有行数限制

# 转换为DataFrame
df = pd.DataFrame(
    adata_small.X.toarray(),
    index=adata_small.obs_names,
    columns=adata_small.var_names
)

# 保存为Excel
df.to_excel('test_data.xlsx', sheet_name='expression_matrix')
```

**方式3：下载公开Excel数据**

Scanpy内置Moignard15数据集（Excel格式）：
```python
import scanpy as sc

# 会下载Excel格式数据
adata = sc.datasets.moignard15()

# 数据缓存在：~/.cache/scanpy-data/moignard15/nbt.3154-S3.xlsx
```

**测试命令：**
```bash
python universal_scrna_reader.py \
  test_data_all_formats/test_expression.xlsx \
  -o test_excel.h5ad \
  --sheet expression_matrix
```

---

### 9. MTX格式（单文件）

**📦 需要创建：**

**方式1：运行create_test_data.py（推荐）**
```bash
python create_test_data.py
# 生成: test_data_all_formats/test_matrix.mtx
#       test_data_all_formats/test_matrix.mtx.gz
```

**方式2：手动创建**
```python
import scanpy as sc
import scipy.io as sio

adata = sc.datasets.pbmc68k_reduced()
sio.mmwrite('matrix.mtx', adata.X.T)  # 转置为基因×细胞
```

**测试命令：**
```bash
python universal_scrna_reader.py \
  test_data_all_formats/test_matrix.mtx \
  -o test_mtx.h5ad
```

---

### 10. HDF5格式 - 通用HDF5格式

**📦 需要创建：**

**方式1：运行create_test_data.py（推荐）**
```bash
python create_test_data.py
# 生成: test_data_all_formats/test_data.hdf5
```

**方式2：手动创建**
```python
import scanpy as sc
import h5py

adata = sc.datasets.pbmc68k_reduced()

with h5py.File('test_data.h5', 'w') as f:
    data_group = f.create_group('data')
    data_group.create_dataset('expression', data=adata.X.toarray())
    data_group.create_dataset('gene_names', data=adata.var_names.values.astype('S'))
    data_group.create_dataset('cell_names', data=adata.obs_names.values.astype('S'))
```

**读取方式：**
```python
# HDF5需要指定key
adata = sc.read_hdf('test_data.h5', key='data')
```

---

### 11. SOFT.GZ格式 - GEO数据库格式

**🌐 需要下载（需要网络）：**

**方式1：使用Scanpy自动下载（推荐）**
```python
import scanpy as sc

# 下载GEO数据（SOFT格式）
adata = sc.datasets.burczynski06()

# 数据会自动缓存到：
# ~/.cache/scanpy-data/burczynski06/GDS1615_full.soft.gz
```

**方式2：手动从GEO下载**
```bash
# 下载GDS1615数据集
wget ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS1nnn/GDS1615/soft/GDS1615_full.soft.gz

# 或使用curl
curl -O ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS1nnn/GDS1615/soft/GDS1615_full.soft.gz
```

**推荐的小规模GEO数据集：**
- **GDS1615**：Burczynski06（127 samples × 22283 genes）~20MB
- **GDS3715**：PBMC数据
- **GSE29617**：单细胞qPCR数据

**GEO数据搜索：**
1. 访问：https://www.ncbi.nlm.nih.gov/geo/
2. 搜索：`single cell RNA`
3. 选择小规模数据集（< 1000样本）
4. 下载 SOFT 格式

**测试命令：**
```bash
# 使用缓存的数据
python universal_scrna_reader.py \
  ~/.cache/scanpy-data/burczynski06/GDS1615_full.soft.gz \
  -o test_soft_gz.h5ad

# 或使用下载的数据
python universal_scrna_reader.py \
  GDS1615_full.soft.gz \
  -o test_soft_gz.h5ad
```

---

### 12. UMI-tools格式

**📦 需要创建：**

**方式1：运行create_test_data.py（推荐）**
```bash
python create_test_data.py
# 生成: test_data_all_formats/umi_tools_counts.tsv.gz
```

**方式2：手动创建**
```python
import scanpy as sc
import pandas as pd

# 读取数据
adata = sc.datasets.pbmc68k_reduced()

# 创建UMI-tools格式（基因×细胞）
df = pd.DataFrame(
    adata.X.T.toarray() if scipy.sparse.issparse(adata.X) else adata.X.T,
    index=adata.var_names,
    columns=adata.obs_names
)

# 保存为UMI-tools格式
df.to_csv('umi_tools_counts.tsv.gz', sep='\t', compression='gzip')
```

**UMI-tools真实输出示例：**

如果有使用UMI-tools的真实数据，输出格式通常是：
```
gene    cell1    cell2    cell3    ...
GAPDH   10       15       12       ...
ACTB    20       18       19       ...
```

**测试命令：**
```bash
python universal_scrna_reader.py \
  test_data_all_formats/umi_tools_counts.tsv.gz \
  -o test_umi_tools.h5ad \
  --transpose
```

**获取真实UMI-tools数据：**
- 使用UMI-tools处理FASTQ数据
- 或从已发表的研究中下载（通常在GEO/SRA的补充数据中）

---

## 🎯 完整获取方案

### 方案A：最简方案（使用Scanpy现有数据）

**只测试Scanpy本地已有的7种格式：**

```bash
bash quick_test.sh
```

**包含格式：**
1. ✅ 10X MTX v3
2. ✅ 10X MTX v2
3. ✅ 10X H5
4. ✅ H5AD
5. ✅ TXT
6. ✅ Zarr

---

### 方案B：推荐方案（自动创建 + Scanpy现有）

**使用create_test_data.py自动创建所有格式：**

```bash
# 一键创建
python create_test_data.py

# 运行完整测试
cd test_data_all_formats
bash run_all_format_tests.sh
```

**包含格式（14种）：**
- ✅ Scanpy现有：7种
- 📦 自动创建：7种（Loom, CSV, TSV, Excel, MTX, HDF5, UMI-tools）

**注意：** SOFT.GZ需要网络连接下载

---

### 方案C：完整方案（包含下载）

```bash
# 1. 创建本地数据
python create_test_data.py

# 2. 下载SOFT.GZ数据
python -c "import scanpy as sc; sc.datasets.burczynski06()"

# 3. 运行所有测试
cd test_data_all_formats
bash run_all_format_tests.sh

# 4. 测试SOFT.GZ
cd ..
python universal_scrna_reader.py \
  ~/.cache/scanpy-data/burczynski06/GDS1615_full.soft.gz \
  -o test_soft_gz.h5ad
```

**包含格式：** 所有14种格式

---

## 📊 数据获取总结表

| 格式 | Scanpy本地 | 需创建 | 需下载 | 获取命令 |
|------|-----------|--------|--------|---------|
| 1. 10X MTX v3 | ✅ | ❌ | ❌ | 直接使用 |
| 2. 10X MTX v2 | ✅ | ❌ | ❌ | 直接使用 |
| 3. 10X H5 | ✅ | ❌ | ❌ | 直接使用 |
| 4. H5AD | ✅ | ❌ | ❌ | 直接使用 |
| 5. Loom | ❌ | ✅ | ❌ | `create_test_data.py` |
| 6. Zarr | ✅ | ❌ | ❌ | 直接使用 |
| 7. CSV | ❌ | ✅ | ❌ | `create_test_data.py` |
| 8. TSV | ✅ + 📦 | ✅ | ❌ | Scanpy有 + 可创建更多 |
| 9. Excel | ❌ | ✅ | 🌐 | `create_test_data.py` 或 `sc.datasets.moignard15()` |
| 10. MTX | ❌ | ✅ | ❌ | `create_test_data.py` |
| 11. HDF5 | ❌ | ✅ | ❌ | `create_test_data.py` |
| 12. SOFT.GZ | ❌ | ❌ | 🌐 | `sc.datasets.burczynski06()` |
| 13. UMI-tools | ❌ | ✅ | ❌ | `create_test_data.py` |

**图例：**
- ✅ Scanpy本地已有
- 📦 需要创建（运行脚本）
- 🌐 需要下载（需要网络）

---

## 🚀 快速获取方案

### 一键获取所有数据

```bash
cd /Users/warm/华大智造/TCGA/gdc

# 运行创建脚本（会自动使用scanpy现有数据 + 创建缺失格式）
python create_test_data.py
```

**输出摘要：**
```
✅ Scanpy现有测试数据:
  1. 10X MTX v3 (压缩):  scanpy/tests/_data/10x_data/3.0.0/filtered_feature_bc_matrix/
  2. 10X MTX v2 (未压缩): scanpy/tests/_data/10x_data/1.2.0/filtered_gene_bc_matrices/hg19_chr21/
  3. 10X H5 v3:          scanpy/tests/_data/10x_data/3.0.0/filtered_feature_bc_matrix.h5
  4. 10X H5 v2:          scanpy/tests/_data/10x_data/1.2.0/filtered_gene_bc_matrices_h5.h5
  5. H5AD:               scanpy/src/scanpy/datasets/10x_pbmc68k_reduced.h5ad
  6. TXT:                scanpy/src/scanpy/datasets/krumsiek11.txt
  7. Zarr:               scanpy/tests/_data/10x-10k-subset.zarr/

📦 新创建的测试数据:
  1. loom: test_data_all_formats/test_data.loom
  2. csv: test_data_all_formats/test_expression.csv
  3. csv_gz: test_data_all_formats/test_expression.csv.gz
  4. tsv: test_data_all_formats/test_expression.tsv
  5. tsv_gz: test_data_all_formats/test_expression.tsv.gz
  6. excel: test_data_all_formats/test_expression.xlsx
  7. mtx: test_data_all_formats/test_matrix.mtx
  8. mtx_gz: test_data_all_formats/test_matrix.mtx.gz
  9. hdf5: test_data_all_formats/test_data.hdf5
  10. umi_tools: test_data_all_formats/umi_tools_counts.tsv.gz
  11. custom_10x_mtx: test_data_all_formats/custom_10x_mtx/

🌐 需要网络下载:
  12. soft_gz: ~/.cache/scanpy-data/burczynski06/GDS1615_full.soft.gz
      (运行: import scanpy as sc; sc.datasets.burczynski06())
```

---

## 📝 各格式详细获取方式

### 5. Zarr格式

**状态：** ✅ Scanpy已有

**直接使用：**
```
/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x-10k-subset.zarr/
```

---

### 7. TSV/TXT格式

**状态：** ✅ Scanpy已有 + 📦 可创建更多

**Scanpy现有：**
```
/Users/warm/华大智造/TCGA/gdc/scanpy/src/scanpy/datasets/krumsiek11.txt
/Users/warm/华大智造/TCGA/gdc/scanpy/src/scanpy/datasets/toggleswitch.txt
```

**创建自定义TSV：**
```bash
python create_test_data.py
# 生成更大的TSV文件: test_data_all_formats/test_expression.tsv
```

---

### 8. Excel格式

**状态：** 📦 需创建 或 🌐 下载

**方式1：创建**
```bash
pip install openpyxl
python create_test_data.py
```

**方式2：使用Scanpy数据集**
```python
import scanpy as sc
adata = sc.datasets.moignard15()  # 会下载Excel格式数据
```

**数据位置：**
```
~/.cache/scanpy-data/moignard15/nbt.3154-S3.xlsx
```

---

### 9. MTX格式

**状态：** 📦 需创建

**创建：**
```bash
python create_test_data.py
```

**生成：**
```
test_data_all_formats/test_matrix.mtx
test_data_all_formats/test_matrix.mtx.gz
```

---

### 10. HDF5格式

**状态：** 📦 需创建

**创建：**
```bash
python create_test_data.py
```

**生成：**
```
test_data_all_formats/test_data.hdf5
```

---

### 11. SOFT.GZ格式

**状态：** 🌐 需下载

**方式1：Scanpy自动下载（推荐）**
```bash
python -c "import scanpy as sc; sc.datasets.burczynski06(); print('✓ 下载完成')"
```

**缓存位置：**
```
~/.cache/scanpy-data/burczynski06/GDS1615_full.soft.gz
```

**方式2：手动下载**
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS1nnn/GDS1615/soft/GDS1615_full.soft.gz
```

**GEO数据库：** https://www.ncbi.nlm.nih.gov/geo/

---

### 12. UMI-tools格式

**状态：** 📦 需创建

**创建：**
```bash
python create_test_data.py
```

**生成：**
```
test_data_all_formats/umi_tools_counts.tsv.gz
```

---

## ⚡ 最快测试方法

### 方法1：只测试Scanpy现有数据（1分钟）

```bash
bash quick_test.sh
```

**测试6种格式：**
1. 10X MTX v3
2. 10X MTX v2  
3. 10X H5
4. H5AD
5. TXT
6. Zarr

### 方法2：完整测试（5分钟）

```bash
# 创建所有格式
python create_test_data.py

# 运行测试
cd test_data_all_formats
bash run_all_format_tests.sh
```

**测试14种格式**（包括创建的格式）

---

## 📍 完整路径速查表

### Scanpy本地数据（绝对路径）

```bash
# 10X MTX v3
/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x_data/3.0.0/filtered_feature_bc_matrix/

# 10X MTX v2
/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x_data/1.2.0/filtered_gene_bc_matrices/hg19_chr21/

# 10X H5 v3
/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x_data/3.0.0/filtered_feature_bc_matrix.h5

# 10X H5 v2
/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x_data/1.2.0/filtered_gene_bc_matrices_h5.h5

# H5AD
/Users/warm/华大智造/TCGA/gdc/scanpy/src/scanpy/datasets/10x_pbmc68k_reduced.h5ad

# TXT
/Users/warm/华大智造/TCGA/gdc/scanpy/src/scanpy/datasets/krumsiek11.txt

# Zarr
/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x-10k-subset.zarr/
```

### 创建后的数据（相对路径）

```bash
# 运行: python create_test_data.py

# Loom
test_data_all_formats/test_data.loom

# CSV
test_data_all_formats/test_expression.csv
test_data_all_formats/test_expression.csv.gz

# TSV
test_data_all_formats/test_expression.tsv
test_data_all_formats/test_expression.tsv.gz

# Excel
test_data_all_formats/test_expression.xlsx

# MTX
test_data_all_formats/test_matrix.mtx
test_data_all_formats/test_matrix.mtx.gz

# HDF5
test_data_all_formats/test_data.hdf5

# UMI-tools
test_data_all_formats/umi_tools_counts.tsv.gz

# 自定义10X
test_data_all_formats/custom_10x_mtx/
```

---

## 🎓 总结

### 数据来源分布

- **Scanpy本地**: 7种格式（直接可用）
- **脚本创建**: 7种格式（运行create_test_data.py）
- **网络下载**: 1种格式（SOFT.GZ，可选）

### 推荐操作步骤

```bash
# 步骤1：快速测试（使用现有数据）
bash quick_test.sh

# 步骤2：创建完整数据（如需要）
python create_test_data.py

# 步骤3：下载GEO数据（可选）
python -c "import scanpy as sc; sc.datasets.burczynski06()"

# 步骤4：验证所有格式
cd test_data_all_formats
bash run_all_format_tests.sh
```

### 预期结果

所有14种格式都成功转换为统一的H5AD格式！

```
✓ 10x_mtx_v3.h5ad
✓ 10x_mtx_v2.h5ad
✓ 10x_h5.h5ad
✓ h5ad.h5ad
✓ loom.h5ad
✓ zarr.h5ad
✓ csv.h5ad
✓ tsv.h5ad
✓ excel.h5ad
✓ mtx.h5ad
✓ hdf5.h5ad
✓ soft_gz.h5ad
✓ umi_tools.h5ad
✓ custom_10x_mtx.h5ad
```

---

## 开始测试

```bash
# 最快方式（推荐）
bash quick_test.sh

# 查看结果
ls -lh quick_test_output/*.h5ad
```

