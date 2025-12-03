# DNB C4 格式特殊说明

## 概述

MGI DNBelab C4 平台使用 dnbc4tools 生成的输出数据采用10X兼容格式，但 `features.tsv` 文件存在格式差异。

---

## 格式差异

### 标准10X格式

**features.tsv** 包含 **2-3列**：

```
ENSG00000243485    MIR1302-2HG    Gene Expression
ENSG00000237613    FAM138A        Gene Expression
ENSG00000186092    OR4F5          Gene Expression
```

**列说明：**
- 第1列：基因ID（Ensembl ID）
- 第2列：基因符号（Gene Symbol）
- 第3列：特征类型（Feature Type，v3+才有）

### DNB C4格式

**features.tsv** 只有 **1列**：

```
Xkr4
Gm1992
Gm19938
Sox17
Lypla1
```

**列说明：**
- 只有基因名（Gene Symbol）
- 没有基因ID
- 没有特征类型

---

## scRNA-DataHub的处理方式

### 自动检测

工具会自动检测features.tsv的列数：

```python
# 读取第一行
first_line = f.readline().strip()
n_cols = len(first_line.split('\t'))

if n_cols == 1:
    # DNB C4格式：使用自定义读取方法
    adata = self._read_dnb_c4_mtx(path, ...)
else:
    # 标准10X格式：使用scanpy标准方法
    adata = sc.read_10x_mtx(path, ...)
```

### 自动补充缺失列

对于DNB C4格式，工具会自动补充：

```python
# 自动添加gene_ids列（与基因名相同）
adata.var['gene_ids'] = gene_names

# 自动添加feature_types列
adata.var['feature_types'] = ['Gene Expression'] * len(gene_names)
```

### 结果对比

| 字段 | 标准10X | DNB C4 | scRNA-DataHub处理后 |
|------|---------|--------|-------------------|
| `var_names` | 基因符号 | 基因符号 | 基因符号 ✅ |
| `var['gene_ids']` | Ensembl ID | ❌ 缺失 | 基因符号（自动填充）✅ |
| `var['feature_types']` | Feature类型 | ❌ 缺失 | 'Gene Expression'（自动填充）✅ |

---

## 使用方法

### 读取DNB C4数据

```bash
# 与标准10X格式完全相同的命令
python src/universal_reader.py \
  CNS1063416_brain/02.count/filter_matrix/ \
  -o CNS1063416_brain.h5ad \
  --sample-id CNS1063416_brain
```

**工具会自动：**
1. 检测features.tsv只有1列
2. 切换到DNB C4读取模式
3. 补充缺失的gene_ids和feature_types列
4. 输出标准的AnnData对象

### Python API使用

```python
from src.universal_reader import UniversalScRNAReader

reader = UniversalScRNAReader(verbose=True)

# 自动检测并处理DNB C4格式
adata = reader.read_auto(
    'CNS1063416_brain/02.count/filter_matrix/'
)

# 查看结果
print(adata.var.columns)
# Index(['gene_ids', 'feature_types'], dtype='object')

print(adata.var.head())
#           gene_ids      feature_types
# Xkr4      Xkr4         Gene Expression
# Gm1992    Gm1992       Gene Expression
# Gm19938   Gm19938      Gene Expression
```

---

## DNB C4 数据特点

### 文件结构

```
filter_matrix/
├── barcodes.tsv.gz      # 细胞条形码（标准格式）
├── features.tsv.gz      # 基因名（只有1列）⚠️
└── matrix.mtx.gz        # 表达矩阵（标准格式）
```

### 与10X的兼容性

| 组件 | 10X格式 | DNB C4格式 | 兼容性 |
|------|---------|-----------|--------|
| **barcodes.tsv** | 标准 | 标准 | ✅ 完全兼容 |
| **matrix.mtx** | 标准 | 标准 | ✅ 完全兼容 |
| **features.tsv** | 2-3列 | 1列 | ⚠️ 需要特殊处理 |

---

## 手动处理DNB C4数据

如果需要手动处理（不使用scRNA-DataHub）：

```python
import pandas as pd
import scipy.io as sio
import gzip
import anndata as ad

# 读取文件
path = 'filter_matrix/'

# 1. 读取barcodes
with gzip.open(f'{path}/barcodes.tsv.gz', 'rt') as f:
    barcodes = [line.strip() for line in f]

# 2. 读取features（只有1列）
with gzip.open(f'{path}/features.tsv.gz', 'rt') as f:
    gene_names = [line.strip() for line in f]

# 3. 读取matrix
matrix = sio.mmread(f'{path}/matrix.mtx.gz').T.tocsr()

# 4. 创建AnnData
adata = ad.AnnData(X=matrix)
adata.obs_names = barcodes
adata.var_names = gene_names

# 5. 补充缺失列（重要！）
adata.var['gene_ids'] = gene_names
adata.var['feature_types'] = ['Gene Expression'] * len(gene_names)

# 6. 处理重复基因名
adata.var_names_make_unique()

print(adata)
```

---

## 常见问题

### Q1: 为什么DNB C4的features.tsv只有1列？

**A:** DNB C4平台（华大智造）的dnbc4tools在生成10X兼容格式时，只输出基因名，不包含Ensembl ID等其他信息。这是该平台的设计选择。

### Q2: 缺少gene_ids会影响分析吗？

**A:** 不会。scRNA-DataHub会自动：
- 用基因名填充gene_ids列
- 添加feature_types列
- 确保与标准10X格式的兼容性

### Q3: 如何区分DNB C4和标准10X格式？

**A:** 工具会自动检测features.tsv的列数：
- 1列 → DNB C4格式
- 2列 → 10X v2格式
- 3列 → 10X v3格式

### Q4: DNB C4数据可以与10X数据合并吗？

**A:** 可以！处理后都是标准AnnData对象：

```python
from src.universal_reader import UniversalScRNAReader
import scanpy as sc

reader = UniversalScRNAReader()

# 读取DNB C4数据
adata_dnb = reader.read_auto('dnb_data/filter_matrix/')
adata_dnb.obs['platform'] = 'DNB_C4'

# 读取10X数据
adata_10x = reader.read_auto('10x_data/filtered_feature_bc_matrix/')
adata_10x.obs['platform'] = '10X'

# 合并
adata_combined = sc.concat([adata_dnb, adata_10x], join='outer')
```

---

## dnbc4tools输出目录结构

### 完整输出

```
CNS1063416_brain/
├── 01.data/
│   └── final_sorted.bam
├── 02.count/
│   ├── filter_matrix/          ← 这个目录
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz     ← 只有1列
│   │   └── matrix.mtx.gz
│   ├── raw_matrix/
│   └── singlecell.csv
└── 03.analysis/
    └── filter_feature.h5ad
```

### 推荐读取位置

```bash
# 方式1：读取filter_matrix（推荐）
python src/universal_reader.py \
  CNS1063416_brain/02.count/filter_matrix/ \
  -o output.h5ad

# 方式2：如果已经有h5ad，直接读取
python src/universal_reader.py \
  CNS1063416_brain/03.analysis/filter_feature.h5ad \
  -o output.h5ad
```

---

## 技术细节

### DNB C4格式读取流程

```
1. 检测features.tsv列数
   ↓
2. 如果是1列 → 使用_read_dnb_c4_mtx()
   ↓
3. 手动读取3个文件：
   - barcodes.tsv.gz
   - features.tsv.gz（1列）
   - matrix.mtx.gz
   ↓
4. 创建AnnData对象
   ↓
5. 补充缺失列：
   - var['gene_ids'] = 基因名
   - var['feature_types'] = 'Gene Expression'
   ↓
6. 处理重复基因名
   ↓
7. 返回标准AnnData对象
```

### 与标准10X的区别

| 步骤 | 标准10X | DNB C4 |
|------|---------|--------|
| **检测** | 调用sc.read_10x_mtx() | 检测列数→自定义读取 |
| **gene_ids** | 从文件读取 | 用基因名填充 |
| **feature_types** | 从文件读取 | 固定为'Gene Expression' |
| **最终结果** | AnnData对象 | AnnData对象（相同结构）|

---

## 测试DNB C4格式

```bash
# 使用DNB C4测试数据
python src/universal_reader.py \
  /storeData/ztron/wangrm/zzTCGA/pipeline/scRNA/results/CNS1063416_brain/02.count/filter_matrix \
  -o CNS1063416_brain.h5ad \
  --sample-id CNS1063416_brain

# 验证结果
python -c "
import scanpy as sc
adata = sc.read_h5ad('CNS1063416_brain.h5ad')
print(adata)
print('\nvar列:')
print(adata.var.columns.tolist())
print('\n前5个基因:')
print(adata.var.head())
"
```

**预期输出：**
```
AnnData object with n_obs × n_vars = XXXX × YYYY
    obs: 'sample_id'
    var: 'gene_ids', 'feature_types'

var列:
['gene_ids', 'feature_types']

前5个基因:
         gene_ids      feature_types
Xkr4     Xkr4         Gene Expression
Gm1992   Gm1992       Gene Expression
Gm19938  Gm19938      Gene Expression
Sox17    Sox17        Gene Expression
Lypla1   Lypla1       Gene Expression
```

---

## 总结

### DNB C4格式支持

- ✅ **自动检测** - 通过列数判断
- ✅ **自动处理** - 专门的读取方法
- ✅ **自动补充** - 添加缺失的列
- ✅ **完全兼容** - 输出标准AnnData对象

### 使用建议

1. **无需特殊操作** - 与10X数据使用相同的命令
2. **自动识别** - 工具会自动识别DNB C4格式
3. **透明处理** - 用户无感知的格式转换
4. **标准输出** - 输出与10X数据完全一致的H5AD

**DNB C4用户：直接使用，无需担心格式差异！** ✅

