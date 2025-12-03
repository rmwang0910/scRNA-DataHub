# scRNA-seq 数据格式与读取方式完整指南

## 概述

单细胞RNA测序（scRNA-seq）数据有多种格式，从原始测序数据（FASTQ）到处理后的表达矩阵（10X格式、H5AD等）。本文档详细介绍各种数据格式及其读取方式。

---

## 1. 原始测序数据格式

### 1.1 FASTQ 格式

**FASTQ** 是单细胞测序的原始数据格式，包含序列和质量分数。

#### 1.1.1 标准FASTQ格式

**文件结构：**
```
@read_id
ATCGATCGATCG...
+
IIIIIIIIIIII...
```

**字段说明：**
- `@read_id`：read标识符
- `ATCG...`：DNA序列
- `+`：分隔符（可选，后面可跟read_id）
- `IIII...`：质量分数（ASCII编码）

**文件命名：**
- `sample_R1.fastq.gz`：Read 1（正向）
- `sample_R2.fastq.gz`：Read 2（反向）
- `sample_I1.fastq.gz`：Index 1（可选）

#### 1.1.2 DNB C4 (DNBelab C Series) FASTQ格式

**DNB C4平台特点：**
- 使用DNB（DNA纳米球）技术
- 分离的cDNA和oligo文库

**文件结构：**
```
# cDNA文库（基因表达数据）
sample_cDNA_R1.fastq.gz    # cDNA Read 1
sample_cDNA_R2.fastq.gz    # cDNA Read 2

# oligo文库（细胞条形码和UMI）
sample_oligo_R1.fastq.gz   # oligo Read 1（包含细胞条形码）
sample_oligo_R2.fastq.gz   # oligo Read 2（包含UMI）
```

**读取方式：**

```python
import gzip

# 读取FASTQ文件
def read_fastq(fastq_file):
    """读取FASTQ文件"""
    reads = []
    with gzip.open(fastq_file, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            separator = f.readline().strip()
            quality = f.readline().strip()
            reads.append({
                'header': header,
                'sequence': sequence,
                'quality': quality
            })
    return reads

# 示例：读取DNB C4的cDNA R1文件
cDNA_R1_reads = read_fastq('sample_cDNA_R1.fastq.gz')
print(f"读取了 {len(cDNA_R1_reads)} 条reads")
```

**使用dnbc4tools处理：**

```bash
dnbc4tools rna run \
    --name sample \
    --cDNAfastq1 sample_cDNA_R1.fastq.gz \
    --cDNAfastq2 sample_cDNA_R2.fastq.gz \
    --oligofastq1 sample_oligo_R1.fastq.gz \
    --oligofastq2 sample_oligo_R2.fastq.gz \
    --genomeDir /path/to/genome \
    --threads 10
```

#### 1.1.3 10X Genomics FASTQ格式

**文件结构：**
```
# 10X Genomics v2/v3格式
sample_S1_L001_R1_001.fastq.gz  # Read 1（包含细胞条形码和UMI）
sample_S1_L001_R2_001.fastq.gz  # Read 2（基因序列）
sample_S1_L001_I1_001.fastq.gz  # Index（样本索引，可选）
```

**特点：**
- R1包含：细胞条形码（16bp）+ UMI（10-12bp）+ 基因序列
- R2包含：完整的基因序列

**处理工具：**
- Cell Ranger（10X官方工具）
- Alevin（Salmon的单细胞版本）
- STARsolo

---

## 2. 表达矩阵格式

### 2.1 10X Genomics格式（最常用）

**10X格式** 是单细胞数据最常用的标准格式，包含三个文件：

#### 2.1.1 文件组成

```
filtered_feature_bc_matrix/
├── barcodes.tsv.gz      # 细胞条形码列表
├── features.tsv.gz      # 基因信息（基因ID、基因名、类型）
└── matrix.mtx.gz        # 表达矩阵（稀疏矩阵格式）
```

#### 2.1.2 文件格式详解

**1. barcodes.tsv.gz**

**内容：** 每行一个细胞条形码

```
AAACCCAAGGAGAGTA-1
AAACCCAAGGCTTGCA-1
AAACCCAAGGTATGGT-1
AAACCCACAAGAAACT-1
...
```

**格式说明：**
- 每行一个条形码
- 格式：`{barcode}-{sample_id}`
- 例如：`AAACCCAAGGAGAGTA-1` 表示样本1的细胞

**2. features.tsv.gz**

**内容：** 基因信息（3列）

```
ENSG00000243485  MIR1302-2HG  Gene Expression
ENSG00000237613  FAM138A      Gene Expression
ENSG00000186092  OR4F5        Gene Expression
ENSG00000238009  AL627309.1   Gene Expression
...
```

**列说明：**
- **第1列**：基因ID（如Ensembl ID）
- **第2列**：基因名称（如基因符号）
- **第3列**：特征类型（通常是"Gene Expression"，也可能是"Antibody Capture"等）

**3. matrix.mtx.gz**

**内容：** 稀疏矩阵格式（Matrix Market格式）

```
%%MatrixMarket matrix coordinate integer general
%metadata_json: {"format_version": 2, "software_version": "3.0.0"}
27998 2700 2286884
1 1 1
1 2 1
1 3 2
2 1 1
...
```

**格式说明：**
- **第1行**：注释行（`%%MatrixMarket matrix coordinate integer general`）
- **第2行**：元数据（JSON格式，可选）
- **第3行**：`{n_genes} {n_cells} {n_entries}`（基因数、细胞数、非零条目数）
- **后续行**：`{gene_index} {cell_index} {count}`（基因索引、细胞索引、表达计数）

**注意：**
- 索引从1开始（不是0）
- 只存储非零值（稀疏矩阵）
- 基因索引对应features.tsv的行号
- 细胞索引对应barcodes.tsv的行号

#### 2.1.3 读取方式

**方法1：使用Scanpy（Python）**

```python
import scanpy as sc

# 读取10X格式数据
adata = sc.read_10x_mtx(
    path='filtered_feature_bc_matrix/',  # 包含三个文件的目录
    var_names='gene_symbols',            # 使用基因符号作为变量名
    cache=True                            # 缓存读取结果
)

# 或者直接指定目录
adata = sc.read_10x_mtx(
    path='filtered_feature_bc_matrix/',
    var_names='gene_symbols'
)

print(adata)
# AnnData object with n_obs × n_vars = 2700 × 27998
#     obs: None
#     var: 'gene_ids', 'feature_types'
#     uns: 'genome'
```

**方法2：使用Seurat（R）**

```r
library(Seurat)

# 读取10X格式数据
data <- Read10X(data.dir = "filtered_feature_bc_matrix/")

# 创建Seurat对象
seurat_obj <- CreateSeuratObject(counts = data, project = "sample1")

# 查看数据
seurat_obj
# An object of class Seurat 
# 27998 features across 2700 samples
```

**方法3：手动读取（Python）**

```python
import gzip
import pandas as pd
import scipy.io
import numpy as np

def read_10x_mtx_manual(matrix_dir):
    """手动读取10X格式数据"""
    
    # 读取barcodes
    with gzip.open(f'{matrix_dir}/barcodes.tsv.gz', 'rt') as f:
        barcodes = [line.strip() for line in f]
    
    # 读取features
    features = pd.read_csv(
        f'{matrix_dir}/features.tsv.gz',
        sep='\t',
        header=None,
        names=['gene_id', 'gene_symbol', 'feature_type']
    )
    
    # 读取matrix
    matrix = scipy.io.mmread(f'{matrix_dir}/matrix.mtx.gz')
    matrix = matrix.T  # 转置：从 (genes × cells) 转为 (cells × genes)
    
    # 创建AnnData对象
    import anndata as ad
    adata = ad.AnnData(X=matrix)
    adata.obs_names = barcodes
    adata.var_names = features['gene_symbol'].values
    adata.var['gene_ids'] = features['gene_id'].values
    adata.var['feature_types'] = features['feature_type'].values
    
    return adata

# 使用
adata = read_10x_mtx_manual('filtered_feature_bc_matrix/')
```

**方法4：读取压缩和非压缩格式**

```python
import scanpy as sc
import os

# 自动检测是否压缩
def read_10x_auto(matrix_dir):
    """自动检测并读取10X格式（支持压缩和非压缩）"""
    
    # 检查文件是否存在（支持.gz压缩）
    barcode_file = f'{matrix_dir}/barcodes.tsv'
    feature_file = f'{matrix_dir}/features.tsv'
    matrix_file = f'{matrix_dir}/matrix.mtx'
    
    if not os.path.exists(barcode_file):
        barcode_file += '.gz'
    if not os.path.exists(feature_file):
        feature_file += '.gz'
    if not os.path.exists(matrix_file):
        matrix_file += '.gz'
    
    # 使用Scanpy读取
    adata = sc.read_10x_mtx(
        path=matrix_dir,
        var_names='gene_symbols',
        cache=True
    )
    
    return adata
```

---

### 2.2 H5AD格式（AnnData）

**H5AD** 是Scanpy/AnnData的标准格式，是单细胞分析中最常用的存储格式。

#### 2.2.1 文件结构

```
data.h5ad
├── X (表达矩阵)
├── obs (细胞元数据)
├── var (基因元数据)
├── obsm (细胞多维数据：PCA、UMAP等)
├── varm (基因多维数据)
├── uns (非结构化数据：参数、颜色等)
└── layers (其他层：counts、normalized等)
```

#### 2.2.2 读取方式

**使用Scanpy读取：**

```python
import scanpy as sc

# 读取H5AD文件
adata = sc.read_h5ad('data.h5ad')

# 查看数据结构
print(adata)
# AnnData object with n_obs × n_vars = 2700 × 27998
#     obs: 'n_genes', 'n_counts', 'leiden', 'SCSA_celltype'
#     var: 'gene_ids', 'n_cells', 'highly_variable'
#     obsm: 'X_pca', 'X_umap'
#     uns: 'pca', 'neighbors', 'umap', 'leiden'

# 查看细胞元数据
print(adata.obs.head())

# 查看基因元数据
print(adata.var.head())

# 查看表达矩阵
print(adata.X.shape)  # (n_cells, n_genes)
print(adata.X[0:5, 0:5])  # 前5个细胞，前5个基因
```

**使用h5py直接读取（高级）：**

```python
import h5py
import numpy as np

# 直接读取H5AD文件（不转换为AnnData）
with h5py.File('data.h5ad', 'r') as f:
    # 读取表达矩阵
    X = f['X'][:]
    
    # 读取细胞信息
    obs_names = [s.decode() for s in f['obs']['_index'][:]]
    
    # 读取基因信息
    var_names = [s.decode() for s in f['var']['_index'][:]]
    
    print(f"表达矩阵形状: {X.shape}")
    print(f"细胞数: {len(obs_names)}")
    print(f"基因数: {len(var_names)}")
```

---

### 2.3 H5格式（10X HDF5）

**H5格式** 是10X Genomics的HDF5格式，包含完整的表达矩阵和元数据。

#### 2.3.1 文件结构

```
filtered_feature_bc_matrix.h5
├── matrix (组)
│   ├── barcodes (数据集)
│   ├── data (数据集：表达值)
│   ├── indices (数据集：基因索引)
│   ├── indptr (数据集：细胞索引指针)
│   └── shape (数据集：矩阵维度)
└── matrix (组)
    └── features (组)
        ├── _all_tag_keys (数据集)
        ├── feature_type (数据集)
        ├── genome (数据集)
        ├── id (数据集：基因ID)
        └── name (数据集：基因名)
```

#### 2.3.2 读取方式

**使用Scanpy读取：**

```python
import scanpy as sc

# 读取10X H5文件
adata = sc.read_10x_h5('filtered_feature_bc_matrix.h5')

print(adata)
```

**使用h5py读取：**

```python
import h5py
import scipy.sparse as sp

def read_10x_h5_manual(h5_file):
    """手动读取10X H5文件"""
    with h5py.File(h5_file, 'r') as f:
        # 读取矩阵信息
        matrix_group = f['matrix']
        
        # 读取稀疏矩阵
        data = matrix_group['data'][:]
        indices = matrix_group['indices'][:]
        indptr = matrix_group['indptr'][:]
        shape = matrix_group['shape'][:]
        
        # 构建稀疏矩阵
        matrix = sp.csr_matrix((data, indices, indptr), shape=shape)
        matrix = matrix.T  # 转置
        
        # 读取barcodes
        barcodes = [s.decode() for s in matrix_group['barcodes'][:]]
        
        # 读取features
        features_group = matrix_group['features']
        gene_ids = [s.decode() for s in features_group['id'][:]]
        gene_names = [s.decode() for s in features_group['name'][:]]
        
        # 创建AnnData
        import anndata as ad
        adata = ad.AnnData(X=matrix)
        adata.obs_names = barcodes
        adata.var_names = gene_names
        adata.var['gene_ids'] = gene_ids
        
        return adata

# 使用
adata = read_10x_h5_manual('filtered_feature_bc_matrix.h5')
```

---

### 2.4 CSV/TSV格式

**CSV/TSV格式** 是简单的表格格式，基因为行，细胞为列。

#### 2.4.1 文件格式

```
gene_id,gene_name,cell1,cell2,cell3,...
GAPDH,GAPDH,100,150,120,...
ACTB,ACTB,200,180,190,...
TP53,TP53,50,60,55,...
...
```

#### 2.4.2 读取方式

**使用Pandas + AnnData：**

```python
import pandas as pd
import anndata as ad

# 读取CSV文件
df = pd.read_csv('expression_matrix.csv', index_col=0)

# 如果第一列是基因ID，第二列是基因名
if 'gene_name' in df.columns:
    gene_names = df['gene_name'].values
    df = df.drop('gene_name', axis=1)
else:
    gene_names = df.index.values

# 创建AnnData对象
adata = ad.AnnData(X=df.T)  # 转置：细胞×基因
adata.var_names = gene_names
adata.obs_names = df.columns

print(adata)
```

**使用Scanpy：**

```python
import scanpy as sc
import pandas as pd

# 读取CSV
df = pd.read_csv('expression_matrix.csv', index_col=0)

# 转换为AnnData
adata = sc.AnnData(df.T)
adata.var_names = df.index
adata.obs_names = df.columns
```

---

### 2.5 Loom格式

**Loom格式** 是专门为单细胞数据设计的HDF5格式。

#### 2.5.1 文件结构

```
data.loom
├── matrix (数据集：表达矩阵)
├── row_attrs (组：基因属性)
│   ├── Gene (数据集：基因名)
│   └── Accession (数据集：基因ID)
└── col_attrs (组：细胞属性)
    ├── CellID (数据集：细胞ID)
    └── ClusterID (数据集：聚类ID)
```

#### 2.5.2 读取方式

**使用Scanpy：**

```python
import scanpy as sc

# 读取Loom文件
adata = sc.read_loom('data.loom')

print(adata)
```

**使用loompy：**

```python
import loompy

# 读取Loom文件
with loompy.connect('data.loom') as ds:
    # 获取表达矩阵
    matrix = ds[:, :]
    
    # 获取基因信息
    gene_names = ds.ra['Gene']
    
    # 获取细胞信息
    cell_ids = ds.ca['CellID']
    
    # 创建AnnData
    import anndata as ad
    adata = ad.AnnData(X=matrix.T)
    adata.var_names = gene_names
    adata.obs_names = cell_ids
```

---

## 3. 其他数据格式

### 3.1 RDS格式（R对象）

**RDS格式** 是R语言的序列化对象格式，常用于存储Seurat对象。

#### 3.3.1 读取方式

**在R中读取：**

```r
library(Seurat)

# 读取RDS文件
seurat_obj <- readRDS('data.rds')

# 转换为10X格式（如果需要）
SaveH5Seurat(seurat_obj, filename = "data.h5seurat")
```

**在Python中读取（需要rpy2）：**

```python
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()

# 读取RDS
ro.r('library(Seurat)')
seurat_obj = ro.r('readRDS("data.rds")')

# 提取表达矩阵
counts = ro.r('GetAssayData')(seurat_obj, slot='counts', assay='RNA')
```

**推荐：先转换为H5AD**

```r
# 在R中转换
library(Seurat)
library(SeuratDisk)

# 读取RDS
seurat_obj <- readRDS('data.rds')

# 转换为H5AD
SaveH5Seurat(seurat_obj, filename = "data.h5seurat")
Convert("data.h5seurat", dest = "h5ad")
```

然后在Python中读取：

```python
import scanpy as sc
adata = sc.read_h5ad('data.h5ad')
```

---

### 3.2 H5Seurat格式

**H5Seurat格式** 是Seurat的HDF5格式。

#### 3.2.1 读取方式

**在R中读取：**

```r
library(Seurat)
library(SeuratDisk)

# 读取H5Seurat
seurat_obj <- LoadH5Seurat("data.h5seurat")

# 转换为H5AD
Convert("data.h5seurat", dest = "h5ad")
```

**在Python中读取（需要先转换）：**

```python
# 需要先在R中转换为H5AD
# 然后使用Scanpy读取
import scanpy as sc
adata = sc.read_h5ad('data.h5ad')
```

---

### 3.3 Zarr格式

**Zarr格式** 是一种云原生的数组存储格式，适合大规模数据。

#### 3.3.1 文件结构

```
data.zarr/
├── .zattrs          # 属性
├── .zgroup          # 组信息
├── X/               # 表达矩阵
├── obs/             # 细胞元数据
├── var/             # 基因元数据
└── ...
```

#### 3.3.2 读取方式

**使用Scanpy读取：**

```python
import scanpy as sc

# 读取Zarr格式
adata = sc.read_zarr('data.zarr')

print(adata)
```

**写入Zarr格式：**

```python
# 保存为Zarr格式（适合大数据）
adata.write_zarr('data.zarr')
```

**特点：**
- 支持云存储（S3、GCS等）
- 支持并行读写
- 适合超大规模数据（百万细胞级别）

---

### 3.4 UMI-tools格式

**UMI-tools格式** 是UMI去重后的count表格式。

#### 3.4.1 读取方式

**使用Scanpy读取：**

```python
import scanpy as sc

# 读取UMI-tools输出的计数表
adata = sc.read_umi_tools('counts.tsv.gz')

print(adata)
```

**文件格式示例：**
```
gene    cell1    cell2    cell3    ...
GAPDH   10       15       12       ...
ACTB    20       18       19       ...
TP53    5        6        5        ...
```

---

### 3.5 Excel格式

**Excel格式** 适合小规模数据或手动整理的数据。

#### 3.5.1 读取方式

**使用Scanpy读取：**

```python
import scanpy as sc

# 读取Excel文件（需要指定sheet名称）
adata = sc.read_excel('data.xlsx', sheet='expression_matrix')

print(adata)
```

**注意：**
- 不推荐用于大数据（Excel行数限制）
- 读写速度慢
- 仅适合小规模数据或示例数据

---

### 3.6 HDF格式

**HDF格式** 是通用的HDF5数据格式（非AnnData专用）。

#### 3.6.1 读取方式

**使用Scanpy读取：**

```python
import scanpy as sc

# 读取HDF文件（需要指定key）
adata = sc.read_hdf('data.h5', key='data')

print(adata)
```

---

### 3.7 SOFT.GZ格式（GEO数据）

**SOFT.GZ格式** 是NCBI GEO数据库的标准格式。

#### 3.7.1 文件结构

SOFT格式包含：
- 样本信息
- 基因表达数据
- 实验描述

#### 3.7.2 读取方式

**使用Scanpy读取：**

```python
import scanpy as sc

# 读取SOFT格式文件（通常从GEO下载）
adata = sc.read('GSE12345_series_matrix.txt.gz')

# 或者明确指定格式
adata = sc.read('GSE12345_family.soft.gz', ext='soft.gz')

print(adata)
```

**应用场景：**
- 直接从GEO数据库下载的微阵列或RNA-seq数据
- 已发表数据的重分析

**参考：**
- GEO SOFT格式文档：https://www.ncbi.nlm.nih.gov/geo/info/soft.html

---

### 3.3 BAM/SAM格式

**BAM/SAM格式** 是比对后的序列数据，需要进一步处理。

#### 3.3.1 处理方式

**使用Cell Ranger：**

```bash
# Cell Ranger count输出10X格式
cellranger count \
    --id=sample1 \
    --transcriptome=/path/to/ref \
    --fastqs=/path/to/fastq \
    --sample=sample1
```

**使用STARsolo：**

```bash
# STARsolo直接输出10X格式
STAR \
    --runMode alignReads \
    --genomeDir /path/to/genome \
    --readFilesIn sample_R1.fastq.gz sample_R2.fastq.gz \
    --readFilesCommand zcat \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist /path/to/barcodes.tsv \
    --outFileNamePrefix sample_
```

---

## 4. 数据格式转换

### 4.1 10X格式 → H5AD

```python
import scanpy as sc

# 读取10X格式
adata = sc.read_10x_mtx('filtered_feature_bc_matrix/')

# 保存为H5AD
adata.write('data.h5ad')
```

### 4.2 H5AD → 10X格式

```python
import scanpy as sc
import pandas as pd
import scipy.io
import gzip

def save_10x_format(adata, output_dir):
    """将AnnData保存为10X格式"""
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存barcodes
    with gzip.open(f'{output_dir}/barcodes.tsv.gz', 'wt') as f:
        for barcode in adata.obs_names:
            f.write(f'{barcode}\n')
    
    # 保存features
    features_df = pd.DataFrame({
        'gene_id': adata.var.get('gene_ids', adata.var_names),
        'gene_symbol': adata.var_names,
        'feature_type': adata.var.get('feature_types', ['Gene Expression'] * len(adata.var_names))
    })
    features_df.to_csv(f'{output_dir}/features.tsv.gz', sep='\t', header=False, index=False, compression='gzip')
    
    # 保存matrix
    matrix = adata.X.T  # 转置为 genes × cells
    scipy.io.mmwrite(f'{output_dir}/matrix.mtx', matrix)
    
    # 压缩matrix
    import subprocess
    subprocess.run(['gzip', f'{output_dir}/matrix.mtx'])

# 使用
adata = sc.read_h5ad('data.h5ad')
save_10x_format(adata, 'output_10x/')
```

### 4.3 Seurat对象 → H5AD

**在R中：**

```r
library(Seurat)
library(SeuratDisk)

# 读取Seurat对象
seurat_obj <- readRDS('data.rds')

# 转换为H5AD
SaveH5Seurat(seurat_obj, filename = "data.h5seurat")
Convert("data.h5seurat", dest = "h5ad")
```

**在Python中读取：**

```python
import scanpy as sc
adata = sc.read_h5ad('data.h5ad')
```

---

## 5. Scanpy支持的所有数据格式总结

### 5.1 Scanpy完整读取函数列表

根据Scanpy源码，以下是所有支持的读取函数：

| 函数名 | 格式 | 说明 | 示例代码 |
|--------|------|------|---------|
| `sc.read()` | 通用 | 自动识别格式 | `adata = sc.read('data.h5ad')` |
| `sc.read_10x_mtx()` | 10X MTX | 10X标准格式（3个文件） | `adata = sc.read_10x_mtx('filtered_feature_bc_matrix/')` |
| `sc.read_10x_h5()` | 10X H5 | 10X HDF5格式 | `adata = sc.read_10x_h5('filtered_feature_bc_matrix.h5')` |
| `sc.read_h5ad()` | H5AD | AnnData标准格式 | `adata = sc.read_h5ad('data.h5ad')` |
| `sc.read_loom()` | Loom | Loom格式 | `adata = sc.read_loom('data.loom')` |
| `sc.read_csv()` | CSV | 逗号分隔 | `adata = sc.read_csv('data.csv')` |
| `sc.read_text()` | TXT/TSV | 文本/制表符分隔 | `adata = sc.read_text('data.tsv', delimiter='\t')` |
| `sc.read_excel()` | Excel | Excel文件 | `adata = sc.read_excel('data.xlsx', sheet='Sheet1')` |
| `sc.read_mtx()` | MTX | 单个mtx文件 | `adata = sc.read_mtx('matrix.mtx')` |
| `sc.read_hdf()` | HDF5 | HDF5格式 | `adata = sc.read_hdf('data.h5', key='data')` |
| `sc.read_umi_tools()` | UMI-tools | UMI-tools计数表 | `adata = sc.read_umi_tools('counts.tsv.gz')` |
| `sc.read_zarr()` | Zarr | Zarr格式 | `adata = sc.read_zarr('data.zarr')` |
| `sc.read_visium()` | Visium | 空间转录组 | `adata = sc.read_visium('spatial_data/')` |

### 5.2 支持的文件扩展名

根据Scanpy源码 `avail_exts` 变量，支持的扩展名：

```python
# 文本格式
text_exts = {"csv", "tsv", "tab", "data", "txt"}

# 所有支持的扩展名
avail_exts = {
    "anndata",      # 通用AnnData格式
    "xlsx",         # Excel格式
    "h5",           # HDF5格式
    "h5ad",         # AnnData HDF5格式
    "zarr",         # Zarr格式
    "mtx",          # Matrix Market格式
    "mtx.gz",       # 压缩的Matrix Market格式
    "soft.gz",      # GEO SOFT格式
    "loom",         # Loom格式
    "csv",          # 逗号分隔
    "tsv",          # 制表符分隔
    "tab",          # 制表符分隔
    "data",         # 通用数据格式
    "txt"           # 文本格式
}
```

**注意：**
- 文本格式（csv, tsv, tab, data, txt）还支持 `.gz` 和 `.bz2` 压缩
- 例如：`data.csv.gz`、`matrix.tsv.bz2` 都是有效的

### 5.3 完整数据格式总结表

| 格式 | 文件扩展名 | 主要用途 | 读取函数 | 优点 | 缺点 | Scanpy支持 |
|------|-----------|---------|---------|------|------|-----------|
| **10X MTX** | `.mtx.gz`, `.tsv.gz` | 标准表达矩阵 | `read_10x_mtx()` | 标准格式、广泛支持 | 需要三个文件 | ✅ 完全支持 |
| **10X H5** | `.h5` | 10X HDF5 | `read_10x_h5()` | 单个文件、完整元数据 | 文件较大 | ✅ 完全支持 |
| **H5AD** | `.h5ad` | Scanpy标准 | `read_h5ad()` | 包含完整分析结果 | Python专用 | ✅ 原生格式 |
| **Zarr** | `.zarr/` | 云原生大数据 | `read_zarr()` | 云存储、并行读写 | 较新的格式 | ✅ 完全支持 |
| **Loom** | `.loom` | 单细胞存储 | `read_loom()` | 专门设计、高效 | 使用较少 | ✅ 完全支持 |
| **CSV** | `.csv`, `.csv.gz` | 简单表格 | `read_csv()` | 易读、通用 | 文件大、无元数据 | ✅ 完全支持 |
| **TSV/TXT** | `.tsv`, `.txt`, `.tab` | 文本表格 | `read_text()` | 通用格式 | 文件大 | ✅ 完全支持 |
| **Excel** | `.xlsx`, `.xls` | Excel表格 | `read_excel()` | 手动编辑友好 | 行数限制、慢 | ✅ 完全支持 |
| **MTX** | `.mtx`, `.mtx.gz` | 稀疏矩阵 | `read_mtx()` | 节省空间 | 需要额外文件 | ✅ 完全支持 |
| **HDF5** | `.h5`, `.hdf5` | 通用HDF5 | `read_hdf()` | 通用格式 | 需要指定key | ✅ 完全支持 |
| **UMI-tools** | `.tsv.gz` | UMI计数表 | `read_umi_tools()` | UMI去重专用 | 特定格式 | ✅ 完全支持 |
| **SOFT.GZ** | `.soft.gz` | GEO数据 | `read()` | GEO直接下载 | 旧格式 | ✅ 完全支持 |
| **Visium** | 目录 | 空间转录组 | `read_visium()` | 包含空间信息 | 已弃用 | ⚠️ 已弃用，用squidpy |
| **RDS** | `.rds` | R对象 | 需要R转换 | R生态完整 | 需要R环境 | ❌ 需要转换 |
| **H5Seurat** | `.h5seurat` | Seurat HDF5 | 需要R转换 | Seurat专用 | 需要转换 | ❌ 需要转换 |
| **FASTQ** | `.fastq.gz` | 原始测序 | 需要比对工具 | 原始数据 | 需要处理 | ❌ 需要预处理 |

---

## 6. Scanpy读取参数详解

### 6.1 read_10x_mtx() 参数详解

```python
sc.read_10x_mtx(
    path='filtered_feature_bc_matrix/',  # 数据目录
    var_names='gene_symbols',            # 使用基因符号 或 'gene_ids'
    make_unique=True,                     # 使基因名唯一
    cache=True,                           # 创建h5ad缓存（加速后续读取）
    cache_compression='gzip',             # 缓存压缩方式
    gex_only=True,                        # 只保留基因表达数据
    prefix=None,                          # 文件名前缀
    compressed=True                       # 是否期望.gz压缩（Cell Ranger v3+）
)
```

**参数说明：**

- **`var_names`**：选择基因名称列
  - `'gene_symbols'`：使用基因符号（如GAPDH、TP53）
  - `'gene_ids'`：使用基因ID（如ENSG00000111640）
  
- **`make_unique`**：处理重复基因名
  - `True`：自动添加后缀（如GAPDH-1, GAPDH-2）
  - `False`：保留重复名称（可能导致错误）

- **`cache`**：是否使用缓存
  - `True`：首次读取时创建 `.h5ad` 缓存，后续读取快10-100倍
  - `False`：每次都从原始文件读取

- **`gex_only`**：是否只保留基因表达数据
  - `True`：过滤掉 Antibody Capture、CRISPR Guide Capture 等
  - `False`：保留所有特征类型

- **`prefix`**：文件名前缀
  - 适用于自定义文件名，如 `sample1_matrix.mtx.gz`，则设置 `prefix='sample1_'`

- **`compressed`**：是否期望文件被压缩
  - `True`：期望 Cell Ranger v3+ 格式（.gz压缩）
  - `False`：期望 STARsolo 输出（未压缩）

**示例：读取STARsolo输出**

```python
# STARsolo输出是未压缩的
adata = sc.read_10x_mtx(
    'Solo.out/Gene/filtered/',
    var_names='gene_symbols',
    compressed=False  # STARsolo输出未压缩
)
```

### 6.2 read_10x_h5() 参数详解

```python
sc.read_10x_h5(
    'filtered_feature_bc_matrix.h5',
    genome=None,      # 基因组名称（多基因组时需要）
    gex_only=True     # 只保留基因表达数据
)
```

**参数说明：**

- **`genome`**：基因组过滤
  - 对于包含多个基因组的h5文件（如人鼠混合样本），需要指定
  - 例如：`genome='GRCh38'` 或 `genome='mm10'`

- **`gex_only`**：特征类型过滤
  - `True`：只保留 `feature_types == 'Gene Expression'`
  - `False`：保留所有类型（包括抗体、CRISPR等）

**示例：读取多基因组数据**

```python
# 人鼠混合样本，只读取人基因组
adata_human = sc.read_10x_h5(
    'filtered_feature_bc_matrix.h5',
    genome='GRCh38',
    gex_only=True
)

# 读取小鼠基因组
adata_mouse = sc.read_10x_h5(
    'filtered_feature_bc_matrix.h5',
    genome='mm10',
    gex_only=True
)
```

### 6.3 read() 通用函数参数详解

```python
sc.read(
    'data.csv',                  # 文件路径
    backed=None,                 # 'r' 或 'r+' 启用backed模式
    sheet=None,                  # Excel/HDF5的sheet名称
    ext=None,                    # 文件扩展名（自动检测）
    delimiter=None,              # 文本文件分隔符
    first_column_names=False,    # 第一列是否是行名
    backup_url=None,             # 备份下载URL
    cache=False,                 # 是否使用缓存
    cache_compression='gzip'     # 缓存压缩方式
)
```

**backed模式说明：**

```python
# 对于超大数据，使用backed模式（不全部加载到内存）
adata = sc.read_h5ad('large_data.h5ad', backed='r')

# 只读取部分数据
subset = adata[:1000, :500]  # 前1000个细胞，前500个基因
subset = subset.to_memory()   # 转为内存模式

# 如果需要修改，使用 'r+' 模式
adata = sc.read_h5ad('data.h5ad', backed='r+')
```

### 6.4 文本格式读取详解

**CSV格式：**

```python
# 标准CSV（逗号分隔）
adata = sc.read_csv('expression.csv')

# 压缩CSV
adata = sc.read_csv('expression.csv.gz')

# 指定第一列为行名
adata = sc.read_csv('expression.csv', first_column_names=True)
```

**TSV格式：**

```python
# 制表符分隔
adata = sc.read_text('expression.tsv', delimiter='\t')

# 或者使用简化写法
adata = sc.read('expression.tsv')  # 自动识别为text格式

# 空格分隔
adata = sc.read_text('expression.txt', delimiter=' ')
```

**自定义分隔符：**

```python
# 使用自定义分隔符
adata = sc.read_text('data.txt', delimiter='|')
```

---

## 7. Scanpy格式选择指南

### 7.1 根据数据来源选择读取方法

| 数据来源 | 文件格式 | 推荐读取方法 | 代码示例 |
|---------|---------|------------|---------|
| **10X Genomics** | MTX目录 | `read_10x_mtx()` | `sc.read_10x_mtx('filtered_feature_bc_matrix/')` |
| **10X Genomics** | H5文件 | `read_10x_h5()` | `sc.read_10x_h5('filtered_feature_bc_matrix.h5')` |
| **DNB C4** | dnbc4tools输出 | `read_10x_mtx()` | `sc.read_10x_mtx('filter_matrix/')` |
| **STARsolo** | MTX目录 | `read_10x_mtx(..., compressed=False)` | `sc.read_10x_mtx('Solo.out/Gene/filtered/', compressed=False)` |
| **Alevin** | Alevin输出 | 转换为10X | 需要alevin-fry转换 |
| **Scanpy分析结果** | H5AD | `read_h5ad()` | `sc.read_h5ad('processed.h5ad')` |
| **Seurat对象** | RDS | 先在R转换 | R: `Convert()` → Python: `read_h5ad()` |
| **GEO数据库** | SOFT.GZ | `read()` | `sc.read('GSE12345_family.soft.gz')` |
| **自定义矩阵** | CSV/TSV | `read_csv()` / `read_text()` | `sc.read_csv('matrix.csv')` |
| **Loom文件** | Loom | `read_loom()` | `sc.read_loom('data.loom')` |
| **大规模数据** | Zarr | `read_zarr()` | `sc.read_zarr('data.zarr')` |

### 7.2 根据数据规模选择格式

| 数据规模 | 推荐格式 | 原因 |
|---------|---------|------|
| **小数据** (< 10K cells) | CSV/TSV/Excel | 易于查看和编辑 |
| **中等数据** (10K-100K cells) | H5AD | 读写速度快，包含完整元数据 |
| **大数据** (100K-1M cells) | H5AD + 稀疏矩阵 | 节省内存，使用backed模式 |
| **超大数据** (> 1M cells) | Zarr | 支持云存储和并行处理 |

### 7.3 快速读取技巧

#### 7.3.1 使用缓存加速读取

```python
# 第一次读取（慢）
adata = sc.read_10x_mtx(
    'filtered_feature_bc_matrix/',
    cache=True  # 创建h5ad缓存
)

# 后续读取（快10-100倍）
adata = sc.read_10x_mtx(
    'filtered_feature_bc_matrix/',
    cache=True  # 从缓存读取
)
```

#### 7.3.2 使用backed模式处理大数据

```python
# 不全部加载到内存
adata = sc.read_h5ad('large_data.h5ad', backed='r')

# 只加载需要的部分
subset = adata[:10000, :].to_memory()  # 前10000个细胞

# 或者分批处理
for i in range(0, adata.n_obs, 10000):
    subset = adata[i:i+10000, :].to_memory()
    # 处理subset...
```

#### 7.3.3 STARsolo输出读取

```python
# STARsolo输出是未压缩的10X格式
adata = sc.read_10x_mtx(
    'Solo.out/Gene/filtered/',
    var_names='gene_symbols',
    compressed=False  # 关键：STARsolo不压缩
)
```

---

## 8. 推荐工作流程

### 8.1 从FASTQ开始

```
FASTQ文件
  ↓
[Cell Ranger / STARsolo / dnbc4tools]
  ↓
10X格式 (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)
  ↓
[Scanpy读取 + cache=True]
  ↓
H5AD格式 (用于后续分析)
```

### 8.2 从10X格式开始

```
10X格式
  ↓
[Scanpy读取 sc.read_10x_mtx()]
  ↓
H5AD格式
  ↓
[分析流程]
  ↓
保存为H5AD (包含所有分析结果)
```

### 8.3 从其他格式转换

```
其他格式 (RDS, CSV, Loom等)
  ↓
[转换为H5AD]
  ↓
[分析流程]
  ↓
保存为H5AD
```

---

## 7. Scanpy格式选择指南

### 7.1 根据数据来源选择读取方法

| 数据来源 | 文件格式 | 推荐读取方法 | 代码示例 |
|---------|---------|------------|---------|
| **Cell Ranger v3+** | MTX + gz | `read_10x_mtx()` | `sc.read_10x_mtx('filtered_feature_bc_matrix/')` |
| **Cell Ranger v2-** | MTX (无gz) | `read_10x_mtx()` | 自动识别legacy格式 |
| **STARsolo** | MTX (无gz) | `read_10x_mtx()` | `sc.read_10x_mtx('Solo.out/Gene/filtered/', compressed=False)` |
| **dnbc4tools** | MTX + gz | `read_10x_mtx()` | `sc.read_10x_mtx('filter_matrix/')` |
| **Alevin-fry** | MTX + gz | `read_10x_mtx()` | 标准格式 |
| **Scanpy分析** | H5AD | `read_h5ad()` | `sc.read_h5ad('processed.h5ad')` |
| **Seurat对象** | RDS | 先在R转换 | R: `Convert()` → Python: `read_h5ad()` |
| **GEO数据库** | SOFT.GZ | `read()` | `sc.read('GSE12345_family.soft.gz')` |
| **自定义矩阵** | CSV/TSV | `read_csv()` / `read_text()` | `sc.read_csv('matrix.csv')` |

### 7.2 完整Scanpy支持格式列表

根据Scanpy源码，支持的所有文件扩展名：

**文本格式：**
- `.csv`, `.csv.gz`, `.csv.bz2`
- `.tsv`, `.tsv.gz`, `.tsv.bz2`
- `.tab`, `.tab.gz`, `.tab.bz2`
- `.txt`, `.txt.gz`, `.txt.bz2`
- `.data`

**HDF5格式：**
- `.h5ad` (AnnData专用)
- `.h5` (通用HDF5)
- `.hdf5`

**其他格式：**
- `.loom`
- `.zarr/` (目录)
- `.mtx`, `.mtx.gz`
- `.soft.gz` (GEO格式)
- `.xlsx`, `.xls` (Excel)
- `.anndata` (旧版本)

---

## 10. 常见问题

### Q1: 如何判断数据格式？

**检查文件：**
```python
import os

def detect_format(data_path):
    """检测数据格式"""
    if os.path.isdir(data_path):
        # 检查是否是10X格式
        files = os.listdir(data_path)
        if 'matrix.mtx' in files or 'matrix.mtx.gz' in files:
            if 'barcodes.tsv' in files or 'barcodes.tsv.gz' in files:
                if 'features.tsv' in files or 'features.tsv.gz' in files:
                    return '10X格式'
    elif data_path.endswith('.h5ad'):
        return 'H5AD格式'
    elif data_path.endswith('.h5'):
        return 'H5格式'
    elif data_path.endswith('.loom'):
        return 'Loom格式'
    elif data_path.endswith('.csv') or data_path.endswith('.tsv'):
        return 'CSV/TSV格式'
    elif data_path.endswith('.rds'):
        return 'RDS格式'
    return '未知格式'

# 使用
format_type = detect_format('data/')
print(f"数据格式: {format_type}")
```

### Q2: 如何处理压缩文件？

**Scanpy自动处理压缩：**
```python
import scanpy as sc

# 自动识别和解压.gz和.bz2压缩
adata = sc.read_10x_mtx('filtered_feature_bc_matrix/')  # 自动处理.gz
adata = sc.read_csv('data.csv.gz')                       # 自动解压
adata = sc.read_text('data.tsv.bz2')                     # 自动解压
```

### Q3: 内存不足怎么办？

**方法1：使用稀疏矩阵**
```python
import scanpy as sc
import scipy.sparse as sp

# 确保使用稀疏矩阵
adata = sc.read_10x_mtx('filtered_feature_bc_matrix/')
print(type(adata.X))  # 应该是 scipy.sparse.csr_matrix

# 如果不是，转换为稀疏矩阵
if not sp.issparse(adata.X):
    adata.X = sp.csr_matrix(adata.X)
```

**方法2：使用backed模式**
```python
# 对于超大数据，使用backed模式
adata = sc.read_h5ad('large_data.h5ad', backed='r')

# 只加载需要的部分到内存
subset = adata[:10000, :].to_memory()
```

**方法3：使用Zarr格式**
```python
# Zarr格式支持按需加载
adata = sc.read_zarr('data.zarr')
```

### Q4: STARsolo输出怎么读取？

**关键：设置compressed=False**
```python
import scanpy as sc

# STARsolo输出是未压缩的10X格式
adata = sc.read_10x_mtx(
    'Solo.out/Gene/filtered/',
    var_names='gene_symbols',
    compressed=False  # 必须设置为False
)
```

### Q5: 读取10X数据时如何选择基因名称？

**选择gene_symbols或gene_ids：**
```python
# 方法1：使用基因符号（推荐）
adata = sc.read_10x_mtx(
    'filtered_feature_bc_matrix/',
    var_names='gene_symbols'  # 使用GAPDH, TP53等
)

# 方法2：使用基因ID
adata = sc.read_10x_mtx(
    'filtered_feature_bc_matrix/',
    var_names='gene_ids'  # 使用ENSG00000111640等
)

# 读取后互换
# 如果想换一下
adata.var['gene_ids_old'] = adata.var['gene_ids']
adata.var['gene_ids'] = adata.var_names
adata.var_names = adata.var['gene_symbols']
```

### Q6: 如何读取包含多种特征类型的数据？

**10X Genomics支持多种特征类型：**
```python
# Gene Expression + Antibody Capture (CITE-seq)
adata = sc.read_10x_h5(
    'filtered_feature_bc_matrix.h5',
    gex_only=False  # 保留所有特征类型
)

# 查看特征类型
print(adata.var['feature_types'].value_counts())
# Gene Expression       20000
# Antibody Capture      50

# 分离不同类型
adata_gex = adata[:, adata.var['feature_types'] == 'Gene Expression'].copy()
adata_adt = adata[:, adata.var['feature_types'] == 'Antibody Capture'].copy()
```

### Q7: 读取速度太慢怎么办？

**使用cache参数：**
```python
# 第一次读取（慢）
adata = sc.read_10x_mtx(
    'filtered_feature_bc_matrix/',
    cache=True  # 创建h5ad缓存
)

# 后续读取（快10-100倍）
# Scanpy会自动从缓存读取
adata = sc.read_10x_mtx(
    'filtered_feature_bc_matrix/',
    cache=True
)

# 缓存文件位置
# 默认在 ~/.cache/scanpy/ 目录下
```

### Q8: 如何读取人鼠混合样本数据？

**指定genome参数：**
```python
# 只读取人类基因组
adata_human = sc.read_10x_h5(
    'filtered_feature_bc_matrix.h5',
    genome='GRCh38'
)

# 只读取小鼠基因组
adata_mouse = sc.read_10x_h5(
    'filtered_feature_bc_matrix.h5',
    genome='mm10'
)

# 或者读取后分离
adata_all = sc.read_10x_h5('filtered_feature_bc_matrix.h5', gex_only=False)
adata_human = adata_all[:, adata_all.var['genome'] == 'GRCh38'].copy()
adata_mouse = adata_all[:, adata_all.var['genome'] == 'mm10'].copy()
```

---

## 8. 多样本数据读取与合并

在实际研究中，通常需要处理多个样本的单细胞数据。本节介绍多样本数据的读取和合并方法。

### 8.1 多样本读取策略

#### 8.1.1 场景分类

| 场景 | 数据来源 | 合并策略 | 批次效应 |
|------|---------|---------|---------|
| **技术重复** | 同一样本多次测序 | 直接合并 | 需要校正 |
| **生物重复** | 不同个体相同条件 | 合并+标记 | 需要校正 |
| **多条件/多组织** | 不同条件/不同组织 | 合并+标记+校正 | 强烈需要校正 |
| **跨平台** | 不同测序平台 | 合并+标记+强校正 | 强烈需要校正 |

---

### 8.2 多样本10X格式数据读取与合并

#### 8.2.1 目录结构示例

```
data/
├── sample1/
│   └── filtered_feature_bc_matrix/
│       ├── barcodes.tsv.gz
│       ├── features.tsv.gz
│       └── matrix.mtx.gz
├── sample2/
│   └── filtered_feature_bc_matrix/
│       ├── barcodes.tsv.gz
│       ├── features.tsv.gz
│       └── matrix.mtx.gz
└── sample3/
    └── filtered_feature_bc_matrix/
        ├── barcodes.tsv.gz
        ├── features.tsv.gz
        └── matrix.mtx.gz
```

#### 8.2.2 方法1：循环读取后合并（推荐）

```python
import scanpy as sc
import pandas as pd

# 定义样本信息
samples = {
    'sample1': {
        'path': 'data/sample1/filtered_feature_bc_matrix/',
        'condition': 'control',
        'batch': 'batch1'
    },
    'sample2': {
        'path': 'data/sample2/filtered_feature_bc_matrix/',
        'condition': 'treatment',
        'batch': 'batch1'
    },
    'sample3': {
        'path': 'data/sample3/filtered_feature_bc_matrix/',
        'condition': 'treatment',
        'batch': 'batch2'
    }
}

# 读取所有样本
adatas = []
for sample_id, sample_info in samples.items():
    # 读取数据
    adata = sc.read_10x_mtx(
        sample_info['path'],
        var_names='gene_symbols',
        cache=True
    )
    
    # 添加样本信息
    adata.obs['sample_id'] = sample_id
    adata.obs['condition'] = sample_info['condition']
    adata.obs['batch'] = sample_info['batch']
    
    # 添加样本前缀到细胞条形码（避免重复）
    adata.obs_names = [f"{sample_id}_{barcode}" for barcode in adata.obs_names]
    
    # 确保基因名唯一
    adata.var_names_make_unique()
    
    adatas.append(adata)
    print(f"✓ 读取 {sample_id}: {adata.n_obs} cells, {adata.n_vars} genes")

# 合并所有样本
adata_combined = sc.concat(
    adatas,
    axis=0,                    # 按细胞维度合并（行合并）
    join='outer',              # 'outer'保留所有基因，'inner'只保留共有基因
    label='sample_batch',      # 添加新的标签列
    keys=[s for s in samples.keys()],  # 标签值
    merge='unique',            # 处理重复的obs/var列
    uns_merge='unique'         # 处理uns字段
)

print(f"\n合并后: {adata_combined.n_obs} cells, {adata_combined.n_vars} genes")
print(f"样本分布:\n{adata_combined.obs['sample_id'].value_counts()}")

# 保存合并结果
adata_combined.write('combined_raw.h5ad')
```

#### 8.2.3 方法2：直接使用样本列表

```python
import scanpy as sc
import os

# 样本目录列表
sample_dirs = [
    'data/sample1/filtered_feature_bc_matrix/',
    'data/sample2/filtered_feature_bc_matrix/',
    'data/sample3/filtered_feature_bc_matrix/'
]

sample_names = ['sample1', 'sample2', 'sample3']

# 批量读取
adatas = []
for sample_dir, sample_name in zip(sample_dirs, sample_names):
    adata = sc.read_10x_mtx(sample_dir, var_names='gene_symbols')
    adata.obs['sample'] = sample_name
    adata.obs_names = [f"{sample_name}_{bc}" for bc in adata.obs_names]
    adatas.append(adata)

# 合并
adata_combined = sc.concat(adatas, join='outer')
```

#### 8.2.4 方法3：自动搜索目录

```python
import scanpy as sc
import os
import glob

def find_and_merge_10x_samples(root_dir, pattern='filtered_feature_bc_matrix'):
    """
    自动搜索并合并所有10X格式样本
    
    参数:
        root_dir: 根目录
        pattern: 搜索模式
    """
    # 搜索所有符合条件的目录
    sample_dirs = glob.glob(f'{root_dir}/**/{pattern}', recursive=True)
    
    if not sample_dirs:
        raise ValueError(f"在 {root_dir} 中未找到符合 {pattern} 的目录")
    
    print(f"找到 {len(sample_dirs)} 个样本:")
    
    adatas = []
    for sample_dir in sample_dirs:
        # 从路径提取样本名
        sample_name = os.path.basename(os.path.dirname(sample_dir))
        
        # 读取数据
        adata = sc.read_10x_mtx(sample_dir, var_names='gene_symbols')
        adata.obs['sample'] = sample_name
        adata.obs_names = [f"{sample_name}_{bc}" for bc in adata.obs_names]
        
        adatas.append(adata)
        print(f"  ✓ {sample_name}: {adata.n_obs} cells")
    
    # 合并
    adata_combined = sc.concat(adatas, join='outer')
    print(f"\n合并完成: {adata_combined.n_obs} cells, {adata_combined.n_vars} genes")
    
    return adata_combined

# 使用
adata_combined = find_and_merge_10x_samples('data/')
```

---

### 8.3 多样本H5AD格式数据读取与合并

#### 8.3.1 从已处理的H5AD文件合并

```python
import scanpy as sc

# H5AD文件列表
h5ad_files = {
    'sample1': 'data/sample1_processed.h5ad',
    'sample2': 'data/sample2_processed.h5ad',
    'sample3': 'data/sample3_processed.h5ad'
}

# 读取所有样本
adatas = []
for sample_id, file_path in h5ad_files.items():
    adata = sc.read_h5ad(file_path)
    
    # 添加样本标记
    adata.obs['sample_id'] = sample_id
    
    # 确保细胞条形码唯一
    adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]
    
    adatas.append(adata)
    print(f"✓ 读取 {sample_id}: {adata.n_obs} cells")

# 合并
adata_combined = sc.concat(adatas, join='outer', merge='unique')
print(f"\n合并后: {adata_combined.n_obs} cells, {adata_combined.n_vars} genes")
```

---

### 8.4 处理基因集不一致的情况

#### 8.4.1 Inner Join：只保留共有基因

```python
# 合并时只保留所有样本共有的基因
adata_combined = sc.concat(adatas, join='inner')

print(f"使用inner join后的基因数: {adata_combined.n_vars}")
```

#### 8.4.2 Outer Join：保留所有基因（推荐）

```python
# 合并时保留所有基因，缺失值用0填充
adata_combined = sc.concat(adatas, join='outer')

print(f"使用outer join后的基因数: {adata_combined.n_vars}")
```

#### 8.4.3 手动对齐基因集

```python
import numpy as np

def align_genes(adatas):
    """
    手动对齐多个样本的基因集
    """
    # 获取所有基因的并集
    all_genes = set()
    for adata in adatas:
        all_genes.update(adata.var_names)
    
    all_genes = sorted(list(all_genes))
    print(f"总共 {len(all_genes)} 个唯一基因")
    
    # 对每个样本进行对齐
    aligned_adatas = []
    for adata in adatas:
        # 创建新的表达矩阵
        new_X = np.zeros((adata.n_obs, len(all_genes)))
        
        # 填充已有基因的表达值
        gene_indices = {gene: i for i, gene in enumerate(all_genes)}
        for i, gene in enumerate(adata.var_names):
            if gene in gene_indices:
                new_X[:, gene_indices[gene]] = adata.X[:, i].toarray().flatten() if scipy.sparse.issparse(adata.X) else adata.X[:, i]
        
        # 创建新的AnnData对象
        import anndata as ad
        new_adata = ad.AnnData(
            X=new_X,
            obs=adata.obs.copy(),
            var=pd.DataFrame(index=all_genes)
        )
        
        aligned_adatas.append(new_adata)
    
    # 合并对齐后的数据
    adata_combined = sc.concat(aligned_adatas, join='outer')
    
    return adata_combined

# 使用
adata_combined = align_genes(adatas)
```

---

### 8.5 批次效应校正

合并多样本数据后，需要进行批次效应校正。

#### 8.5.1 使用Harmony（推荐）

```python
import scanpy as sc

# 读取合并后的数据
adata = sc.read_h5ad('combined_raw.h5ad')

# 标准化和归一化
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 识别高变基因
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='sample_id')

# PCA
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)

# 使用Harmony进行批次校正
import harmonypy
harmony_out = harmonypy.run_harmony(
    adata.obsm['X_pca'],
    adata.obs,
    'batch',              # 批次列名
    max_iter_harmony=20
)

# 保存校正后的PCA
adata.obsm['X_pca_harmony'] = harmony_out.Z_corr.T

# 使用校正后的PCA进行下游分析
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40, use_rep='X_pca_harmony')
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# 可视化
sc.pl.umap(adata, color=['batch', 'sample_id', 'leiden'], ncols=3)
```

#### 8.5.2 使用Scanpy内置的Combat

```python
import scanpy as sc

# 读取数据
adata = sc.read_h5ad('combined_raw.h5ad')

# 标准化
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Combat批次校正
sc.pp.combat(adata, key='batch')

# 后续分析
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
```

#### 8.5.3 使用BBKNN

```python
import scanpy as sc
import bbknn

# 读取数据
adata = sc.read_h5ad('combined_raw.h5ad')

# 预处理
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='batch')
sc.tl.pca(adata, n_comps=50)

# BBKNN批次校正
bbknn.bbknn(adata, batch_key='batch', n_pcs=40)

# 后续分析
sc.tl.umap(adata)
sc.tl.leiden(adata)

# 可视化
sc.pl.umap(adata, color=['batch', 'leiden'], ncols=2)
```

#### 8.5.4 使用Scanorama

```python
import scanpy as sc
import scanorama

# 读取数据
adata = sc.read_h5ad('combined_raw.h5ad')

# 预处理
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='batch')

# 按batch分割
adatas = [adata[adata.obs['batch'] == b].copy() for b in adata.obs['batch'].unique()]

# Scanorama整合
scanorama.integrate_scanpy(adatas, dimred=50)

# 合并整合后的数据
adata_integrated = sc.concat(adatas)

# 后续分析
sc.pp.neighbors(adata_integrated, use_rep='X_scanorama')
sc.tl.umap(adata_integrated)
sc.tl.leiden(adata_integrated)
```

---

### 8.6 多样本合并完整示例

```python
import scanpy as sc
import pandas as pd
import numpy as np
import os

def merge_multiple_samples(sample_info_dict, output_dir='merged_output'):
    """
    完整的多样本合并流程
    
    参数:
        sample_info_dict: 样本信息字典
            {
                'sample1': {
                    'path': '路径',
                    'condition': '条件',
                    'batch': '批次',
                    'patient_id': '病人ID' (可选)
                }
            }
        output_dir: 输出目录
    """
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 70)
    print("多样本数据合并流程")
    print("=" * 70)
    
    # 步骤1：读取所有样本
    print("\n步骤1：读取样本数据...")
    adatas = []
    for sample_id, info in sample_info_dict.items():
        # 读取数据
        adata = sc.read_10x_mtx(
            info['path'],
            var_names='gene_symbols',
            cache=True
        )
        
        # 添加元数据
        adata.obs['sample_id'] = sample_id
        for key, value in info.items():
            if key != 'path':
                adata.obs[key] = value
        
        # 确保条形码唯一
        adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]
        adata.var_names_make_unique()
        
        adatas.append(adata)
        print(f"  ✓ {sample_id}: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # 步骤2：合并样本
    print("\n步骤2：合并样本...")
    adata_combined = sc.concat(
        adatas,
        axis=0,
        join='outer',
        merge='unique'
    )
    print(f"  合并后: {adata_combined.n_obs} cells, {adata_combined.n_vars} genes")
    
    # 保存原始合并数据
    adata_combined.write(f'{output_dir}/01_raw_merged.h5ad')
    print(f"  ✓ 保存原始合并数据")
    
    # 步骤3：质量控制
    print("\n步骤3：质量控制...")
    # 计算QC指标
    sc.pp.calculate_qc_metrics(
        adata_combined,
        percent_top=None,
        log1p=False,
        inplace=True
    )
    
    # 过滤低质量细胞
    n_cells_before = adata_combined.n_obs
    sc.pp.filter_cells(adata_combined, min_genes=200)
    sc.pp.filter_genes(adata_combined, min_cells=3)
    n_cells_after = adata_combined.n_obs
    print(f"  过滤后: {n_cells_after} cells (过滤掉 {n_cells_before - n_cells_after} cells)")
    
    # 步骤4：标准化
    print("\n步骤4：标准化...")
    sc.pp.normalize_total(adata_combined, target_sum=1e4)
    sc.pp.log1p(adata_combined)
    
    # 步骤5：识别高变基因（考虑批次）
    print("\n步骤5：识别高变基因...")
    sc.pp.highly_variable_genes(
        adata_combined,
        n_top_genes=2000,
        batch_key='batch'  # 考虑批次效应
    )
    print(f"  识别到 {adata_combined.var['highly_variable'].sum()} 个高变基因")
    
    # 步骤6：PCA
    print("\n步骤6：PCA降维...")
    sc.tl.pca(adata_combined, n_comps=50, use_highly_variable=True)
    
    # 步骤7：批次校正（使用Harmony）
    print("\n步骤7：批次校正（Harmony）...")
    try:
        import harmonypy
        harmony_out = harmonypy.run_harmony(
            adata_combined.obsm['X_pca'],
            adata_combined.obs,
            'batch',
            max_iter_harmony=20
        )
        adata_combined.obsm['X_pca_harmony'] = harmony_out.Z_corr.T
        use_rep = 'X_pca_harmony'
        print("  ✓ Harmony批次校正完成")
    except ImportError:
        print("  ⚠️  Harmony未安装，跳过批次校正")
        use_rep = 'X_pca'
    
    # 步骤8：聚类和可视化
    print("\n步骤8：聚类和可视化...")
    sc.pp.neighbors(adata_combined, n_neighbors=15, n_pcs=40, use_rep=use_rep)
    sc.tl.umap(adata_combined)
    sc.tl.leiden(adata_combined, resolution=0.5)
    
    # 保存最终结果
    adata_combined.write(f'{output_dir}/02_merged_processed.h5ad')
    print(f"\n✓ 处理完成，结果保存在 {output_dir}/")
    
    # 生成可视化
    print("\n生成可视化...")
    import matplotlib.pyplot as plt
    
    # 按批次、样本、聚类可视化
    sc.pl.umap(
        adata_combined,
        color=['batch', 'sample_id', 'condition', 'leiden'],
        ncols=2,
        show=False
    )
    plt.savefig(f'{output_dir}/umap_overview.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 样本统计
    print("\n样本统计:")
    print(adata_combined.obs['sample_id'].value_counts())
    
    return adata_combined

# 使用示例
sample_info = {
    'sample1': {
        'path': 'data/sample1/filtered_feature_bc_matrix/',
        'condition': 'control',
        'batch': 'batch1',
        'patient_id': 'patient1'
    },
    'sample2': {
        'path': 'data/sample2/filtered_feature_bc_matrix/',
        'condition': 'treatment',
        'batch': 'batch1',
        'patient_id': 'patient2'
    },
    'sample3': {
        'path': 'data/sample3/filtered_feature_bc_matrix/',
        'condition': 'treatment',
        'batch': 'batch2',
        'patient_id': 'patient3'
    }
}

# 执行合并
adata_merged = merge_multiple_samples(sample_info, output_dir='merged_output')
```

---

### 8.7 使用Seurat合并多样本（R语言）

```r
library(Seurat)
library(harmony)

# 读取多个10X样本
sample_dirs <- c(
  "data/sample1/filtered_feature_bc_matrix/",
  "data/sample2/filtered_feature_bc_matrix/",
  "data/sample3/filtered_feature_bc_matrix/"
)

sample_names <- c("sample1", "sample2", "sample3")

# 读取并创建Seurat对象
seurat_list <- list()
for (i in seq_along(sample_dirs)) {
  # 读取数据
  data <- Read10X(data.dir = sample_dirs[i])
  
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(
    counts = data,
    project = sample_names[i],
    min.cells = 3,
    min.features = 200
  )
  
  # 添加样本信息
  seurat_obj$sample <- sample_names[i]
  
  seurat_list[[sample_names[i]]] <- seurat_obj
}

# 合并所有样本
seurat_merged <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.ids = sample_names,
  project = "merged_project"
)

# 标准化和PCA
seurat_merged <- NormalizeData(seurat_merged)
seurat_merged <- FindVariableFeatures(seurat_merged, nfeatures = 2000)
seurat_merged <- ScaleData(seurat_merged)
seurat_merged <- RunPCA(seurat_merged, npcs = 50)

# 使用Harmony批次校正
seurat_merged <- RunHarmony(
  seurat_merged,
  group.by.vars = "sample",
  reduction = "pca",
  assay.use = "RNA"
)

# 聚类和可视化
seurat_merged <- RunUMAP(seurat_merged, reduction = "harmony", dims = 1:40)
seurat_merged <- FindNeighbors(seurat_merged, reduction = "harmony", dims = 1:40)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.5)

# 可视化
DimPlot(seurat_merged, reduction = "umap", group.by = c("sample", "seurat_clusters"))
```

---

### 8.8 批次效应评估

合并和校正后，需要评估批次效应是否被有效去除。

#### 8.8.1 可视化评估

```python
import scanpy as sc
import matplotlib.pyplot as plt

# 读取数据
adata = sc.read_h5ad('merged_output/02_merged_processed.h5ad')

# 多角度可视化
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# 校正前后对比
sc.pl.pca(adata, color='batch', ax=axes[0, 0], show=False, title='PCA - 校正前')
sc.pl.pca(adata, color='batch', use_rep='X_pca_harmony', ax=axes[0, 1], show=False, title='PCA - 校正后')
sc.pl.umap(adata, color='batch', ax=axes[0, 2], show=False, title='UMAP - 校正后')

# 按样本、条件、聚类可视化
sc.pl.umap(adata, color='sample_id', ax=axes[1, 0], show=False)
sc.pl.umap(adata, color='condition', ax=axes[1, 1], show=False)
sc.pl.umap(adata, color='leiden', ax=axes[1, 2], show=False, legend_loc='on data')

plt.tight_layout()
plt.savefig('batch_correction_evaluation.png', dpi=300, bbox_inches='tight')
plt.close()
```

#### 8.8.2 定量评估

```python
import numpy as np

def calculate_batch_mixing_score(adata, batch_key='batch', cluster_key='leiden'):
    """
    计算批次混合分数
    分数越高表示批次效应去除越好
    """
    from sklearn.neighbors import NearestNeighbors
    
    # 使用UMAP坐标
    X = adata.obsm['X_umap']
    batches = adata.obs[batch_key].values
    
    # 对每个细胞找最近的k个邻居
    k = 50
    nbrs = NearestNeighbors(n_neighbors=k+1).fit(X)
    distances, indices = nbrs.kneighbors(X)
    
    # 计算每个细胞的邻居中其他批次的比例
    mixing_scores = []
    for i, batch in enumerate(batches):
        neighbor_batches = batches[indices[i, 1:]]  # 排除自己
        different_batch_ratio = np.mean(neighbor_batches != batch)
        mixing_scores.append(different_batch_ratio)
    
    return np.mean(mixing_scores)

# 计算混合分数
mixing_score = calculate_batch_mixing_score(adata)
print(f"批次混合分数: {mixing_score:.3f}")
print("(分数越高表示批次效应去除越好，建议 > 0.3)")
```

---

### 8.9 多样本合并最佳实践

#### 8.9.1 推荐流程

```
1. 读取样本 → 添加元数据 → 确保条形码唯一
   ↓
2. 合并样本（join='outer'）
   ↓
3. 质量控制（统一标准）
   ↓
4. 标准化和归一化
   ↓
5. 识别高变基因（batch-aware）
   ↓
6. PCA降维
   ↓
7. 批次校正（Harmony/Combat/BBKNN）
   ↓
8. 聚类和可视化
   ↓
9. 评估批次效应去除效果
```

#### 8.9.2 注意事项

1. **条形码唯一性**：
   - 不同样本可能有相同的条形码序列
   - 必须添加样本前缀：`{sample_id}_{barcode}`

2. **基因集处理**：
   - 推荐使用 `join='outer'` 保留所有基因
   - 缺失值会自动填充为0

3. **元数据管理**：
   - 添加完整的样本信息（sample_id, batch, condition等）
   - 便于后续分析和追溯

4. **批次校正选择**：
   - Harmony：最常用，效果好，速度快（推荐）
   - Combat：经典方法，但可能过度校正
   - BBKNN：保留局部结构
   - Scanorama：大规模数据整合

5. **质量控制**：
   - 对所有样本使用统一的QC标准
   - 避免样本间的系统性差异

---

## 9. 总结

单细胞RNA测序数据有多种格式，最常用的是：

1. **10X格式**：标准表达矩阵格式（三个文件）
2. **H5AD格式**：Scanpy分析的标准格式
3. **FASTQ格式**：原始测序数据（需要处理）

**推荐工作流程：**
- 从FASTQ → 使用Cell Ranger/dnbc4tools处理 → 10X格式
- 从10X格式 → 使用Scanpy读取 → H5AD格式
- 后续分析都在H5AD格式上进行

**多样本合并流程：**
- 循环读取 → 添加元数据 → 确保唯一性 → sc.concat合并
- 质量控制 → 标准化 → 批次校正（Harmony推荐）
- 聚类分析 → 评估批次效应

**关键工具：**
- **Scanpy**：Python单细胞分析（推荐）
- **Seurat**：R单细胞分析
- **Cell Ranger**：10X官方处理工具
- **dnbc4tools**：MGI DNB C4平台处理工具
- **Harmony**：批次效应校正（推荐）

