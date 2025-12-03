# API参考文档

## UniversalScRNAReader类

### 初始化

```python
from src.universal_reader import UniversalScRNAReader

reader = UniversalScRNAReader(verbose=True)
```

**参数：**
- `verbose` (bool): 是否打印详细信息，默认True

---

## 核心方法

### read_auto()

自动检测格式并读取数据。

```python
adata = reader.read_auto(input_path, **kwargs)
```

**参数：**
- `input_path` (str): 输入文件或目录路径
- `**kwargs`: 传递给特定读取函数的参数

**返回：**
- `AnnData`: AnnData对象

**示例：**
```python
adata = reader.read_auto('filtered_feature_bc_matrix/')
```

---

### read_10x_mtx()

读取10X MTX格式数据。

```python
adata = reader.read_10x_mtx(
    path,
    var_names='gene_symbols',
    make_unique=True,
    cache=True,
    gex_only=True,
    compressed=True
)
```

**参数：**
- `path` (str): 包含matrix.mtx, barcodes.tsv, features.tsv的目录
- `var_names` (str): 'gene_symbols' 或 'gene_ids'
- `make_unique` (bool): 是否使基因名唯一
- `cache` (bool): 是否使用缓存
- `gex_only` (bool): 是否只保留Gene Expression
- `compressed` (bool): 文件是否压缩

**返回：**
- `AnnData`: AnnData对象

---

### read_10x_h5()

读取10X H5格式数据。

```python
adata = reader.read_10x_h5(path, genome=None, gex_only=True)
```

**参数：**
- `path` (str): H5文件路径
- `genome` (str, optional): 基因组名称
- `gex_only` (bool): 是否只保留Gene Expression

**返回：**
- `AnnData`: AnnData对象

---

### read_h5ad()

读取H5AD格式数据。

```python
adata = reader.read_h5ad(path, backed=None)
```

**参数：**
- `path` (str): H5AD文件路径
- `backed` (str, optional): 'r' 或 'r+' 启用backed模式

**返回：**
- `AnnData`: AnnData对象

---

### 其他读取方法

```python
# Loom格式
adata = reader.read_loom(path)

# CSV格式
adata = reader.read_csv(path, first_column_names=True, transpose=False)

# TSV/TXT格式
adata = reader.read_text(path, delimiter='\t', transpose=False)

# Excel格式
adata = reader.read_excel(path, sheet='Sheet1', transpose=False)

# Zarr格式
adata = reader.read_zarr(path)

# MTX单文件
adata = reader.read_mtx(path)

# HDF5格式
adata = reader.read_hdf5(path, key='data')

# SOFT.GZ格式
adata = reader.read_soft_gz(path)

# UMI-tools格式
adata = reader.read_umi_tools(path)
```

---

## 工具方法

### detect_format()

自动检测数据格式。

```python
format_type = reader.detect_format(input_path)
```

**参数：**
- `input_path` (str): 输入路径

**返回：**
- `str`: 格式类型

**示例：**
```python
format_type = reader.detect_format('filtered_feature_bc_matrix/')
print(f"检测到格式: {format_type}")  # 输出: 10x_mtx
```

---

### validate_adata()

验证和检查AnnData对象。

```python
stats = reader.validate_adata(adata)
```

**参数：**
- `adata` (AnnData): AnnData对象

**返回：**
- `dict`: 包含统计信息的字典

**返回字段：**
```python
{
    'n_obs': 细胞数,
    'n_vars': 基因数,
    'is_sparse': 是否稀疏矩阵,
    'has_raw': 是否有raw数据,
    'obs_keys': 细胞元数据列,
    'var_keys': 基因元数据列,
    'obsm_keys': 多维数组键,
    'uns_keys': 非结构化数据键,
    'layers_keys': 数据层键
}
```

---

### standardize_adata()

标准化AnnData对象。

```python
adata = reader.standardize_adata(
    adata,
    ensure_sparse=True,
    ensure_unique_names=True
)
```

**参数：**
- `adata` (AnnData): 输入的AnnData对象
- `ensure_sparse` (bool): 确保X是稀疏矩阵
- `ensure_unique_names` (bool): 确保基因名唯一

**返回：**
- `AnnData`: 标准化后的AnnData对象

---

### save_h5ad()

保存为H5AD格式。

```python
reader.save_h5ad(adata, output_path, compression='gzip')
```

**参数：**
- `adata` (AnnData): AnnData对象
- `output_path` (str): 输出文件路径
- `compression` (str): 压缩方式，'gzip', 'lzf', 或 None

---

### process_single_sample()

一键处理单个样本。

```python
adata = reader.process_single_sample(
    input_path,
    output_path=None,
    sample_id=None,
    **read_kwargs
)
```

**参数：**
- `input_path` (str): 输入路径
- `output_path` (str, optional): 输出H5AD路径
- `sample_id` (str, optional): 样本ID
- `**read_kwargs`: 传递给读取函数的参数

**返回：**
- `AnnData`: 标准化后的AnnData对象

**执行流程：**
1. 自动检测格式
2. 读取数据
3. 添加样本ID（如果提供）
4. 验证数据
5. 标准化
6. 保存（如果提供output_path）

---

## 属性

### supported_formats

```python
formats = reader.supported_formats
```

**返回：**
- `dict`: 支持的格式字典

```python
{
    '10x_mtx': '10X Genomics MTX格式 (3个文件)',
    '10x_h5': '10X Genomics H5格式',
    'h5ad': 'AnnData H5AD格式',
    'loom': 'Loom格式',
    ...
}
```

---

## 完整示例

### 示例1：基本使用

```python
from src.universal_reader import UniversalScRNAReader

# 创建读取器
reader = UniversalScRNAReader(verbose=True)

# 读取数据
adata = reader.read_auto('filtered_feature_bc_matrix/')

# 验证
stats = reader.validate_adata(adata)
print(f"细胞数: {stats['n_obs']}")
print(f"基因数: {stats['n_vars']}")

# 标准化
adata = reader.standardize_adata(adata)

# 保存
reader.save_h5ad(adata, 'output.h5ad')
```

### 示例2：批量处理

```python
from src.universal_reader import UniversalScRNAReader
from pathlib import Path

reader = UniversalScRNAReader()

# 批量处理
sample_dirs = Path('data/').glob('sample*/filtered_feature_bc_matrix/')

for sample_dir in sample_dirs:
    sample_id = sample_dir.parent.parent.name
    
    adata = reader.process_single_sample(
        input_path=str(sample_dir),
        output_path=f'{sample_id}.h5ad',
        sample_id=sample_id
    )
```

### 示例3：格式转换

```python
# CSV → H5AD
reader.process_single_sample(
    'matrix.csv',
    'output.h5ad',
    transpose=True
)

# Loom → H5AD
reader.process_single_sample(
    'data.loom',
    'output.h5ad'
)

# 10X MTX → H5AD
reader.process_single_sample(
    'filtered_feature_bc_matrix/',
    'output.h5ad',
    var_names='gene_symbols',
    cache=True
)
```

---

## 返回值说明

所有读取方法都返回标准的AnnData对象：

```python
AnnData object with n_obs × n_vars
    X: 表达矩阵 (稀疏矩阵)
    obs: 细胞元数据
        - sample_id: 样本ID（如果提供）
    var: 基因元数据
        - gene_ids: 基因ID
        - feature_types: 特征类型
    uns: 非结构化数据
    obsm: 多维数组
    layers: 其他数据层
```

---

## 错误处理

所有方法都会抛出标准Python异常：

- `ValueError`: 参数错误或格式不支持
- `FileNotFoundError`: 文件不存在
- `IOError`: 读写错误
- `KeyError`: 缺少必需字段

---

## 更多信息

- [完整使用教程](user_guide.md)
- [数据格式详解](data_formats.md)
- [示例代码](../examples/)

