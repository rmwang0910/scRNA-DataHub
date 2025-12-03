# 通用单细胞数据读取器使用指南

## 概述

`universal_scrna_reader.py` 是一个通用的单细胞RNA-seq数据读取工具，支持所有常见的单细胞数据格式，并统一转换为标准的H5AD格式。

## 支持的格式

### 完整格式列表

| 格式 | 文件/目录 | 自动检测 | 特殊参数 |
|------|----------|---------|---------|
| **10X MTX (v3+)** | `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz` | ✅ | `--no-compressed` (STARsolo) |
| **10X MTX (v2-)** | `matrix.mtx`, `barcodes.tsv`, `genes.tsv` | ✅ | 自动识别 |
| **10X H5** | `*.h5` | ✅ | `--genome` (多基因组) |
| **H5AD** | `*.h5ad` | ✅ | `--backed` (大数据) |
| **Loom** | `*.loom` | ✅ | - |
| **Zarr** | `*.zarr/` | ✅ | - |
| **CSV** | `*.csv` | ✅ | `--transpose` |
| **TSV/TXT** | `*.tsv`, `*.txt` | ✅ | `--delimiter`, `--transpose` |
| **Excel** | `*.xlsx`, `*.xls` | ✅ | `--sheet` |
| **MTX** | `*.mtx` | ✅ | - |
| **HDF5** | `*.h5` | ✅ | - |
| **SOFT.GZ** | `*.soft.gz` | ✅ | - |
| **UMI-tools** | `*.tsv.gz` | ❌ | 手动指定 |

---

## 安装

### 依赖包

```bash
# 基础依赖
pip install scanpy anndata pandas numpy scipy

# 可选依赖（用于批次校正）
pip install harmonypy bbknn scanorama

# Zarr支持
pip install zarr

# Loom支持
pip install loompy
```

---

## 快速开始

### 1. 基本用法

```bash
# 自动检测格式并读取
python universal_scrna_reader.py input_data/ -o output.h5ad
```

### 2. 作为Python模块使用

```python
from universal_scrna_reader import UniversalScRNAReader

# 创建读取器
reader = UniversalScRNAReader(verbose=True)

# 自动读取任意格式
adata = reader.read_auto('filtered_feature_bc_matrix/')

# 保存为H5AD
reader.save_h5ad(adata, 'output.h5ad')
```

---

## 详细用法

### 1. 读取10X Genomics数据

#### 1.1 Cell Ranger v3+ 输出（压缩）

```bash
python universal_scrna_reader.py \
    filtered_feature_bc_matrix/ \
    -o sample1.h5ad \
    --sample-id sample1
```

#### 1.2 STARsolo输出（未压缩）

```bash
python universal_scrna_reader.py \
    Solo.out/Gene/filtered/ \
    -o sample1.h5ad \
    --no-compressed \
    --sample-id sample1
```

#### 1.3 10X H5格式

```bash
python universal_scrna_reader.py \
    filtered_feature_bc_matrix.h5 \
    -o sample1.h5ad \
    --sample-id sample1
```

#### 1.4 人鼠混合样本

```bash
# 只读取人类基因组
python universal_scrna_reader.py \
    filtered_feature_bc_matrix.h5 \
    -o human_cells.h5ad \
    --genome GRCh38

# 只读取小鼠基因组
python universal_scrna_reader.py \
    filtered_feature_bc_matrix.h5 \
    -o mouse_cells.h5ad \
    --genome mm10
```

---

### 2. 读取文本格式数据

#### 2.1 CSV格式

```bash
# 行是基因，列是细胞（需要转置）
python universal_scrna_reader.py \
    expression_matrix.csv \
    -o output.h5ad \
    --transpose

# 行是细胞，列是基因（不需要转置）
python universal_scrna_reader.py \
    expression_matrix.csv \
    -o output.h5ad
```

#### 2.2 TSV格式

```bash
python universal_scrna_reader.py \
    expression_matrix.tsv \
    -o output.h5ad \
    --delimiter '\t' \
    --transpose
```

#### 2.3 自定义分隔符

```bash
python universal_scrna_reader.py \
    data.txt \
    -o output.h5ad \
    --delimiter '|'
```

---

### 3. 读取其他格式

#### 3.1 Loom格式

```bash
python universal_scrna_reader.py \
    data.loom \
    -o output.h5ad
```

#### 3.2 Zarr格式

```bash
python universal_scrna_reader.py \
    data.zarr/ \
    -o output.h5ad
```

#### 3.3 Excel格式

```bash
python universal_scrna_reader.py \
    data.xlsx \
    -o output.h5ad \
    --sheet Sheet1
```

---

## Python API使用

### 1. 基础使用

```python
from universal_scrna_reader import UniversalScRNAReader

# 创建读取器
reader = UniversalScRNAReader(verbose=True)

# 自动检测并读取
adata = reader.read_auto('filtered_feature_bc_matrix/')

# 验证数据
stats = reader.validate_adata(adata)

# 标准化
adata = reader.standardize_adata(adata)

# 保存
reader.save_h5ad(adata, 'output.h5ad')
```

### 2. 读取特定格式

```python
from universal_scrna_reader import UniversalScRNAReader

reader = UniversalScRNAReader(verbose=True)

# 10X MTX格式
adata = reader.read_10x_mtx(
    'filtered_feature_bc_matrix/',
    var_names='gene_symbols',
    cache=True,
    compressed=True
)

# 10X H5格式
adata = reader.read_10x_h5(
    'filtered_feature_bc_matrix.h5',
    genome='GRCh38',
    gex_only=True
)

# H5AD格式（大数据用backed模式）
adata = reader.read_h5ad(
    'large_data.h5ad',
    backed='r'  # 只读，不全部加载到内存
)

# CSV格式
adata = reader.read_csv(
    'expression.csv',
    first_column_names=True,
    transpose=True
)

# TSV格式
adata = reader.read_text(
    'expression.tsv',
    delimiter='\t',
    transpose=True
)

# Loom格式
adata = reader.read_loom('data.loom')

# Zarr格式
adata = reader.read_zarr('data.zarr')
```

### 3. 完整处理流程

```python
from universal_scrna_reader import UniversalScRNAReader

# 创建读取器
reader = UniversalScRNAReader(verbose=True)

# 一键处理：读取、验证、标准化、保存
adata = reader.process_single_sample(
    input_path='filtered_feature_bc_matrix/',
    output_path='sample1.h5ad',
    sample_id='sample1',
    var_names='gene_symbols',
    cache=True
)

# 此时adata已经是标准化的AnnData对象
print(adata)
```

---

## 测试脚本

### 运行测试

```bash
# 运行完整测试（创建测试数据并测试所有格式）
python test_universal_reader.py
```

**测试内容：**
1. 创建测试数据（100 cells × 500 genes）
2. 保存为7种不同格式
3. 分别读取并验证
4. 统一转换为H5AD格式
5. 验证一致性

**测试输出：**
```
创建测试数据: 100 cells × 500 genes
测试1: 10X MTX格式 (压缩)
  ✓ 读取成功
测试2: 10X MTX格式 (未压缩)
  ✓ 读取成功
测试3: H5AD格式
  ✓ 读取成功
...
结果验证：
  ✓ 所有格式的维度一致: 100 cells × 500 genes
  ✓ 所有读取结果都是AnnData对象
  ✓ 所有格式都可以统一转换为H5AD格式
```

---

## 常用示例

### 示例1：读取Cell Ranger输出

```python
from universal_scrna_reader import UniversalScRNAReader

reader = UniversalScRNAReader()

# Cell Ranger v3+ 标准输出
adata = reader.process_single_sample(
    input_path='sample1/outs/filtered_feature_bc_matrix/',
    output_path='sample1.h5ad',
    sample_id='sample1'
)
```

### 示例2：读取dnbc4tools输出

```python
from universal_scrna_reader import UniversalScRNAReader

reader = UniversalScRNAReader()

# dnbc4tools输出（10X格式）
adata = reader.process_single_sample(
    input_path='CNS1063416_brain/02.count/filter_matrix/',
    output_path='CNS1063416_brain.h5ad',
    sample_id='CNS1063416_brain'
)
```

### 示例3：读取STARsolo输出

```python
from universal_scrna_reader import UniversalScRNAReader

reader = UniversalScRNAReader()

# STARsolo输出（未压缩）
adata = reader.process_single_sample(
    input_path='Solo.out/Gene/filtered/',
    output_path='sample1.h5ad',
    sample_id='sample1',
    compressed=False  # 关键参数
)
```

### 示例4：读取CSV矩阵

```python
from universal_scrna_reader import UniversalScRNAReader

reader = UniversalScRNAReader()

# CSV文件（行是基因，列是细胞）
adata = reader.process_single_sample(
    input_path='expression_matrix.csv',
    output_path='sample1.h5ad',
    sample_id='sample1',
    transpose=True  # 转置为 cells × genes
)
```

### 示例5：处理大数据（backed模式）

```python
from universal_scrna_reader import UniversalScRNAReader

reader = UniversalScRNAReader()

# 使用backed模式（不全部加载到内存）
adata = reader.read_h5ad(
    'large_data.h5ad',
    backed='r'
)

# 只处理部分数据
subset = adata[:10000, :].to_memory()

# 或者分批处理
for i in range(0, adata.n_obs, 10000):
    batch = adata[i:i+10000, :].to_memory()
    # 处理batch...
```

---

## 高级功能

### 1. 批量处理多个文件

```python
from universal_scrna_reader import UniversalScRNAReader
from pathlib import Path

reader = UniversalScRNAReader(verbose=True)

# 批量处理
input_files = [
    'sample1/filtered_feature_bc_matrix/',
    'sample2/filtered_feature_bc_matrix/',
    'sample3/filtered_feature_bc_matrix/'
]

sample_ids = ['sample1', 'sample2', 'sample3']

for input_file, sample_id in zip(input_files, sample_ids):
    output_file = f'{sample_id}.h5ad'
    
    adata = reader.process_single_sample(
        input_path=input_file,
        output_path=output_file,
        sample_id=sample_id
    )
    
    print(f"✓ {sample_id} 处理完成\n")
```

### 2. 格式检测

```python
from universal_scrna_reader import UniversalScRNAReader

reader = UniversalScRNAReader(verbose=False)

# 检测格式
format_type = reader.detect_format('filtered_feature_bc_matrix/')
print(f"检测到的格式: {format_type}")

# 支持的格式列表
print("支持的格式:")
for fmt, desc in reader.supported_formats.items():
    print(f"  - {fmt}: {desc}")
```

### 3. 数据验证

```python
from universal_scrna_reader import UniversalScRNAReader
import scanpy as sc

reader = UniversalScRNAReader(verbose=True)

# 读取数据
adata = sc.read_h5ad('data.h5ad')

# 验证和统计
stats = reader.validate_adata(adata)

# 查看统计信息
print(f"细胞数: {stats['n_obs']}")
print(f"基因数: {stats['n_vars']}")
print(f"是否稀疏: {stats['is_sparse']}")
print(f"元数据列: {stats['obs_keys']}")
```

---

## 命令行参数

### 完整参数列表

```bash
python universal_scrna_reader.py --help
```

**主要参数：**

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `input` | 输入文件或目录路径 | 必需 |
| `-o, --output` | 输出H5AD文件路径 | None |
| `--sample-id` | 样本ID | None |
| `--format` | 强制指定格式 | 自动检测 |
| `--var-names` | 10X: gene_symbols 或 gene_ids | gene_symbols |
| `--genome` | 10X H5: 基因组名称 | None |
| `--no-gex-only` | 10X: 保留所有特征类型 | False |
| `--no-compressed` | 10X MTX: 未压缩（STARsolo） | False |
| `--delimiter` | 文本: 分隔符 | '\t' |
| `--transpose` | 文本: 转置矩阵 | False |
| `--sheet` | Excel: sheet名称 | Sheet1 |
| `--backed` | H5AD: backed模式 | None |
| `--no-cache` | 10X: 不使用缓存 | False |
| `--quiet` | 静默模式 | False |

---

## 实际应用场景

### 场景1：处理Cell Ranger输出

```bash
# 标准Cell Ranger v3+输出
python universal_scrna_reader.py \
    /data/cellranger_output/sample1/outs/filtered_feature_bc_matrix/ \
    -o /output/sample1.h5ad \
    --sample-id sample1 \
    --var-names gene_symbols
```

### 场景2：处理dnbc4tools输出

```bash
# DNB C4平台输出
python universal_scrna_reader.py \
    /data/dnbc4tools_output/CNS1063416_brain/02.count/filter_matrix/ \
    -o /output/CNS1063416_brain.h5ad \
    --sample-id CNS1063416_brain
```

### 场景3：处理STARsolo输出

```bash
# STARsolo输出（注意：未压缩）
python universal_scrna_reader.py \
    /data/star_output/Solo.out/Gene/filtered/ \
    -o /output/sample1.h5ad \
    --sample-id sample1 \
    --no-compressed
```

### 场景4：处理GEO数据

```bash
# 从GEO下载的SOFT格式
python universal_scrna_reader.py \
    GSE12345_family.soft.gz \
    -o GSE12345.h5ad
```

### 场景5：处理自定义CSV矩阵

```bash
# 行是基因，列是细胞
python universal_scrna_reader.py \
    expression_matrix.csv \
    -o sample1.h5ad \
    --transpose \
    --sample-id sample1
```

---

## 完整Python脚本示例

### 示例1：批量处理多个10X样本

```python
#!/usr/bin/env python3
from universal_scrna_reader import UniversalScRNAReader
from pathlib import Path

# 创建读取器
reader = UniversalScRNAReader(verbose=True)

# 定义样本信息
samples = [
    {
        'input': '/data/sample1/outs/filtered_feature_bc_matrix/',
        'output': '/output/sample1.h5ad',
        'sample_id': 'sample1'
    },
    {
        'input': '/data/sample2/outs/filtered_feature_bc_matrix/',
        'output': '/output/sample2.h5ad',
        'sample_id': 'sample2'
    },
    {
        'input': '/data/sample3/outs/filtered_feature_bc_matrix/',
        'output': '/output/sample3.h5ad',
        'sample_id': 'sample3'
    }
]

# 批量处理
for sample_info in samples:
    print(f"\n处理 {sample_info['sample_id']}...")
    
    adata = reader.process_single_sample(
        input_path=sample_info['input'],
        output_path=sample_info['output'],
        sample_id=sample_info['sample_id']
    )
    
    print(f"✓ {sample_info['sample_id']} 完成")

print("\n所有样本处理完成！")
```

### 示例2：混合格式批量处理

```python
#!/usr/bin/env python3
from universal_scrna_reader import UniversalScRNAReader

reader = UniversalScRNAReader(verbose=True)

# 不同格式的样本
samples = [
    {
        'input': 'sample1/filtered_feature_bc_matrix/',
        'format': '10x_mtx',
        'kwargs': {'compressed': True}
    },
    {
        'input': 'sample2/filtered_feature_bc_matrix.h5',
        'format': '10x_h5',
        'kwargs': {}
    },
    {
        'input': 'sample3/expression.csv',
        'format': 'csv',
        'kwargs': {'transpose': True}
    }
]

adatas = []
for i, sample in enumerate(samples):
    sample_id = f"sample{i+1}"
    
    # 读取
    adata = reader.read_auto(
        sample['input'],
        **sample['kwargs']
    )
    
    # 添加样本信息
    adata.obs['sample_id'] = sample_id
    
    # 标准化
    adata = reader.standardize_adata(adata)
    
    # 保存
    reader.save_h5ad(adata, f'{sample_id}.h5ad')
    
    adatas.append(adata)

print(f"\n✓ 处理了 {len(adatas)} 个样本")
```

---

## 输出格式说明

所有输入格式都会统一转换为标准的**AnnData/H5AD格式**，包含以下结构：

```python
AnnData object with n_obs × n_vars
    obs: 细胞元数据
        - sample_id: 样本ID（如果提供）
        - （其他可能的列）
    var: 基因元数据
        - gene_ids: 基因ID（如果有）
        - feature_types: 特征类型（如果有）
        - （其他可能的列）
    uns: 非结构化数据
    obsm: 多维数组（如PCA、UMAP等）
    layers: 其他数据层
```

**统一特点：**
- ✅ 所有数据都是 `AnnData` 对象
- ✅ 矩阵维度统一为 `(n_cells, n_genes)`
- ✅ 使用稀疏矩阵存储（节省内存）
- ✅ 基因名和细胞条形码唯一
- ✅ 可以直接用于Scanpy分析

---

## 故障排除

### Q1: 读取10X数据时报错 "index 0 is out of bounds"

**原因：** 可能是压缩设置不对

**解决：**
```bash
# 如果是STARsolo输出，添加 --no-compressed
python universal_scrna_reader.py \
    Solo.out/Gene/filtered/ \
    -o output.h5ad \
    --no-compressed
```

### Q2: 内存不足

**解决方案1：使用backed模式**
```python
adata = reader.read_h5ad('large_data.h5ad', backed='r')
```

**解决方案2：使用Zarr格式**
```bash
# 先转换为Zarr
python universal_scrna_reader.py input.h5ad -o output.zarr
```

### Q3: 基因名重复

**自动处理：**
```python
# 读取器会自动处理重复基因名（添加后缀-1, -2等）
adata = reader.standardize_adata(adata, ensure_unique_names=True)
```

### Q4: CSV/TSV格式维度不对

**检查是否需要转置：**
```bash
# 如果行是基因，列是细胞，需要转置
python universal_scrna_reader.py data.csv -o output.h5ad --transpose
```

---

## 性能优化

### 1. 使用缓存加速读取

```python
# 第一次读取时创建缓存
adata = reader.read_10x_mtx('filtered_feature_bc_matrix/', cache=True)

# 后续读取速度提升10-100倍
adata = reader.read_10x_mtx('filtered_feature_bc_matrix/', cache=True)
```

### 2. 使用backed模式处理大数据

```python
# 不全部加载到内存
adata = reader.read_h5ad('large_data.h5ad', backed='r')

# 只加载需要的部分
subset = adata[:10000, :].to_memory()
```

### 3. 使用稀疏矩阵

```python
# 读取器默认会转换为稀疏矩阵
adata = reader.standardize_adata(adata, ensure_sparse=True)
```

---

## 总结

### 关键特性

1. ✅ **自动格式检测**：无需手动指定格式
2. ✅ **统一输出**：所有格式统一转换为AnnData/H5AD
3. ✅ **完整验证**：自动验证数据完整性
4. ✅ **智能标准化**：稀疏矩阵、唯一基因名等
5. ✅ **灵活使用**：支持命令行和Python API
6. ✅ **性能优化**：缓存、backed模式、稀疏矩阵

### 推荐工作流程

```
任意格式数据
  ↓
[universal_scrna_reader.py]
  ↓
统一的H5AD格式
  ↓
[Scanpy/OmicVerse分析]
  ↓
分析结果
```

### 文件清单

- `universal_scrna_reader.py`：主脚本（读取器类）
- `test_universal_reader.py`：测试脚本
- `通用单细胞数据读取器使用指南.md`：本文档

开始使用：
```bash
python test_universal_reader.py  # 运行测试
python universal_scrna_reader.py your_data/ -o output.h5ad  # 实际使用
```

