# 快速开始

## 5分钟入门

### 1. 安装

```bash
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub
pip install -r requirements.txt
```

### 2. 读取数据

```bash
# 读取10X Genomics数据
python src/universal_reader.py \
  filtered_feature_bc_matrix/ \
  -o output.h5ad
```

### 3. 查看结果

```python
import scanpy as sc

# 读取转换后的数据
adata = sc.read_h5ad('output.h5ad')

# 查看数据
print(adata)
# AnnData object with n_obs × n_vars = 2700 × 32738
#     obs: 'sample_id'
#     var: 'gene_ids', 'feature_types'
```

## 常用场景

### 场景1：处理10X Genomics数据

```bash
# Cell Ranger v3+ 输出
python src/universal_reader.py \
  sample1/outs/filtered_feature_bc_matrix/ \
  -o sample1.h5ad \
  --sample-id sample1
```

### 场景2：处理DNB C4数据

```bash
# dnbc4tools输出
python src/universal_reader.py \
  CNS1063416_brain/02.count/filter_matrix/ \
  -o CNS1063416_brain.h5ad \
  --sample-id CNS1063416_brain
```

### 场景3：处理STARsolo输出

```bash
# STARsolo输出（注意：未压缩）
python src/universal_reader.py \
  Solo.out/Gene/filtered/ \
  -o sample1.h5ad \
  --no-compressed
```

### 场景4：处理CSV矩阵

```bash
# CSV文件（行是基因，列是细胞）
python src/universal_reader.py \
  expression_matrix.csv \
  -o sample1.h5ad \
  --transpose
```

## Python API使用

```python
from src.universal_reader import UniversalScRNAReader

# 创建读取器
reader = UniversalScRNAReader(verbose=True)

# 方法1：一键处理
adata = reader.process_single_sample(
    input_path='filtered_feature_bc_matrix/',
    output_path='sample1.h5ad',
    sample_id='sample1'
)

# 方法2：分步处理
adata = reader.read_auto('filtered_feature_bc_matrix/')
stats = reader.validate_adata(adata)
adata = reader.standardize_adata(adata)
reader.save_h5ad(adata, 'sample1.h5ad')
```

## 测试工具

```bash
# 快速测试（使用scanpy内置数据）
bash scripts/quick_test.sh

# 完整测试（创建所有格式）
python scripts/create_test_data.py
cd test_data_all_formats
bash run_all_format_tests.sh
```

## 下一步

- 查看[完整使用教程](user_guide.md)
- 了解[所有支持的格式](data_formats.md)
- 查看[更多示例](../examples/)
- 遇到问题查看[常见问题](faq.md)

