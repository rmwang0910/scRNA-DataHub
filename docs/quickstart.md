# 快速开始

5分钟快速上手 scRNA-DataHub

---

## 第1步：安装（2分钟）

### 方式A：使用 Conda（推荐）⭐

```bash
# 进入项目目录
cd scRNA-DataHub

# 创建 conda 环境
conda env create -f environment.yml

# 激活环境
conda activate scrna-datahub

# 验证安装
python src/universal_reader.py --help
```

### 方式B：使用 Python venv

```bash
# 进入项目目录
cd scRNA-DataHub

# 创建虚拟环境
python -m venv venv

# 激活环境
source venv/bin/activate  # Linux/macOS

# 安装依赖
pip install -r requirements.txt

# 验证安装
python src/universal_reader.py --help
```

---

## 第2步：测试（2分钟）

### 快速测试（推荐）

```bash
# 进入测试目录
cd scripts/test_data_all_formats

# 运行快速测试（测试 6 种格式）
bash QUICK_START.sh
```

### 完整测试

```bash
# 进入测试目录
cd scripts/test_data_all_formats

# 运行完整测试（测试 17 种格式）
bash run_all_format_tests_simple.sh
```

---

## 第3步：使用（1分钟）

### 命令行使用

```bash
# 读取 10X Genomics 数据
python src/universal_reader.py \
  filtered_feature_bc_matrix/ \
  -o output.h5ad

# 读取 DNB C4 数据
python src/universal_reader.py \
  02.count/filter_matrix/ \
  -o output.h5ad

# 读取 CSV 矩阵
python src/universal_reader.py \
  matrix.csv \
  -o output.h5ad \
  --transpose
```

### Python API 使用

```python
from src.universal_reader import UniversalScRNAReader

# 创建读取器
reader = UniversalScRNAReader(verbose=True)

# 读取数据（自动检测格式）
adata = reader.read_auto('filtered_feature_bc_matrix/')

# 保存为 H5AD
reader.save_h5ad(adata, 'output.h5ad')
```

### 验证结果

```python
import scanpy as sc

# 读取转换后的数据
adata = sc.read_h5ad('output.h5ad')

# 查看数据信息
print(adata)
# AnnData object with n_obs × n_vars = 2700 × 32738
#     obs: 'sample_id'
#     var: 'gene_ids', 'feature_types'
```

---

## ✅ 完成！

现在你已经可以：

- ✅ 读取 17 种单细胞数据格式
- ✅ 统一转换为 H5AD 格式
- ✅ 用于 Scanpy/Seurat 分析

---

## 常用场景

### 场景1：10X Genomics 数据

```bash
# Cell Ranger v3+ 输出
python src/universal_reader.py \
  sample1/outs/filtered_feature_bc_matrix/ \
  -o sample1.h5ad \
  --sample-id sample1
```

### 场景2：DNB C4 数据

```bash
# dnbc4tools 输出
python src/universal_reader.py \
  CNS1063416_brain/02.count/filter_matrix/ \
  -o CNS1063416_brain.h5ad \
  --sample-id CNS1063416_brain
```

### 场景3：STARsolo 输出

```bash
# STARsolo 输出（未压缩）
python src/universal_reader.py \
  Solo.out/Gene/filtered/ \
  -o sample1.h5ad \
  --no-compressed
```

### 场景4：批量处理

```python
# 查看示例代码
python examples/batch_processing.py
```

---

## 📚 下一步

### 深入学习

- [完整使用教程](user_guide.md) - 详细的使用说明
- [API 文档](api_reference.md) - 完整的 API 参考
- [数据格式详解](data_formats.md) - 所有支持的格式
- [示例代码](../examples/) - 更多实际案例

### 遇到问题？

- [常见问题](faq.md) - FAQ
- [故障排除](troubleshooting.md) - 常见错误解决方案
- [GitHub Issues](https://github.com/yourusername/scRNA-DataHub/issues) - 提交问题

---

## 🎉 开始分析！

```bash
# 一键转换
python src/universal_reader.py your_data/ -o output.h5ad

# 开始 Scanpy 分析
python -c "import scanpy as sc; adata = sc.read_h5ad('output.h5ad'); print(adata)"
```

**祝您分析顺利！** 🚀
