# 🚀 scRNA-DataHub 快速启动指南

## 欢迎使用 scRNA-DataHub！

这是一个3步快速启动指南，让您在5分钟内开始使用。

---

## 第1步：安装（1分钟）

### 方式A：使用Conda（推荐）⭐

```bash
# 进入项目目录
cd /Users/warm/华大智造/TCGA/gdc/scRNA-DataHub

# 创建独立的conda环境
conda env create -f environment.yml

# 激活环境
conda activate scrna-datahub

# 验证安装
python src/universal_reader.py --help
```

### 方式B：使用Python venv

```bash
# 进入项目目录
cd /Users/warm/华大智造/TCGA/gdc/scRNA-DataHub

# 创建虚拟环境
python -m venv venv

# 激活环境
source venv/bin/activate  # Linux/macOS

# 安装依赖
pip install -r requirements.txt

# 验证安装
python src/universal_reader.py --help
```

**验证成功标志：**

如果看到帮助信息，说明安装成功！✅

**重要提示：** 推荐使用conda环境，可以完全隔离依赖，避免与其他项目冲突。

---

## 第2步：测试（2分钟）

### 选项A：快速测试（推荐）

```bash
# 使用scanpy内置数据快速测试
bash scripts/quick_test.sh
```

**测试6种格式，生成6个H5AD文件** ✅

### 选项B：完整测试

```bash
# 创建所有格式的测试数据
python scripts/create_test_data.py

# 运行完整测试
cd test_data_all_formats
bash run_all_format_tests.sh
```

**测试12+种格式，生成12+个H5AD文件** ✅

---

## 第3步：使用（2分钟）

### 命令行使用

```bash
# 读取您的数据
python src/universal_reader.py \
  your_data/ \
  -o output.h5ad
```

### Python API使用

```python
from src.universal_reader import UniversalScRNAReader

# 创建读取器
reader = UniversalScRNAReader(verbose=True)

# 读取数据
adata = reader.read_auto('your_data/')

# 保存为H5AD
reader.save_h5ad(adata, 'output.h5ad')
```

---

## ✅ 完成！

现在您已经可以：

- ✅ 读取12+种单细胞数据格式
- ✅ 统一转换为H5AD格式
- ✅ 用于Scanpy/Seurat分析

---

## 📚 下一步

### 学习更多

- [完整使用教程](docs/user_guide.md) - 详细的使用说明
- [API文档](docs/api_reference.md) - 完整的API参考
- [数据格式详解](docs/data_formats.md) - 所有支持的格式
- [示例代码](examples/) - 更多实际案例

### 常用场景

1. **读取10X Genomics数据**
   ```bash
   python src/universal_reader.py \
     filtered_feature_bc_matrix/ \
     -o sample1.h5ad \
     --sample-id sample1
   ```

2. **读取DNB C4数据**
   ```bash
   python src/universal_reader.py \
     02.count/filter_matrix/ \
     -o sample1.h5ad \
     --sample-id sample1
   ```

3. **读取STARsolo输出**
   ```bash
   python src/universal_reader.py \
     Solo.out/Gene/filtered/ \
     -o sample1.h5ad \
     --no-compressed
   ```

4. **读取CSV矩阵**
   ```bash
   python src/universal_reader.py \
     matrix.csv \
     -o sample1.h5ad \
     --transpose
   ```

5. **批量处理**
   ```python
   python examples/batch_processing.py
   ```

---

## ❓ 遇到问题？

1. **查看文档**: [docs/faq.md](docs/faq.md)
2. **查看示例**: [examples/](examples/)
3. **提交Issue**: [GitHub Issues](https://github.com/yourusername/scRNA-DataHub/issues)

---

## 🌟 开始您的单细胞数据之旅！

```bash
# 一键测试
bash scripts/quick_test.sh

# 处理您的数据
python src/universal_reader.py your_data/ -o output.h5ad

# 继续分析
python -c "import scanpy as sc; adata = sc.read_h5ad('output.h5ad'); print(adata)"
```

**祝您分析顺利！** 🎉

