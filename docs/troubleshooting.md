# 故障排除指南

## 常见错误和解决方案

### 1. 参数错误

#### 错误信息
```
TypeError: UniversalScRNAReader.read_10x_mtx() got an unexpected keyword argument 'delimiter'
```

#### 原因
不同格式的读取方法支持不同的参数。例如：
- `delimiter` 只用于文本格式（CSV/TSV）
- `compressed` 只用于10X MTX格式
- `genome` 只用于10X H5格式

#### 解决方案

**已修复：** v1.0.0+ 版本会自动过滤不适用的参数。

**手动解决：** 确保只传递对应格式支持的参数

```bash
# ✓ 正确：10X MTX格式使用compressed参数
python src/universal_reader.py \
  filtered_feature_bc_matrix/ \
  -o output.h5ad \
  --no-compressed

# ✓ 正确：文本格式使用delimiter参数  
python src/universal_reader.py \
  data.tsv \
  -o output.h5ad \
  --delimiter '\t'

# ✗ 错误：10X格式不支持delimiter
python src/universal_reader.py \
  filtered_feature_bc_matrix/ \
  -o output.h5ad \
  --delimiter '\t'  # 会被自动忽略（v1.0.0+）
```

#### 各格式支持的参数

| 格式 | 支持的参数 |
|------|-----------|
| **10X MTX** | `var_names`, `make_unique`, `cache`, `gex_only`, `prefix`, `compressed` |
| **10X H5** | `genome`, `gex_only` |
| **H5AD** | `backed` |
| **CSV** | `first_column_names`, `transpose` |
| **TSV/TXT** | `delimiter`, `first_column_names`, `transpose` |
| **Excel** | `sheet`, `transpose` |
| **HDF5** | `key` |
| **其他** | 无额外参数 |

---

### 2. 格式检测错误

#### 错误信息
```
ValueError: 未知或不支持的格式: unknown
```

#### 原因
自动格式检测失败或格式不支持。

#### 解决方案

```bash
# 方式1：检查文件/目录是否存在
ls -la your_data/

# 方式2：手动指定格式（未来版本支持）
# python src/universal_reader.py data/ -o output.h5ad --format 10x_mtx

# 方式3：检查是否是支持的格式
python -c "
from src.universal_reader import UniversalScRNAReader
reader = UniversalScRNAReader(False)
print(reader.supported_formats)
"
```

---

### 3. 文件权限错误

#### 错误信息
```
PermissionError: [Errno 13] Permission denied
```

#### 解决方案

```bash
# 检查文件权限
ls -l your_data/

# 修改权限
chmod -R 755 your_data/

# 或以当前用户运行
chown -R $USER:$USER your_data/
```

---

### 4. 内存不足

#### 错误信息
```
MemoryError: Unable to allocate array
```

#### 解决方案

**方式1：使用backed模式**
```python
from src.universal_reader import UniversalScRNAReader

reader = UniversalScRNAReader()
adata = reader.read_h5ad('large_data.h5ad', backed='r')

# 只加载部分到内存
subset = adata[:10000, :].to_memory()
```

**方式2：增加系统内存或使用swap**
```bash
# Linux创建swap
sudo fallocate -l 32G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

**方式3：分批处理**
```python
# 读取大文件时分批处理
for i in range(0, adata.n_obs, 10000):
    batch = adata[i:i+10000, :].to_memory()
    # 处理batch...
```

---

### 5. 依赖包缺失

#### 错误信息
```
ModuleNotFoundError: No module named 'scanpy'
```

#### 解决方案

```bash
# 确保环境已激活
conda activate scrna-datahub

# 如果还是报错，重新安装依赖
pip install -r requirements.txt

# 或安装特定包
pip install scanpy anndata pandas
```

---

### 6. 环境问题

#### 错误信息
```
command not found: conda
```

#### 解决方案

**安装Conda：**
```bash
# 下载Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# 安装
bash Miniconda3-latest-Linux-x86_64.sh

# 初始化
conda init bash
source ~/.bashrc
```

**或使用venv替代：**
```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

---

### 7. 数据读取失败

#### 错误信息
```
ValueError: Unable to read file
OSError: Unable to open file
```

#### 解决方案

**检查清单：**

1. **文件是否存在**
   ```bash
   ls -la your_data/
   ```

2. **文件是否完整**
   ```bash
   # 检查10X格式必需文件
   ls your_data/matrix.mtx*
   ls your_data/barcodes.tsv*
   ls your_data/features.tsv*
   ```

3. **文件是否损坏**
   ```bash
   # 尝试解压缩测试
   gzip -t your_data/*.gz
   ```

4. **路径是否正确**
   ```bash
   # 使用绝对路径
   python src/universal_reader.py \
     /absolute/path/to/data/ \
     -o output.h5ad
   ```

---

### 8. STARsolo输出读取失败

#### 症状
```
Error reading 10X data
```

#### 原因
STARsolo输出是**未压缩**的，但默认期望压缩文件。

#### 解决方案

```bash
# 添加 --no-compressed 参数
python src/universal_reader.py \
  Solo.out/Gene/filtered/ \
  -o output.h5ad \
  --no-compressed
```

---

### 9. 基因名重复警告

#### 警告信息
```
UserWarning: Variable names are not unique. To make them unique, call `.var_names_make_unique`.
```

#### 解决方案

**自动处理：** 工具会自动调用 `var_names_make_unique()`

**手动处理：**
```python
import scanpy as sc

adata = sc.read_h5ad('data.h5ad')
adata.var_names_make_unique()
```

---

### 10. CSV/TSV维度错误

#### 症状
维度颠倒（应该是1000 cells × 500 genes，结果是500 cells × 1000 genes）

#### 原因
矩阵需要转置

#### 解决方案

```bash
# 添加 --transpose 参数
python src/universal_reader.py \
  matrix.csv \
  -o output.h5ad \
  --transpose
```

---

## 调试技巧

### 1. 启用详细输出

```bash
# 命令行：默认已启用verbose

# Python API：
reader = UniversalScRNAReader(verbose=True)
```

### 2. 检查格式检测

```python
from src.universal_reader import UniversalScRNAReader

reader = UniversalScRNAReader(verbose=False)
format_type = reader.detect_format('your_data/')
print(f"检测到的格式: {format_type}")
```

### 3. 验证数据

```python
reader = UniversalScRNAReader()
adata = reader.read_auto('your_data/')

# 验证数据
stats = reader.validate_adata(adata)
print(stats)
```

### 4. 测试小数据集

```bash
# 先用小数据测试
python src/universal_reader.py \
  /Users/warm/华大智造/TCGA/gdc/scanpy/src/scanpy/datasets/10x_pbmc68k_reduced.h5ad \
  -o test.h5ad
```

---

## 性能问题

### 1. 读取速度慢

#### 解决方案

**使用缓存：**
```bash
python src/universal_reader.py \
  filtered_feature_bc_matrix/ \
  -o output.h5ad
  # 第一次会创建缓存，第二次读取快10-100倍
```

**使用SSD：**
- 将数据存储在SSD上而不是HDD

### 2. 内存占用高

#### 解决方案

**确保使用稀疏矩阵：**
```python
reader = UniversalScRNAReader()
adata = reader.read_auto('data/')
adata = reader.standardize_adata(adata, ensure_sparse=True)
```

**使用backed模式：**
```python
adata = reader.read_h5ad('large_data.h5ad', backed='r')
```

---

## 环境相关问题

### 1. Conda环境激活失败

```bash
# 初始化conda
conda init bash  # 或 zsh

# 重新加载
source ~/.bashrc

# 再次尝试
conda activate scrna-datahub
```

### 2. 包版本冲突

```bash
# 删除环境重新创建
conda deactivate
conda env remove -n scrna-datahub
conda env create -f environment.yml
```

### 3. pip安装失败

```bash
# 升级pip
pip install --upgrade pip setuptools wheel

# 清除缓存
pip cache purge

# 重新安装
pip install -r requirements.txt
```

---

## 获取帮助

### 自助排查

1. 查看 [常见问题FAQ](faq.md)
2. 查看 [安装指南](installation.md)
3. 查看 [环境管理指南](environment_management.md)

### 社区支持

1. **GitHub Issues**: https://github.com/yourusername/scRNA-DataHub/issues
2. **GitHub Discussions**: https://github.com/yourusername/scRNA-DataHub/discussions

### 提交Bug报告

请包含以下信息：

1. 错误信息（完整的traceback）
2. 使用的命令或代码
3. 数据格式和来源
4. 环境信息：
   ```bash
   python --version
   pip list
   # 或
   conda list
   ```
5. 操作系统信息

---

## 预防措施

### 1. 使用环境隔离

✅ **推荐：**
```bash
conda env create -f environment.yml
conda activate scrna-datahub
```

❌ **不推荐：**
```bash
pip install -r requirements.txt  # 全局安装
```

### 2. 定期更新

```bash
# 更新环境
conda activate scrna-datahub
conda update --all

# 或
pip install --upgrade -r requirements.txt
```

### 3. 测试数据验证

```bash
# 先用测试数据验证工具可用
bash scripts/quick_test.sh

# 再处理实际数据
python src/universal_reader.py your_data/ -o output.h5ad
```

---

## 总结

### 快速检查清单

- [ ] 环境已激活（`conda activate scrna-datahub`）
- [ ] 依赖已安装（`pip list | grep scanpy`）
- [ ] 文件路径正确（`ls your_data/`）
- [ ] 参数使用正确（查看上面的参数表）
- [ ] 快速测试通过（`bash scripts/quick_test.sh`）

### 最常见的3个问题

1. **环境未激活** → `conda activate scrna-datahub`
2. **参数不适用** → 已在v1.0.0+自动过滤
3. **STARsolo数据** → 添加 `--no-compressed`

**遇到问题？** 查看本文档或提交Issue！

