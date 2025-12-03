# 版本兼容性说明

## Scanpy版本兼容性

scRNA-DataHub 支持 Scanpy 1.9.0 及以上版本，并会自动适配不同版本的API差异。

---

## 版本特性

### Scanpy 1.9.3+

**新增功能：**
- ✅ `read_10x_mtx()` 支持 `compressed` 参数
- ✅ 可以明确指定文件是否压缩

**使用示例：**
```python
import scanpy as sc

# Cell Ranger v3+（压缩）
adata = sc.read_10x_mtx('dir/', compressed=True)

# STARsolo（未压缩）
adata = sc.read_10x_mtx('dir/', compressed=False)
```

### Scanpy 1.9.0 - 1.9.2

**限制：**
- ❌ `read_10x_mtx()` 不支持 `compressed` 参数
- ℹ️ 会自动检测文件是否压缩

**自动适配：**

scRNA-DataHub会自动检测scanpy版本并适配：

```python
# 工具会自动判断
if supports_compressed:
    # 新版本：使用compressed参数
    adata = sc.read_10x_mtx(path, compressed=True)
else:
    # 旧版本：自动检测
    adata = sc.read_10x_mtx(path)
```

---

## 版本检查

### 查看当前版本

```bash
# 命令行
python -c "import scanpy; print(scanpy.__version__)"

# Python
import scanpy as sc
print(f"Scanpy version: {sc.__version__}")
```

### 升级Scanpy

```bash
# 升级到最新版本
pip install --upgrade scanpy

# 或使用conda
conda update scanpy
```

---

## 兼容性矩阵

| Scanpy版本 | scRNA-DataHub | compressed参数 | 推荐度 |
|-----------|--------------|---------------|--------|
| **1.10.0+** | ✅ 完全兼容 | ✅ 支持 | ⭐⭐⭐⭐⭐ |
| **1.9.3+** | ✅ 完全兼容 | ✅ 支持 | ⭐⭐⭐⭐⭐ |
| **1.9.0-1.9.2** | ✅ 兼容 | ❌ 不支持（自动适配） | ⭐⭐⭐⭐ |
| **<1.9.0** | ⚠️ 未测试 | ❌ 不支持 | ⭐⭐ |

---

## 推荐配置

### 最新稳定版（推荐）

```bash
# 安装最新stable版本
pip install 'scanpy>=1.10.0'
```

### 最小要求版本

```bash
# 最低要求
pip install 'scanpy>=1.9.0'
```

---

## 特殊情况处理

### 情况1：使用旧版本Scanpy

如果必须使用旧版本scanpy：

```bash
# 不要使用 --no-compressed 参数
python src/universal_reader.py \
  your_data/ \
  -o output.h5ad

# 工具会自动检测文件格式
```

### 情况2：升级Scanpy后的变化

```python
# 升级前（scanpy 1.9.2）
# compressed参数会被忽略

# 升级后（scanpy 1.9.3+）
# compressed参数会被使用
```

**注意：** 工具会自动适配，无需修改代码。

---

## 依赖版本要求

### 核心依赖

```
scanpy>=1.9.0
anndata>=0.8.0
pandas>=1.3.0
numpy>=1.20.0
scipy>=1.7.0
h5py>=3.0.0
```

### Python版本

```
python>=3.8
```

**推荐组合：**
- Python 3.10 + Scanpy 1.10.0

---

## 自动兼容性检测

scRNA-DataHub 包含自动兼容性检测：

```python
import inspect
import scanpy as sc

# 检查是否支持compressed参数
sig = inspect.signature(sc.read_10x_mtx)
supports_compressed = 'compressed' in sig.parameters

if supports_compressed:
    print("✓ 支持compressed参数")
else:
    print("⚠️  不支持compressed参数，将使用自动检测")
```

这个检测是透明的，用户无需关心。

---

## 故障排除

### 问题：compressed参数报错

**症状：**
```
TypeError: read_10x_mtx() got an unexpected keyword argument 'compressed'
```

**原因：**
- Scanpy版本 < 1.9.3
- 或代码未更新

**解决方案：**

**方式1：升级Scanpy（推荐）**
```bash
pip install --upgrade 'scanpy>=1.9.3'
```

**方式2：更新代码**

确保使用最新版本的 `universal_reader.py`（已包含版本检测）。

**方式3：不使用compressed参数**
```bash
# 让scanpy自动检测
python src/universal_reader.py your_data/ -o output.h5ad
# 不要添加 --no-compressed
```

---

## 版本更新日志

### scRNA-DataHub v1.0.1（当前）

**改进：**
- ✅ 自动检测scanpy版本
- ✅ 兼容compressed参数的有无
- ✅ 向后兼容旧版本scanpy

### scRNA-DataHub v1.0.0

**要求：**
- Scanpy >= 1.9.3（带compressed参数）

---

## 总结

### 推荐配置

```bash
# 创建环境时指定版本
conda create -n scrna-datahub python=3.10
conda activate scrna-datahub
pip install 'scanpy>=1.10.0' anndata pandas numpy scipy h5py
```

### 兼容性保证

scRNA-DataHub v1.0.1+ 会自动适配：
- ✅ Scanpy 1.9.0+
- ✅ 自动检测API变化
- ✅ 向后兼容

**无需担心版本问题！** ✅

