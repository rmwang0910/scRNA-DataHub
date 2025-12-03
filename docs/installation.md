# 安装指南

## 系统要求

- Python 3.8 或更高版本
- 支持的操作系统：Linux, macOS, Windows (WSL)
- 推荐内存：8GB 以上

## 安装方式

### 方式1：使用Conda环境（强烈推荐）⭐

**优点：** 完全隔离，避免包冲突，环境可复现

```bash
# 克隆仓库
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub

# 创建独立的conda环境
conda env create -f environment.yml

# 激活环境
conda activate scrna-datahub

# 验证安装
python src/universal_reader.py --help
python -c "from src.universal_reader import UniversalScRNAReader; print('✓ 安装成功')"
```

**卸载：**
```bash
conda deactivate
conda env remove -n scrna-datahub
```

### 方式2：使用Python venv（推荐）

**优点：** 不需要conda，使用Python标准库

```bash
# 克隆仓库
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub

# 创建Python虚拟环境
python -m venv venv

# 激活虚拟环境
source venv/bin/activate  # Linux/macOS
# 或
venv\Scripts\activate     # Windows

# 安装依赖
pip install -r requirements.txt

# 验证安装
python src/universal_reader.py --help
```

**卸载：**
```bash
deactivate
rm -rf venv/
```

### 方式3：使用pip全局安装（不推荐）

**注意：** 可能与其他Python包冲突，仅用于测试

```bash
# 克隆仓库
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub

# 安装
pip install -e .

# 或只安装依赖
pip install -r requirements.txt
```

### 方式4：开发模式（开发者）

```bash
# 克隆仓库
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub

# 创建conda环境
conda env create -f environment.yml
conda activate scrna-datahub

# 安装开发依赖
pip install -e ".[dev]"

# 安装pre-commit hooks
pip install pre-commit
pre-commit install
```

---

## 环境管理最佳实践

### 推荐工作流程

```bash
# 1. 创建环境（只需一次）
conda env create -f environment.yml

# 2. 每次使用前激活
conda activate scrna-datahub

# 3. 使用工具
python src/universal_reader.py your_data/ -o output.h5ad

# 4. 使用完毕后退出
conda deactivate
```

### 环境列表

```bash
# 查看所有conda环境
conda env list

# 查看当前环境
conda info --envs

# 更新环境
conda env update -f environment.yml --prune
```

### 环境导出

```bash
# 导出当前环境
conda env export > environment_backup.yml

# 或只导出手动安装的包
conda env export --from-history > environment_minimal.yml
```

## 验证安装

```bash
# 测试命令行工具
python src/universal_reader.py --help

# 测试Python API
python -c "from src.universal_reader import UniversalScRNAReader; print('✓ 安装成功')"

# 运行快速测试
bash scripts/quick_test.sh
```

## 可选依赖

### 格式支持

```bash
# Loom格式
pip install loompy

# Zarr格式
pip install zarr

# Excel格式
pip install openpyxl
```

### 批次校正

```bash
# Harmony
pip install harmonypy

# BBKNN
pip install bbknn

# Scanorama
pip install scanorama
```

## 开发环境设置

```bash
# 克隆仓库
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub

# 安装开发依赖
pip install -e ".[dev]"

# 安装pre-commit hooks
pre-commit install
```

## 故障排除

### 问题1：安装失败

```bash
# 升级pip
pip install --upgrade pip setuptools wheel

# 重新安装
pip install -e .
```

### 问题2：缺少依赖

```bash
# 安装所有依赖
pip install -r requirements.txt
pip install loompy zarr openpyxl
```

### 问题3：权限问题

```bash
# 使用--user安装
pip install --user -e .
```

## 卸载

```bash
# pip安装的
pip uninstall scrna-datahub

# conda环境
conda deactivate
conda env remove -n scrna-datahub
```

