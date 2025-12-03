# 安装指南

## 系统要求

- Python 3.8 或更高版本
- 支持的操作系统：Linux, macOS, Windows (WSL)
- 推荐内存：8GB 以上

## 安装方式

### 方式1：使用pip（推荐）

```bash
# 克隆仓库
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub

# 安装
pip install -e .

# 或安装完整版本（包含所有可选依赖）
pip install -e ".[full]"
```

### 方式2：使用conda

```bash
# 克隆仓库
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub

# 创建conda环境
conda env create -f environment.yml

# 激活环境
conda activate scrna-datahub
```

### 方式3：只安装依赖

```bash
# 克隆仓库
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub

# 安装核心依赖
pip install -r requirements.txt
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

