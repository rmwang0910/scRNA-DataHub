# 环境管理指南

## 概述

本文档详细说明如何创建和管理scRNA-DataHub的隔离环境，避免与其他Python项目的依赖冲突。

---

## 为什么需要环境隔离？

### 常见问题

1. **依赖冲突**：不同项目需要不同版本的包
2. **系统污染**：全局安装包影响系统Python
3. **难以复现**：无法重现相同的环境
4. **权限问题**：系统Python需要sudo权限

### 解决方案

使用虚拟环境隔离每个项目的依赖。

---

## 方式1：Conda环境（强烈推荐）⭐⭐⭐⭐⭐

### 优点

- ✅ 完全隔离
- ✅ 跨平台一致
- ✅ 包管理强大
- ✅ 支持非Python依赖
- ✅ 环境可复现

### 创建环境

```bash
cd /Users/warm/华大智造/TCGA/gdc/scRNA-DataHub

# 创建环境（根据environment.yml）
conda env create -f environment.yml

# 环境名称：scrna-datahub
# Python版本：3.8+
# 包含所有依赖
```

### 激活环境

```bash
# 激活环境
conda activate scrna-datahub

# 验证环境
which python
# 应显示：/path/to/miniconda3/envs/scrna-datahub/bin/python

# 查看已安装的包
conda list
```

### 使用环境

```bash
# 确保环境已激活
conda activate scrna-datahub

# 运行工具
python src/universal_reader.py your_data/ -o output.h5ad

# 运行测试
bash scripts/quick_test.sh

# 使用Python API
python
>>> from src.universal_reader import UniversalScRNAReader
>>> reader = UniversalScRNAReader()
```

### 退出环境

```bash
# 退出当前环境
conda deactivate
```

### 环境管理

```bash
# 查看所有环境
conda env list

# 更新环境
conda env update -f environment.yml --prune

# 导出环境
conda env export > environment_backup.yml

# 删除环境
conda deactivate
conda env remove -n scrna-datahub
```

---

## 方式2：Python venv（推荐）⭐⭐⭐⭐

### 优点

- ✅ Python标准库，无需额外安装
- ✅ 轻量级
- ✅ 快速创建

### 创建环境

```bash
cd /Users/warm/华大智造/TCGA/gdc/scRNA-DataHub

# 创建虚拟环境
python -m venv venv

# 或指定Python版本
python3.10 -m venv venv
```

### 激活环境

```bash
# Linux/macOS
source venv/bin/activate

# Windows
venv\Scripts\activate

# 验证
which python
# 应显示：/path/to/scRNA-DataHub/venv/bin/python
```

### 安装依赖

```bash
# 确保环境已激活
pip install --upgrade pip

# 安装依赖
pip install -r requirements.txt

# 验证安装
pip list
```

### 使用环境

```bash
# 确保环境已激活（提示符前会有(venv)）
(venv) $ python src/universal_reader.py --help

# 如果没有激活，先激活
source venv/bin/activate
```

### 退出环境

```bash
# 退出虚拟环境
deactivate
```

### 环境管理

```bash
# 导出依赖
pip freeze > requirements_freeze.txt

# 删除环境
deactivate
rm -rf venv/

# 重新创建
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

---

## 方式3：Pipenv（可选）⭐⭐⭐

### 优点

- ✅ 自动管理虚拟环境
- ✅ 锁定依赖版本（Pipfile.lock）
- ✅ 安全性检查

### 安装Pipenv

```bash
pip install pipenv
```

### 创建环境

```bash
cd /Users/warm/华大智造/TCGA/gdc/scRNA-DataHub

# 创建环境并安装依赖
pipenv install -r requirements.txt

# 激活环境
pipenv shell
```

### 使用

```bash
# 运行命令（自动使用环境）
pipenv run python src/universal_reader.py --help

# 或进入shell
pipenv shell
python src/universal_reader.py --help
```

---

## 环境对比

| 特性 | Conda | venv | Pipenv |
|------|-------|------|--------|
| **安装难度** | 需要安装conda | Python内置 | 需要安装pipenv |
| **隔离程度** | 完全隔离 | Python包隔离 | Python包隔离 |
| **跨平台** | ✅ 优秀 | ✅ 良好 | ✅ 良好 |
| **包管理** | conda + pip | pip | pip |
| **依赖锁定** | ✅ | ❌ | ✅ |
| **速度** | 较慢 | 快 | 中等 |
| **推荐度** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ |

---

## 常见场景

### 场景1：日常使用

```bash
# 每次使用前
cd /Users/warm/华大智造/TCGA/gdc/scRNA-DataHub
conda activate scrna-datahub

# 使用工具
python src/universal_reader.py ...

# 使用完毕
conda deactivate
```

### 场景2：多项目切换

```bash
# 项目A
conda activate project-a
# 工作...
conda deactivate

# 切换到scRNA-DataHub
conda activate scrna-datahub
# 工作...
conda deactivate

# 项目B
conda activate project-b
# 工作...
```

### 场景3：环境共享

```bash
# 导出环境配置
conda env export > environment.yml

# 在其他机器上重现
conda env create -f environment.yml
```

---

## 故障排除

### 问题1：conda命令找不到

**解决：** 安装Miniconda或Anaconda

```bash
# 下载Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# 安装
bash Miniconda3-latest-Linux-x86_64.sh

# 重新加载shell
source ~/.bashrc
```

### 问题2：环境创建失败

**解决：** 检查environment.yml文件

```bash
# 查看文件内容
cat environment.yml

# 手动创建环境
conda create -n scrna-datahub python=3.10
conda activate scrna-datahub
pip install -r requirements.txt
```

### 问题3：包冲突

**解决：** 删除环境重新创建

```bash
conda deactivate
conda env remove -n scrna-datahub
conda env create -f environment.yml
```

### 问题4：环境激活失败

**解决：** 初始化conda

```bash
# 初始化conda（只需一次）
conda init bash  # 或 zsh, fish等

# 重新加载shell
source ~/.bashrc  # 或 ~/.zshrc
```

---

## 最佳实践

### 1. 始终使用环境

❌ **不推荐：**
```bash
pip install scanpy  # 安装到系统Python
```

✅ **推荐：**
```bash
conda activate scrna-datahub
pip install scanpy  # 安装到隔离环境
```

### 2. 定期更新环境

```bash
# 更新环境
conda activate scrna-datahub
conda update --all

# 或更新特定包
pip install --upgrade scanpy
```

### 3. 保存环境快照

```bash
# 定期导出环境
conda env export > environment_$(date +%Y%m%d).yml
```

### 4. 项目切换

```bash
# 使用前激活
conda activate scrna-datahub

# 使用后退出
conda deactivate
```

---

## 快速参考

### Conda常用命令

```bash
# 创建环境
conda env create -f environment.yml

# 激活环境
conda activate scrna-datahub

# 退出环境
conda deactivate

# 查看环境
conda env list

# 删除环境
conda env remove -n scrna-datahub

# 更新环境
conda env update -f environment.yml

# 导出环境
conda env export > environment.yml
```

### venv常用命令

```bash
# 创建环境
python -m venv venv

# 激活环境（Linux/macOS）
source venv/bin/activate

# 激活环境（Windows）
venv\Scripts\activate

# 退出环境
deactivate

# 删除环境
rm -rf venv/
```

---

## 推荐配置

### 推荐方案：Conda

```bash
# 一次性设置
cd /Users/warm/华大智造/TCGA/gdc/scRNA-DataHub
conda env create -f environment.yml

# 每次使用
conda activate scrna-datahub
python src/universal_reader.py ...
conda deactivate
```

### 添加到shell配置（可选）

```bash
# 添加到 ~/.bashrc 或 ~/.zshrc
alias scrna-datahub='cd /Users/warm/华大智造/TCGA/gdc/scRNA-DataHub && conda activate scrna-datahub'

# 使用
scrna-datahub  # 自动进入项目目录并激活环境
```

---

## 总结

### 推荐使用Conda环境的理由

1. ✅ **完全隔离** - 避免包冲突
2. ✅ **可复现** - environment.yml保证一致性
3. ✅ **易管理** - conda命令简单直观
4. ✅ **跨平台** - Linux/macOS/Windows一致
5. ✅ **专业** - 科学计算标准

### 基本工作流程

```bash
# 1. 一次性创建环境
conda env create -f environment.yml

# 2. 每次使用前激活
conda activate scrna-datahub

# 3. 使用工具
python src/universal_reader.py ...

# 4. 使用完毕退出
conda deactivate
```

**开始使用：**
```bash
conda env create -f environment.yml
conda activate scrna-datahub
bash scripts/quick_test.sh
```

