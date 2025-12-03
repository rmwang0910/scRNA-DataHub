# 贡献指南

感谢您对 scRNA-DataHub 的关注！我们欢迎各种形式的贡献。

## 如何贡献

### 报告Bug

如果您发现了bug，请在 [GitHub Issues](https://github.com/yourusername/scRNA-DataHub/issues) 创建issue，并包含：

1. Bug描述
2. 复现步骤
3. 预期行为
4. 实际行为
5. 环境信息（Python版本、操作系统等）
6. 错误日志（如果有）

### 提出新功能

如果您有新功能建议，请创建issue并描述：

1. 功能描述
2. 使用场景
3. 预期实现方式
4. 示例代码（如果有）

### 提交代码

1. **Fork 仓库**
   ```bash
   # 在GitHub上fork本仓库
   ```

2. **克隆您的fork**
   ```bash
   git clone https://github.com/yourusername/scRNA-DataHub.git
   cd scRNA-DataHub
   ```

3. **创建分支**
   ```bash
   git checkout -b feature/your-feature-name
   ```

4. **进行更改**
   - 遵循代码风格
   - 添加测试
   - 更新文档

5. **运行测试**
   ```bash
   pytest tests/
   ```

6. **提交更改**
   ```bash
   git add .
   git commit -m "Add feature: your feature description"
   ```

7. **推送到GitHub**
   ```bash
   git push origin feature/your-feature-name
   ```

8. **创建Pull Request**
   - 在GitHub上创建PR
   - 描述您的更改
   - 等待review

## 代码规范

### Python代码风格

- 遵循 PEP 8
- 使用 Black 格式化代码
- 使用类型提示（Type Hints）
- 添加docstrings

```python
def example_function(param1: str, param2: int) -> bool:
    """
    函数简短描述
    
    参数:
        param1: 参数1描述
        param2: 参数2描述
        
    返回:
        返回值描述
    """
    pass
```

### 测试

- 为新功能添加测试
- 确保所有测试通过
- 测试覆盖率 > 80%

```python
def test_new_feature():
    """测试新功能"""
    # 测试代码
    assert expected == actual
```

### 文档

- 更新相关文档
- 添加使用示例
- 更新CHANGELOG.md

## 开发环境设置

```bash
# 克隆仓库
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub

# 创建虚拟环境
conda env create -f environment.yml
conda activate scrna-datahub

# 安装开发依赖
pip install -e ".[dev]"

# 运行测试
pytest tests/
```

## 行为准则

- 尊重他人
- 建设性反馈
- 专注于问题本身
- 保持友好和专业

## 问题和支持

- 技术问题：[GitHub Issues](https://github.com/yourusername/scRNA-DataHub/issues)
- 功能建议：[GitHub Discussions](https://github.com/yourusername/scRNA-DataHub/discussions)

感谢您的贡献！

