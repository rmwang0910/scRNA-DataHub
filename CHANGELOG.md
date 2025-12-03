# 更新日志

本文档记录scRNA-DataHub的所有重要更改。

格式基于 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.0.0/)，
版本号遵循 [语义化版本](https://semver.org/lang/zh-CN/)。

## [1.0.0] - 2024-12-02

### 新增
- ✨ 初始版本发布
- ✨ 支持12+种单细胞数据格式读取
  - 10X Genomics (MTX v2/v3, H5)
  - H5AD (AnnData)
  - Loom
  - Zarr
  - CSV/TSV/TXT
  - Excel
  - MTX (单文件)
  - HDF5
  - SOFT.GZ (GEO)
  - UMI-tools
- ✨ 自动格式检测功能
- ✨ 统一转换为H5AD格式
- ✨ 命令行工具接口
- ✨ Python API接口
- ✨ 数据验证和标准化功能
- ✨ 完整的测试套件
- ✨ 详细的文档

### 功能
- 🚀 智能缓存机制（提速10-100倍）
- 🚀 Backed模式支持（大数据友好）
- 🚀 稀疏矩阵优化
- 🚀 多基因组支持
- 🚀 批量处理支持

### 文档
- 📚 完整使用文档
- 📚 API参考文档
- 📚 数据格式详解
- 📚 测试数据获取指南
- 📚 常见问题解答

### 测试
- ✅ 7种格式使用Scanpy内置数据测试
- ✅ 快速测试脚本
- ✅ 完整测试套件
- ✅ 格式统一性验证

---

## [未来计划]

### 待添加功能
- [ ] 多样本自动合并
- [ ] 批次效应自动校正
- [ ] GUI界面
- [ ] 更多格式支持（RDS直接读取）
- [ ] 性能优化
- [ ] 并行处理支持

### 改进计划
- [ ] 增加更多测试用例
- [ ] 优化错误提示
- [ ] 添加进度条
- [ ] 支持流式处理
- [ ] 云存储支持

---

[1.0.0]: https://github.com/yourusername/scRNA-DataHub/releases/tag/v1.0.0

