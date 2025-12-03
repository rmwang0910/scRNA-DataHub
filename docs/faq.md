# 常见问题

## 安装相关

### Q1: 如何安装scRNA-DataHub？

```bash
git clone https://github.com/yourusername/scRNA-DataHub.git
cd scRNA-DataHub
pip install -r requirements.txt
```

详见 [安装指南](installation.md)

### Q2: 缺少某些依赖怎么办？

```bash
# 安装完整依赖
pip install loompy zarr openpyxl harmonypy bbknn scanorama
```

---

## 使用相关

### Q3: 如何自动检测数据格式？

工具会自动检测格式，无需手动指定：

```bash
python src/universal_reader.py your_data/ -o output.h5ad
```

### Q4: STARsolo输出读取失败？

STARsolo输出是**未压缩**的，需要添加 `--no-compressed` 参数：

```bash
python src/universal_reader.py \
  Solo.out/Gene/filtered/ \
  -o output.h5ad \
  --no-compressed
```

### Q5: 如何处理人鼠混合样本？

使用 `--genome` 参数指定基因组：

```bash
# 只读人类基因组
python src/universal_reader.py \
  data.h5 \
  -o human.h5ad \
  --genome GRCh38

# 只读小鼠基因组
python src/universal_reader.py \
  data.h5 \
  -o mouse.h5ad \
  --genome mm10
```

### Q6: CSV/TSV文件维度不对？

如果行是基因、列是细胞，需要转置：

```bash
python src/universal_reader.py \
  matrix.csv \
  -o output.h5ad \
  --transpose
```

### Q7: 内存不足怎么办？

使用backed模式：

```python
from src.universal_reader import UniversalScRNAReader

reader = UniversalScRNAReader()
adata = reader.read_h5ad('large_data.h5ad', backed='r')

# 只加载部分数据到内存
subset = adata[:10000, :].to_memory()
```

### Q8: 如何批量处理多个样本？

```python
from src.universal_reader import UniversalScRNAReader

reader = UniversalScRNAReader()

samples = ['sample1', 'sample2', 'sample3']
for sample in samples:
    reader.process_single_sample(
        input_path=f'{sample}/filtered_feature_bc_matrix/',
        output_path=f'{sample}.h5ad',
        sample_id=sample
    )
```

---

## 数据格式相关

### Q9: 支持哪些数据格式？

支持12+种格式：
- 10X Genomics (MTX, H5)
- H5AD, Loom, Zarr
- CSV, TSV, Excel
- MTX, HDF5, SOFT.GZ, UMI-tools

详见 [数据格式详解](data_formats.md)

### Q10: RDS/H5Seurat格式如何处理？

需要先在R中转换：

```r
library(SeuratDisk)
Convert("data.rds", dest = "h5ad")
```

然后读取转换后的H5AD文件。

### Q11: 如何获取测试数据？

参考 [测试数据获取指南](test_data_guide.md)

---

## 性能相关

### Q12: 读取速度太慢？

使用缓存机制：

```python
# 第一次会创建缓存
adata = reader.read_10x_mtx('dir/', cache=True)

# 后续读取快10-100倍
adata = reader.read_10x_mtx('dir/', cache=True)
```

### Q13: 如何处理超大数据（>100万细胞）？

1. 使用Zarr格式
2. 使用backed模式
3. 分批处理

```python
# 使用backed模式
adata = reader.read_h5ad('huge_data.h5ad', backed='r')

# 分批处理
for i in range(0, adata.n_obs, 100000):
    batch = adata[i:i+100000, :].to_memory()
    # 处理batch...
```

---

## 错误处理

### Q14: "index 0 is out of bounds" 错误？

可能是压缩设置不对。如果是STARsolo输出，添加 `--no-compressed`。

### Q15: "Variable names are not unique" 警告？

工具会自动处理重复基因名，添加后缀-1, -2等。

### Q16: 基因名不匹配？

确保使用正确的var_names设置：

```bash
# 使用基因符号
python src/universal_reader.py dir/ -o out.h5ad --var-names gene_symbols

# 使用基因ID
python src/universal_reader.py dir/ -o out.h5ad --var-names gene_ids
```

---

## 其他问题

### Q17: 如何贡献代码？

查看 [贡献指南](../CONTRIBUTING.md)

### Q18: 如何报告bug？

在 [GitHub Issues](https://github.com/yourusername/scRNA-DataHub/issues) 创建issue

### Q19: 在哪里提问？

- GitHub Issues: Bug报告
- GitHub Discussions: 使用问题和讨论
- Email: 技术支持

### Q20: 项目更新频率？

查看 [更新日志](../CHANGELOG.md)

---

## 联系我们

- GitHub: https://github.com/yourusername/scRNA-DataHub
- Issues: https://github.com/yourusername/scRNA-DataHub/issues
- Email: your.email@example.com

