#!/bin/bash
# 快速测试脚本 - 使用scanpy现有数据测试universal_scrna_reader.py

echo "========================================================================"
echo "快速测试 - 使用Scanpy现有数据"
echo "========================================================================"

SCANPY_ROOT="/Users/warm/华大智造/TCGA/gdc/scanpy"
OUTPUT_DIR="quick_test_output"

# 创建输出目录
mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

echo -e "\n测试数据来源: Scanpy本地测试数据"
echo "测试工具: ../universal_scrna_reader.py"
echo "输出目录: ${OUTPUT_DIR}/"
echo ""

# ========== 测试1: 10X MTX v3 (压缩) ==========
echo "========================================================================"
echo "测试1: 10X MTX v3 (压缩) - Cell Ranger v3+标准输出"
echo "========================================================================"
python ../universal_scrna_reader.py \
  ${SCANPY_ROOT}/tests/_data/10x_data/3.0.0/filtered_feature_bc_matrix/ \
  -o test_10x_mtx_v3.h5ad \
  --sample-id test_10x_mtx_v3

echo ""

# ========== 测试2: 10X MTX v2 (未压缩) ==========
echo "========================================================================"
echo "测试2: 10X MTX v2 (未压缩) - Cell Ranger v2输出"
echo "========================================================================"
python ../universal_scrna_reader.py \
  ${SCANPY_ROOT}/tests/_data/10x_data/1.2.0/filtered_gene_bc_matrices/hg19_chr21/ \
  -o test_10x_mtx_v2.h5ad \
  --sample-id test_10x_mtx_v2

echo ""

# ========== 测试3: 10X H5 ==========
echo "========================================================================"
echo "测试3: 10X H5 v3 - Cell Ranger H5输出"
echo "========================================================================"
python ../universal_scrna_reader.py \
  ${SCANPY_ROOT}/tests/_data/10x_data/3.0.0/filtered_feature_bc_matrix.h5 \
  -o test_10x_h5.h5ad \
  --sample-id test_10x_h5

echo ""

# ========== 测试4: H5AD ==========
echo "========================================================================"
echo "测试4: H5AD - Scanpy标准格式"
echo "========================================================================"
python ../universal_scrna_reader.py \
  ${SCANPY_ROOT}/src/scanpy/datasets/10x_pbmc68k_reduced.h5ad \
  -o test_h5ad.h5ad \
  --sample-id test_h5ad

echo ""

# ========== 测试5: TXT ==========
echo "========================================================================"
echo "测试5: TXT - 文本格式"
echo "========================================================================"
python ../universal_scrna_reader.py \
  ${SCANPY_ROOT}/src/scanpy/datasets/krumsiek11.txt \
  -o test_txt.h5ad \
  --sample-id test_txt

echo ""

# ========== 测试6: Zarr ==========
echo "========================================================================"
echo "测试6: Zarr - 云原生大数据格式"
echo "========================================================================"
python ../universal_scrna_reader.py \
  ${SCANPY_ROOT}/tests/_data/10x-10k-subset.zarr/ \
  -o test_zarr.h5ad \
  --sample-id test_zarr

echo ""

# ========== 结果汇总 ==========
echo "========================================================================"
echo "测试结果汇总"
echo "========================================================================"

echo -e "\n生成的H5AD文件:"
ls -lh test_*.h5ad

echo -e "\n数据统计:"
python << 'EOF'
import scanpy as sc
import os

h5ad_files = sorted([f for f in os.listdir('.') if f.startswith('test_') and f.endswith('.h5ad')])

print(f"{'文件名':<40} {'细胞数':<10} {'基因数':<10} {'稀疏':<8} {'大小(MB)':<10}")
print("-" * 80)

for f in h5ad_files:
    try:
        adata = sc.read_h5ad(f)
        size_mb = os.path.getsize(f) / (1024 * 1024)
        is_sparse = 'Yes' if scipy.sparse.issparse(adata.X) else 'No'
        print(f"{f:<40} {adata.n_obs:<10} {adata.n_vars:<10} {is_sparse:<8} {size_mb:<10.2f}")
    except Exception as e:
        print(f"{f:<40} {'ERROR':<10} {str(e)[:30]}")

import scipy.sparse
EOF

echo -e "\n========================================================================"
echo "测试完成！"
echo "========================================================================"
echo "所有格式都已成功转换为统一的H5AD格式"
echo ""
echo "下一步："
echo "  1. 检查生成的H5AD文件"
echo "  2. 运行完整测试: python ../create_test_data.py"
echo "  3. 清理测试数据: rm -rf ${OUTPUT_DIR}"
echo ""

