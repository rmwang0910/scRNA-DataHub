#!/bin/bash
# 测试DNB C4格式读取

echo "========================================================================"
echo "测试DNB C4格式读取"
echo "========================================================================"

DNB_DATA="/storeData/ztron/wangrm/zzTCGA/pipeline/scRNA/results/CNS1063416_brain/02.count/filter_matrix"
OUTPUT="CNS1063416_brain.h5ad"

cd /storeData/ztron/wangrm/tools/scRNA-DataHub

echo -e "\n1. 检查DNB C4数据文件..."
echo "数据目录: $DNB_DATA"

if [ ! -d "$DNB_DATA" ]; then
    echo "❌ 数据目录不存在"
    exit 1
fi

echo "文件列表:"
ls -lh "$DNB_DATA"

echo -e "\n2. 检查features.tsv格式..."
echo "前5行features.tsv:"
zcat "$DNB_DATA/features.tsv.gz" | head -5

echo -e "\n列数检查:"
n_cols=$(zcat "$DNB_DATA/features.tsv.gz" | head -1 | awk -F'\t' '{print NF}')
echo "features.tsv列数: $n_cols"

if [ "$n_cols" -eq 1 ]; then
    echo "✓ 确认为DNB C4格式（只有1列基因名）"
else
    echo "⚠️  不是标准DNB C4格式（有${n_cols}列）"
fi

echo -e "\n3. 使用scRNA-DataHub读取..."
python src/universal_reader.py \
  "$DNB_DATA" \
  -o "$OUTPUT" \
  --sample-id CNS1063416_brain

if [ $? -eq 0 ]; then
    echo -e "\n✓ 读取成功！"
    
    echo -e "\n4. 验证输出..."
    python << 'EOF'
import scanpy as sc

adata = sc.read_h5ad('CNS1063416_brain.h5ad')

print(f"数据维度: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"\nvar列: {list(adata.var.columns)}")
print(f"\n前5个基因:")
print(adata.var.head())

# 检查必需列
assert 'gene_ids' in adata.var.columns, "缺少gene_ids列"
assert 'feature_types' in adata.var.columns, "缺少feature_types列"

print("\n✓ 所有必需列都存在")
print("✓ DNB C4格式读取成功！")
EOF

else
    echo -e "\n❌ 读取失败"
    exit 1
fi

echo -e "\n========================================================================"
echo "测试完成"
echo "========================================================================"

