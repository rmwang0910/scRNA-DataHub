#!/usr/bin/env python3
"""
多样本合并示例
演示如何读取多个样本并合并为一个AnnData对象
"""

import sys
sys.path.insert(0, '../src')

from universal_reader import UniversalScRNAReader
import scanpy as sc
import pandas as pd


def merge_multiple_samples():
    """读取并合并多个样本"""
    
    print("=" * 70)
    print("多样本合并示例")
    print("=" * 70)
    
    reader = UniversalScRNAReader(verbose=True)
    
    # 定义样本信息
    samples = {
        'sample1': {
            'path': 'data/sample1/filtered_feature_bc_matrix/',
            'condition': 'control',
            'batch': 'batch1'
        },
        'sample2': {
            'path': 'data/sample2/filtered_feature_bc_matrix/',
            'condition': 'treatment',
            'batch': 'batch1'
        },
        'sample3': {
            'path': 'data/sample3/filtered_feature_bc_matrix/',
            'condition': 'treatment',
            'batch': 'batch2'
        }
    }
    
    # 步骤1：读取所有样本
    print("\n步骤1：读取样本...")
    adatas = []
    
    for sample_id, info in samples.items():
        # 读取数据
        adata = reader.read_auto(info['path'])
        
        # 添加元数据
        adata.obs['sample_id'] = sample_id
        adata.obs['condition'] = info['condition']
        adata.obs['batch'] = info['batch']
        
        # 确保条形码唯一
        adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]
        
        # 标准化
        adata = reader.standardize_adata(adata)
        
        adatas.append(adata)
        print(f"  ✓ {sample_id}: {adata.n_obs} cells")
    
    # 步骤2：合并样本
    print("\n步骤2：合并样本...")
    adata_combined = sc.concat(
        adatas,
        axis=0,           # 按细胞合并
        join='outer',     # 保留所有基因
        merge='unique'    # 处理重复列
    )
    
    print(f"  合并后: {adata_combined.n_obs} cells × {adata_combined.n_vars} genes")
    
    # 步骤3：质量控制
    print("\n步骤3：质量控制...")
    sc.pp.filter_cells(adata_combined, min_genes=200)
    sc.pp.filter_genes(adata_combined, min_cells=3)
    
    # 步骤4：保存
    print("\n步骤4：保存合并数据...")
    adata_combined.write('merged_samples.h5ad')
    
    print("\n✓ 多样本合并完成")
    
    return adata_combined


def merge_with_batch_correction():
    """多样本合并 + 批次校正"""
    
    print("\n" + "=" * 70)
    print("多样本合并 + 批次校正示例")
    print("=" * 70)
    
    # 先合并
    adata_combined = merge_multiple_samples()
    
    # 预处理
    print("\n预处理...")
    sc.pp.normalize_total(adata_combined, target_sum=1e4)
    sc.pp.log1p(adata_combined)
    sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000, batch_key='batch')
    sc.tl.pca(adata_combined, n_comps=50)
    
    # 批次校正（使用Harmony）
    print("\n批次校正（Harmony）...")
    try:
        import harmonypy
        harmony_out = harmonypy.run_harmony(
            adata_combined.obsm['X_pca'],
            adata_combined.obs,
            'batch',
            max_iter_harmony=20
        )
        adata_combined.obsm['X_pca_harmony'] = harmony_out.Z_corr.T
        use_rep = 'X_pca_harmony'
        print("✓ Harmony校正完成")
    except ImportError:
        print("⚠️  Harmony未安装，跳过批次校正")
        use_rep = 'X_pca'
    
    # 聚类和可视化
    print("\n聚类和可视化...")
    sc.pp.neighbors(adata_combined, use_rep=use_rep)
    sc.tl.umap(adata_combined)
    sc.tl.leiden(adata_combined)
    
    # 保存
    adata_combined.write('merged_batch_corrected.h5ad')
    
    # 可视化
    sc.pl.umap(adata_combined, color=['batch', 'condition', 'leiden'], save='_merged.png')
    
    print("\n✓ 批次校正完成")


if __name__ == '__main__':
    print("多样本合并示例\n")
    
    # 运行示例（需要有相应的数据文件）
    # merge_multiple_samples()
    # merge_with_batch_correction()
    
    print("查看代码了解详细用法")

