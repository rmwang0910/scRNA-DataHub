#!/usr/bin/env python3
"""
批量处理示例
演示如何批量处理多个单细胞样本
"""

import sys
sys.path.insert(0, '../src')

from universal_reader import UniversalScRNAReader
from pathlib import Path


def example_batch_processing():
    """批量处理多个样本"""
    
    print("=" * 70)
    print("批量处理示例")
    print("=" * 70)
    
    reader = UniversalScRNAReader(verbose=True)
    
    # 定义样本列表
    samples = [
        {
            'input': 'data/sample1/filtered_feature_bc_matrix/',
            'output': 'output/sample1.h5ad',
            'sample_id': 'sample1'
        },
        {
            'input': 'data/sample2/filtered_feature_bc_matrix/',
            'output': 'output/sample2.h5ad',
            'sample_id': 'sample2'
        },
        {
            'input': 'data/sample3/filtered_feature_bc_matrix/',
            'output': 'output/sample3.h5ad',
            'sample_id': 'sample3'
        }
    ]
    
    # 批量处理
    for sample in samples:
        print(f"\n处理 {sample['sample_id']}...")
        
        adata = reader.process_single_sample(
            input_path=sample['input'],
            output_path=sample['output'],
            sample_id=sample['sample_id']
        )
        
        print(f"✓ {sample['sample_id']} 完成")
    
    print("\n✓ 所有样本处理完成")


def example_auto_discovery():
    """自动发现和批量处理"""
    
    print("\n" + "=" * 70)
    print("自动发现样本并批量处理")
    print("=" * 70)
    
    reader = UniversalScRNAReader()
    
    # 自动搜索所有10X样本目录
    data_root = Path('data/')
    sample_dirs = data_root.glob('*/outs/filtered_feature_bc_matrix/')
    
    for sample_dir in sample_dirs:
        # 从路径提取样本名
        sample_id = sample_dir.parent.parent.name
        output_path = f'output/{sample_id}.h5ad'
        
        print(f"\n处理 {sample_id}...")
        
        adata = reader.process_single_sample(
            input_path=str(sample_dir),
            output_path=output_path,
            sample_id=sample_id
        )
        
        print(f"✓ {sample_id}: {adata.n_obs} cells × {adata.n_vars} genes")
    
    print("\n✓ 自动发现并处理完成")


def example_mixed_formats():
    """处理混合格式的样本"""
    
    print("\n" + "=" * 70)
    print("处理混合格式样本")
    print("=" * 70)
    
    reader = UniversalScRNAReader()
    
    # 不同格式的样本
    samples = [
        {
            'id': 'sample1',
            'path': 'data/sample1/filtered_feature_bc_matrix/',  # 10X MTX
            'kwargs': {}
        },
        {
            'id': 'sample2',
            'path': 'data/sample2/filtered_feature_bc_matrix.h5',  # 10X H5
            'kwargs': {}
        },
        {
            'id': 'sample3',
            'path': 'data/sample3/expression.csv',  # CSV
            'kwargs': {'transpose': True}
        }
    ]
    
    for sample in samples:
        print(f"\n处理 {sample['id']} ({reader.detect_format(sample['path'])})...")
        
        adata = reader.read_auto(sample['path'], **sample['kwargs'])
        adata.obs['sample_id'] = sample['id']
        adata = reader.standardize_adata(adata)
        
        reader.save_h5ad(adata, f"output/{sample['id']}.h5ad")
        
        print(f"✓ {sample['id']} 完成")
    
    print("\n✓ 混合格式处理完成")


if __name__ == '__main__':
    print("scRNA-DataHub 批量处理示例\n")
    
    # 运行示例
    # example_batch_processing()
    # example_auto_discovery()
    # example_mixed_formats()
    
    print("查看代码了解详细用法")

