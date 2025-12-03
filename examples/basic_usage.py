#!/usr/bin/env python3
"""
基础使用示例
演示如何使用scRNA-DataHub读取单细胞数据
"""

import sys
sys.path.insert(0, '../src')

from universal_reader import UniversalScRNAReader


def example_basic_usage():
    """基础使用示例"""
    
    print("=" * 70)
    print("示例1：基本使用")
    print("=" * 70)
    
    # 创建读取器
    reader = UniversalScRNAReader(verbose=True)
    
    # 自动读取任意格式
    adata = reader.read_auto('your_data/')
    
    # 验证数据
    stats = reader.validate_adata(adata)
    
    # 标准化
    adata = reader.standardize_adata(adata)
    
    # 保存
    reader.save_h5ad(adata, 'output.h5ad')
    
    print("\n✓ 完成")


def example_10x_genomics():
    """读取10X Genomics数据"""
    
    print("\n" + "=" * 70)
    print("示例2：读取10X Genomics数据")
    print("=" * 70)
    
    reader = UniversalScRNAReader()
    
    # Cell Ranger v3+ 输出
    adata = reader.read_10x_mtx(
        'sample1/outs/filtered_feature_bc_matrix/',
        var_names='gene_symbols',
        cache=True
    )
    
    # 添加样本信息
    adata.obs['sample_id'] = 'sample1'
    
    # 保存
    reader.save_h5ad(adata, 'sample1.h5ad')
    
    print("\n✓ 完成")


def example_starsolo():
    """读取STARsolo输出"""
    
    print("\n" + "=" * 70)
    print("示例3：读取STARsolo输出")
    print("=" * 70)
    
    reader = UniversalScRNAReader()
    
    # STARsolo输出（未压缩）
    adata = reader.read_10x_mtx(
        'Solo.out/Gene/filtered/',
        var_names='gene_symbols',
        compressed=False  # 关键：STARsolo不压缩
    )
    
    reader.save_h5ad(adata, 'starsolo_output.h5ad')
    
    print("\n✓ 完成")


def example_csv():
    """读取CSV矩阵"""
    
    print("\n" + "=" * 70)
    print("示例4：读取CSV矩阵")
    print("=" * 70)
    
    reader = UniversalScRNAReader()
    
    # CSV文件（行是基因，列是细胞）
    adata = reader.read_csv(
        'expression_matrix.csv',
        first_column_names=True,
        transpose=True  # 转置为 cells × genes
    )
    
    reader.save_h5ad(adata, 'from_csv.h5ad')
    
    print("\n✓ 完成")


def example_one_liner():
    """一行代码处理"""
    
    print("\n" + "=" * 70)
    print("示例5：一行代码处理")
    print("=" * 70)
    
    reader = UniversalScRNAReader()
    
    # 一键处理：读取 → 验证 → 标准化 → 保存
    adata = reader.process_single_sample(
        input_path='filtered_feature_bc_matrix/',
        output_path='sample1.h5ad',
        sample_id='sample1',
        var_names='gene_symbols',
        cache=True
    )
    
    print(f"\n最终数据: {adata}")
    print("✓ 完成")


if __name__ == '__main__':
    print("scRNA-DataHub 基础使用示例")
    print("=" * 70)
    
    # 运行所有示例
    # 注意：需要有相应的数据文件
    
    # example_basic_usage()
    # example_10x_genomics()
    # example_starsolo()
    # example_csv()
    example_one_liner()

