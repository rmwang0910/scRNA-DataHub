#!/usr/bin/env python3
"""
基础使用示例
演示如何使用 scRNA-DataHub 读取单细胞数据
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
    
    # 自动读取任意格式（会自动检测格式）
    adata = reader.read_auto('your_data/')
    
    # 验证数据
    stats = reader.validate_adata(adata)
    print(f"\n数据统计: {stats}")
    
    # 标准化
    adata = reader.standardize_adata(adata)
    
    # 保存
    reader.save_h5ad(adata, 'output.h5ad')
    
    print("\n✓ 完成！文件已保存到 output.h5ad")


def example_10x_genomics():
    """读取 10X Genomics 数据"""
    
    print("\n" + "=" * 70)
    print("示例2：读取 10X Genomics 数据")
    print("=" * 70)
    
    reader = UniversalScRNAReader(verbose=True)
    
    # Cell Ranger v3+ 输出目录
    adata = reader.read_10x_mtx(
        'sample1/outs/filtered_feature_bc_matrix/',
        var_names='gene_symbols',
        cache=True
    )
    
    # 添加样本信息
    adata.obs['sample_id'] = 'sample1'
    
    # 保存
    reader.save_h5ad(adata, 'sample1.h5ad')
    
    print("\n✓ 完成！")


def example_dnb_c4():
    """读取 DNB C4 数据"""
    
    print("\n" + "=" * 70)
    print("示例3：读取 DNB C4 数据")
    print("=" * 70)
    
    reader = UniversalScRNAReader(verbose=True)
    
    # dnbc4tools 输出目录（与 10X 格式兼容）
    adata = reader.read_10x_mtx(
        'CNS1063416_brain/02.count/filter_matrix/',
        var_names='gene_symbols',
        cache=True
    )
    
    # 添加样本信息
    adata.obs['sample_id'] = 'CNS1063416_brain'
    
    reader.save_h5ad(adata, 'dnb_c4_output.h5ad')
    
    print("\n✓ 完成！")


def example_starsolo():
    """读取 STARsolo 输出"""
    
    print("\n" + "=" * 70)
    print("示例4：读取 STARsolo 输出")
    print("=" * 70)
    
    reader = UniversalScRNAReader(verbose=True)
    
    # STARsolo 输出（注意：未压缩）
    adata = reader.read_10x_mtx(
        'Solo.out/Gene/filtered/',
        var_names='gene_symbols',
        compressed=False  # 关键：STARsolo 不压缩
    )
    
    reader.save_h5ad(adata, 'starsolo_output.h5ad')
    
    print("\n✓ 完成！")


def example_csv():
    """读取 CSV 矩阵"""
    
    print("\n" + "=" * 70)
    print("示例5：读取 CSV 矩阵")
    print("=" * 70)
    
    reader = UniversalScRNAReader(verbose=True)
    
    # CSV 文件（行是基因，列是细胞）
    adata = reader.read_csv(
        'expression_matrix.csv',
        first_column_names=True,
        transpose=True  # 转置为 cells × genes
    )
    
    reader.save_h5ad(adata, 'from_csv.h5ad')
    
    print("\n✓ 完成！")


def example_one_liner():
    """一行代码处理"""
    
    print("\n" + "=" * 70)
    print("示例6：一行代码处理（推荐）")
    print("=" * 70)
    
    reader = UniversalScRNAReader(verbose=True)
    
    # 一键处理：读取 → 验证 → 标准化 → 保存
    adata = reader.process_single_sample(
        input_path='filtered_feature_bc_matrix/',
        output_path='sample1.h5ad',
        sample_id='sample1',
        var_names='gene_symbols',
        cache=True
    )
    
    print(f"\n最终数据: {adata}")
    print("✓ 完成！")


def example_auto_detect():
    """自动检测格式"""
    
    print("\n" + "=" * 70)
    print("示例7：自动检测格式（最简单）")
    print("=" * 70)
    
    reader = UniversalScRNAReader(verbose=True)
    
    # 支持的所有格式都可以自动检测
    test_files = [
        'filtered_feature_bc_matrix/',  # 10X MTX
        'filtered_feature_bc_matrix.h5', # 10X H5
        'data.h5ad',                      # H5AD
        'data.loom',                      # Loom
        'data.zarr',                      # Zarr
        'matrix.csv',                     # CSV
        'matrix.tsv',                     # TSV
        'matrix.xlsx',                    # Excel
    ]
    
    for file_path in test_files:
        print(f"\n处理: {file_path}")
        try:
            # 自动检测并读取
            adata = reader.read_auto(file_path)
            print(f"✓ 成功读取: {adata.n_obs} cells × {adata.n_vars} genes")
        except FileNotFoundError:
            print(f"⊘ 文件不存在（仅用于演示）")
        except Exception as e:
            print(f"✗ 错误: {e}")


if __name__ == '__main__':
    print("=" * 70)
    print("scRNA-DataHub 基础使用示例")
    print("=" * 70)
    print("\n本文件包含 7 个示例，演示不同的使用场景")
    print("注意：运行示例需要有相应的数据文件\n")
    
    # 取消注释以运行不同的示例
    
    # example_basic_usage()        # 示例1：基本使用
    # example_10x_genomics()       # 示例2：10X Genomics
    # example_dnb_c4()             # 示例3：DNB C4
    # example_starsolo()           # 示例4：STARsolo
    # example_csv()                # 示例5：CSV
    # example_one_liner()          # 示例6：一行代码（推荐）
    example_auto_detect()          # 示例7：自动检测（最简单）
    
    print("\n" + "=" * 70)
    print("更多示例请查看:")
    print("  - batch_processing.py  （批量处理）")
    print("  - multi_sample.py      （多样本整合）")
    print("=" * 70)
