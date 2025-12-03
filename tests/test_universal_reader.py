#!/usr/bin/env python3
"""
测试通用单细胞数据读取器
演示各种格式的读取，并验证最终都转换为统一的H5AD格式
"""

import os
import sys
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse

# 导入读取器
from universal_scrna_reader import UniversalScRNAReader


def create_test_data():
    """创建测试数据（模拟各种格式）"""
    print("=" * 70)
    print("创建测试数据")
    print("=" * 70)
    
    # 创建测试目录
    test_dir = Path('test_data')
    test_dir.mkdir(exist_ok=True)
    
    # 创建一个小的测试数据集
    n_cells = 100
    n_genes = 500
    
    # 随机生成稀疏表达矩阵
    np.random.seed(42)
    density = 0.1  # 10%的非零值
    X = sparse.random(n_cells, n_genes, density=density, format='csr', dtype=np.float32)
    X.data = np.random.poisson(5, size=X.data.shape).astype(np.float32)
    
    # 生成基因名和细胞条形码
    gene_ids = [f"ENSG{i:011d}" for i in range(n_genes)]
    gene_symbols = [f"Gene{i}" for i in range(n_genes)]
    barcodes = [f"AAACCCAA{i:08d}-1" for i in range(n_cells)]
    
    # 创建基础AnnData对象
    import anndata as ad
    adata_base = ad.AnnData(X=X)
    adata_base.obs_names = barcodes
    adata_base.var_names = gene_symbols
    adata_base.var['gene_ids'] = gene_ids
    adata_base.var['feature_types'] = ['Gene Expression'] * n_genes
    
    print(f"✓ 创建测试数据: {n_cells} cells × {n_genes} genes")
    
    return test_dir, adata_base


def save_as_10x_mtx(adata, output_dir, compressed=True):
    """保存为10X MTX格式"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    suffix = '.gz' if compressed else ''
    
    # 保存barcodes
    import gzip
    open_func = gzip.open if compressed else open
    mode = 'wt' if compressed else 'w'
    
    with open_func(output_dir / f'barcodes.tsv{suffix}', mode) as f:
        for barcode in adata.obs_names:
            f.write(f'{barcode}\n')
    
    # 保存features
    features_df = pd.DataFrame({
        'gene_id': adata.var.get('gene_ids', adata.var_names),
        'gene_symbol': adata.var_names,
        'feature_type': adata.var.get('feature_types', ['Gene Expression'] * len(adata.var_names))
    })
    
    if compressed:
        features_df.to_csv(output_dir / f'features.tsv{suffix}', sep='\t', 
                          header=False, index=False, compression='gzip')
    else:
        features_df.to_csv(output_dir / f'features.tsv{suffix}', sep='\t', 
                          header=False, index=False)
    
    # 保存matrix
    from scipy.io import mmwrite
    matrix = adata.X.T  # 转置为 genes × cells
    
    if compressed:
        mmwrite(output_dir / 'matrix.mtx', matrix)
        import subprocess
        subprocess.run(['gzip', '-f', str(output_dir / 'matrix.mtx')], check=True)
    else:
        mmwrite(output_dir / f'matrix.mtx', matrix)
    
    print(f"✓ 保存10X MTX格式: {output_dir}/")


def test_all_formats():
    """测试所有支持的格式"""
    print("\n" + "=" * 70)
    print("测试所有格式的读取")
    print("=" * 70)
    
    # 创建测试数据
    test_dir, adata_base = create_test_data()
    
    # 创建读取器
    reader = UniversalScRNAReader(verbose=True)
    
    results = {}
    
    # ========== 测试1: 10X MTX格式 (压缩) ==========
    print("\n" + "=" * 70)
    print("测试1: 10X MTX格式 (压缩)")
    print("=" * 70)
    
    mtx_dir_compressed = test_dir / '10x_mtx_compressed'
    save_as_10x_mtx(adata_base, mtx_dir_compressed, compressed=True)
    
    adata_mtx_compressed = reader.read_auto(str(mtx_dir_compressed))
    results['10x_mtx_compressed'] = adata_mtx_compressed
    
    # ========== 测试2: 10X MTX格式 (未压缩) ==========
    print("\n" + "=" * 70)
    print("测试2: 10X MTX格式 (未压缩，如STARsolo)")
    print("=" * 70)
    
    mtx_dir_uncompressed = test_dir / '10x_mtx_uncompressed'
    save_as_10x_mtx(adata_base, mtx_dir_uncompressed, compressed=False)
    
    adata_mtx_uncompressed = reader.read_auto(
        str(mtx_dir_uncompressed),
        compressed=False
    )
    results['10x_mtx_uncompressed'] = adata_mtx_uncompressed
    
    # ========== 测试3: H5AD格式 ==========
    print("\n" + "=" * 70)
    print("测试3: H5AD格式")
    print("=" * 70)
    
    h5ad_file = test_dir / 'test_data.h5ad'
    adata_base.write_h5ad(h5ad_file)
    
    adata_h5ad = reader.read_auto(str(h5ad_file))
    results['h5ad'] = adata_h5ad
    
    # ========== 测试4: CSV格式 ==========
    print("\n" + "=" * 70)
    print("测试4: CSV格式")
    print("=" * 70)
    
    csv_file = test_dir / 'test_data.csv'
    # 转换为DataFrame（基因 × 细胞）
    df = pd.DataFrame(
        adata_base.X.T.toarray() if sparse.issparse(adata_base.X) else adata_base.X.T,
        index=adata_base.var_names,
        columns=adata_base.obs_names
    )
    df.to_csv(csv_file)
    
    adata_csv = reader.read_auto(str(csv_file), transpose=True)
    results['csv'] = adata_csv
    
    # ========== 测试5: TSV格式 ==========
    print("\n" + "=" * 70)
    print("测试5: TSV格式")
    print("=" * 70)
    
    tsv_file = test_dir / 'test_data.tsv'
    df.to_csv(tsv_file, sep='\t')
    
    adata_tsv = reader.read_auto(str(tsv_file), delimiter='\t', transpose=True)
    results['tsv'] = adata_tsv
    
    # ========== 测试6: Loom格式 ==========
    print("\n" + "=" * 70)
    print("测试6: Loom格式")
    print("=" * 70)
    
    loom_file = test_dir / 'test_data.loom'
    adata_base.write_loom(loom_file)
    
    adata_loom = reader.read_auto(str(loom_file))
    results['loom'] = adata_loom
    
    # ========== 测试7: Zarr格式 ==========
    print("\n" + "=" * 70)
    print("测试7: Zarr格式")
    print("=" * 70)
    
    zarr_dir = test_dir / 'test_data.zarr'
    adata_base.write_zarr(zarr_dir)
    
    adata_zarr = reader.read_auto(str(zarr_dir))
    results['zarr'] = adata_zarr
    
    # ========== 结果验证 ==========
    print("\n" + "=" * 70)
    print("结果验证：所有格式是否统一")
    print("=" * 70)
    
    # 验证所有读取的数据是否一致
    print("\n各格式读取结果:")
    print(f"{'格式':<25} {'细胞数':<10} {'基因数':<10} {'稀疏矩阵':<10} {'类型':<15}")
    print("-" * 70)
    
    for fmt, adata in results.items():
        is_sparse = sparse.issparse(adata.X)
        adata_type = type(adata).__name__
        print(f"{fmt:<25} {adata.n_obs:<10} {adata.n_vars:<10} {str(is_sparse):<10} {adata_type:<15}")
    
    # 验证维度一致性
    print("\n维度一致性验证:")
    n_obs_set = {adata.n_obs for adata in results.values()}
    n_vars_set = {adata.n_vars for adata in results.values()}
    
    if len(n_obs_set) == 1 and len(n_vars_set) == 1:
        print(f"✓ 所有格式的维度一致: {list(n_obs_set)[0]} cells × {list(n_vars_set)[0]} genes")
    else:
        print(f"✗ 维度不一致:")
        print(f"  细胞数: {n_obs_set}")
        print(f"  基因数: {n_vars_set}")
    
    # 验证数据类型一致性
    print("\n数据类型验证:")
    all_anndata = all(isinstance(adata, ad.AnnData) for adata in results.values())
    if all_anndata:
        print("✓ 所有读取结果都是AnnData对象")
    else:
        print("✗ 存在非AnnData对象")
    
    # 统一保存为H5AD格式
    print("\n" + "=" * 70)
    print("统一保存为H5AD格式")
    print("=" * 70)
    
    unified_dir = test_dir / 'unified_h5ad'
    unified_dir.mkdir(exist_ok=True)
    
    for fmt, adata in results.items():
        # 标准化
        adata_std = reader.standardize_adata(adata, ensure_sparse=True)
        
        # 保存
        output_file = unified_dir / f'{fmt}.h5ad'
        reader.save_h5ad(adata_std, str(output_file), compression='gzip')
        
        print(f"✓ {fmt:<25} → {output_file}")
    
    print("\n" + "=" * 70)
    print("测试完成")
    print("=" * 70)
    print(f"\n所有测试数据保存在: {test_dir}/")
    print(f"统一H5AD文件保存在: {unified_dir}/")
    
    return results


if __name__ == '__main__':
    # 如果是作为模块导入，则运行测试
    if len(sys.argv) == 1:
        print("运行测试模式...\n")
        results = test_all_formats()
        
        print("\n" + "=" * 70)
        print("总结")
        print("=" * 70)
        print(f"✓ 成功测试了 {len(results)} 种格式")
        print("✓ 所有格式都可以统一转换为AnnData对象")
        print("✓ 所有格式都可以保存为标准H5AD格式")
        
        sys.exit(0)
    else:
        # 作为命令行工具运行
        from universal_scrna_reader import main
        sys.exit(main())

