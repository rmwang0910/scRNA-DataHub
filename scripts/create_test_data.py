#!/usr/bin/env python3
"""
创建所有格式的测试数据
基于scanpy现有测试数据生成其他格式
"""

import os
import sys
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.io as sio


def create_all_test_data(output_dir='test_data_all_formats'):
    """创建所有格式的测试数据"""
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    print("=" * 70)
    print("创建所有格式的测试数据")
    print("=" * 70)
    
    # ========== 获取基础数据 ==========
    print("\n1. 读取基础数据...")
    
    # 使用scanpy自带的PBMC数据（较小，适合测试）
    adata_base = sc.datasets.pbmc68k_reduced()
    print(f"   基础数据: {adata_base.n_obs} cells × {adata_base.n_vars} genes")
    
    # 为了测试方便，取一个子集
    np.random.seed(42)
    cell_indices = np.random.choice(adata_base.n_obs, size=min(1000, adata_base.n_obs), replace=False)
    gene_indices = np.random.choice(adata_base.n_vars, size=min(500, adata_base.n_vars), replace=False)
    
    adata_subset = adata_base[cell_indices, :][:, gene_indices].copy()
    print(f"   测试子集: {adata_subset.n_obs} cells × {adata_subset.n_vars} genes")
    
    results = {}
    
    # ========== 1-2. 10X格式（使用scanpy现有数据）==========
    print("\n2. 10X格式 - 使用scanpy测试数据")
    
    scanpy_test_data = Path('/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x_data')
    
    # 10X MTX v3 (压缩)
    test_10x_v3 = scanpy_test_data / '3.0.0/filtered_feature_bc_matrix'
    if test_10x_v3.exists():
        print(f"   ✓ 10X MTX v3 (压缩): {test_10x_v3}")
        results['10x_mtx_v3_compressed'] = test_10x_v3
    
    # 10X MTX v2 (未压缩)
    test_10x_v2 = scanpy_test_data / '1.2.0/filtered_gene_bc_matrices/hg19_chr21'
    if test_10x_v2.exists():
        print(f"   ✓ 10X MTX v2 (未压缩): {test_10x_v2}")
        results['10x_mtx_v2_uncompressed'] = test_10x_v2
    
    # 10X H5 v3
    test_10x_h5_v3 = scanpy_test_data / '3.0.0/filtered_feature_bc_matrix.h5'
    if test_10x_h5_v3.exists():
        print(f"   ✓ 10X H5 v3: {test_10x_h5_v3}")
        results['10x_h5_v3'] = test_10x_h5_v3
    
    # ========== 3. H5AD格式（使用scanpy现有数据）==========
    print("\n3. H5AD格式 - 使用scanpy自带数据")
    
    scanpy_h5ad = Path('/Users/warm/华大智造/TCGA/gdc/scanpy/src/scanpy/datasets/10x_pbmc68k_reduced.h5ad')
    if scanpy_h5ad.exists():
        print(f"   ✓ H5AD: {scanpy_h5ad}")
        results['h5ad'] = scanpy_h5ad
    
    # ========== 4. Loom格式 ==========
    print("\n4. 创建Loom格式...")
    
    loom_file = output_dir / 'test_data.loom'
    try:
        adata_subset.write_loom(loom_file)
        print(f"   ✓ Loom: {loom_file}")
        results['loom'] = loom_file
    except Exception as e:
        print(f"   ✗ Loom创建失败: {e}")
    
    # ========== 5. Zarr格式（使用scanpy现有数据）==========
    print("\n5. Zarr格式 - 使用scanpy测试数据")
    
    scanpy_zarr = Path('/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x-10k-subset.zarr')
    if scanpy_zarr.exists():
        print(f"   ✓ Zarr: {scanpy_zarr}")
        results['zarr'] = scanpy_zarr
    else:
        # 如果不存在，创建新的
        zarr_dir = output_dir / 'test_data.zarr'
        try:
            adata_subset.write_zarr(zarr_dir)
            print(f"   ✓ Zarr: {zarr_dir}")
            results['zarr'] = zarr_dir
        except Exception as e:
            print(f"   ✗ Zarr创建失败: {e}")
    
    # ========== 6. CSV格式 ==========
    print("\n6. 创建CSV格式...")
    
    csv_file = output_dir / 'test_expression.csv'
    
    # 转换为DataFrame（基因×细胞）
    df = pd.DataFrame(
        adata_subset.X.T.toarray() if sp.issparse(adata_subset.X) else adata_subset.X.T,
        index=adata_subset.var_names,
        columns=adata_subset.obs_names
    )
    df.to_csv(csv_file)
    print(f"   ✓ CSV: {csv_file}")
    print(f"     格式: 基因(行) × 细胞(列)")
    results['csv'] = csv_file
    
    # 压缩版本
    csv_gz_file = output_dir / 'test_expression.csv.gz'
    df.to_csv(csv_gz_file, compression='gzip')
    print(f"   ✓ CSV.GZ: {csv_gz_file}")
    results['csv_gz'] = csv_gz_file
    
    # ========== 7. TSV格式 ==========
    print("\n7. 创建TSV格式...")
    
    tsv_file = output_dir / 'test_expression.tsv'
    df.to_csv(tsv_file, sep='\t')
    print(f"   ✓ TSV: {tsv_file}")
    results['tsv'] = tsv_file
    
    # 压缩版本
    tsv_gz_file = output_dir / 'test_expression.tsv.gz'
    df.to_csv(tsv_gz_file, sep='\t', compression='gzip')
    print(f"   ✓ TSV.GZ: {tsv_gz_file}")
    results['tsv_gz'] = tsv_gz_file
    
    # ========== 8. Excel格式 ==========
    print("\n8. 创建Excel格式...")
    
    excel_file = output_dir / 'test_expression.xlsx'
    try:
        # 取一个更小的子集（Excel有限制）
        df_small = df.iloc[:1000, :50] if df.shape[0] > 1000 else df
        df_small.to_excel(excel_file, sheet_name='expression_matrix')
        print(f"   ✓ Excel: {excel_file}")
        print(f"     注意: 只保留了 {df_small.shape[0]} genes × {df_small.shape[1]} cells")
        results['excel'] = excel_file
    except Exception as e:
        print(f"   ✗ Excel创建失败: {e}")
        print(f"     提示: 可能需要安装 openpyxl: pip install openpyxl")
    
    # ========== 9. MTX格式（单文件）==========
    print("\n9. 创建MTX格式（单文件）...")
    
    mtx_file = output_dir / 'test_matrix.mtx'
    
    # 转置为 基因×细胞
    matrix = adata_subset.X.T
    sio.mmwrite(mtx_file, matrix)
    print(f"   ✓ MTX: {mtx_file}")
    results['mtx'] = mtx_file
    
    # 压缩版本
    import subprocess
    import shutil
    if shutil.which('gzip'):
        mtx_gz_file = output_dir / 'test_matrix.mtx.gz'
        shutil.copy(mtx_file, output_dir / 'test_matrix_temp.mtx')
        subprocess.run(['gzip', '-f', str(output_dir / 'test_matrix_temp.mtx')])
        shutil.move(output_dir / 'test_matrix_temp.mtx.gz', mtx_gz_file)
        print(f"   ✓ MTX.GZ: {mtx_gz_file}")
        results['mtx_gz'] = mtx_gz_file
    
    # ========== 10. HDF5格式 ==========
    print("\n10. 创建HDF5格式...")
    
    hdf5_file = output_dir / 'test_data.hdf5'
    
    import h5py
    with h5py.File(hdf5_file, 'w') as f:
        # 创建主数据集
        data_group = f.create_group('data')
        data_group.create_dataset(
            'expression',
            data=adata_subset.X.toarray() if sp.issparse(adata_subset.X) else adata_subset.X
        )
        data_group.create_dataset('gene_names', data=adata_subset.var_names.values.astype('S'))
        data_group.create_dataset('cell_names', data=adata_subset.obs_names.values.astype('S'))
        
        # 添加属性
        f.attrs['n_obs'] = adata_subset.n_obs
        f.attrs['n_vars'] = adata_subset.n_vars
    
    print(f"   ✓ HDF5: {hdf5_file}")
    print(f"     读取方式: sc.read_hdf('{hdf5_file}', key='data')")
    results['hdf5'] = hdf5_file
    
    # ========== 11. SOFT.GZ格式（使用scanpy下载）==========
    print("\n11. SOFT.GZ格式 - 使用scanpy数据集")
    
    try:
        # 下载Burczynski06数据（SOFT格式）
        print("   下载GEO数据（需要网络连接）...")
        adata_soft = sc.datasets.burczynski06()
        soft_file = Path.home() / '.cache/scanpy-data/burczynski06/GDS1615_full.soft.gz'
        if soft_file.exists():
            print(f"   ✓ SOFT.GZ: {soft_file}")
            results['soft_gz'] = soft_file
        else:
            print(f"   ⚠️  SOFT.GZ下载失败，请检查网络连接")
    except Exception as e:
        print(f"   ✗ SOFT.GZ下载失败: {e}")
        print(f"     提示: 需要网络连接从GEO下载")
    
    # ========== 12. UMI-tools格式 ==========
    print("\n12. 创建UMI-tools格式...")
    
    umi_tools_file = output_dir / 'umi_tools_counts.tsv.gz'
    
    # UMI-tools格式：基因×细胞计数表
    df_umi = pd.DataFrame(
        adata_subset.X.T.toarray() if sp.issparse(adata_subset.X) else adata_subset.X.T,
        index=adata_subset.var_names,
        columns=adata_subset.obs_names
    )
    
    df_umi.to_csv(umi_tools_file, sep='\t', compression='gzip')
    print(f"   ✓ UMI-tools: {umi_tools_file}")
    results['umi_tools'] = umi_tools_file
    
    # ========== 额外：创建自己的10X格式数据 ==========
    print("\n13. 创建自定义10X MTX格式...")
    
    custom_10x_dir = output_dir / 'custom_10x_mtx'
    custom_10x_dir.mkdir(exist_ok=True)
    
    # barcodes
    import gzip
    with gzip.open(custom_10x_dir / 'barcodes.tsv.gz', 'wt') as f:
        for bc in adata_subset.obs_names:
            f.write(f'{bc}\n')
    
    # features
    features_df = pd.DataFrame({
        'gene_id': adata_subset.var.get('gene_ids', adata_subset.var_names),
        'gene_symbol': adata_subset.var_names,
        'feature_type': ['Gene Expression'] * adata_subset.n_vars
    })
    features_df.to_csv(custom_10x_dir / 'features.tsv.gz', sep='\t',
                       header=False, index=False, compression='gzip')
    
    # matrix
    matrix = adata_subset.X.T  # 转置为 基因×细胞
    sio.mmwrite(custom_10x_dir / 'matrix.mtx', matrix)
    
    # 压缩matrix
    import subprocess
    subprocess.run(['gzip', '-f', str(custom_10x_dir / 'matrix.mtx')])
    
    print(f"   ✓ 自定义10X MTX: {custom_10x_dir}/")
    results['custom_10x_mtx'] = custom_10x_dir
    
    # ========== 生成测试数据清单 ==========
    print("\n" + "=" * 70)
    print("测试数据清单")
    print("=" * 70)
    
    manifest = []
    
    # Scanpy现有数据
    print("\n✅ Scanpy现有测试数据:")
    print(f"  1. 10X MTX v3 (压缩):  scanpy/tests/_data/10x_data/3.0.0/filtered_feature_bc_matrix/")
    print(f"  2. 10X MTX v2 (未压缩): scanpy/tests/_data/10x_data/1.2.0/filtered_gene_bc_matrices/hg19_chr21/")
    print(f"  3. 10X H5 v3:          scanpy/tests/_data/10x_data/3.0.0/filtered_feature_bc_matrix.h5")
    print(f"  4. 10X H5 v2:          scanpy/tests/_data/10x_data/1.2.0/filtered_gene_bc_matrices_h5.h5")
    print(f"  5. H5AD:               scanpy/src/scanpy/datasets/10x_pbmc68k_reduced.h5ad")
    print(f"  6. TXT:                scanpy/src/scanpy/datasets/krumsiek11.txt")
    print(f"  7. Zarr:               scanpy/tests/_data/10x-10k-subset.zarr/")
    
    manifest.extend([
        "# Scanpy现有数据",
        "10x_mtx_v3=/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x_data/3.0.0/filtered_feature_bc_matrix/",
        "10x_mtx_v2=/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x_data/1.2.0/filtered_gene_bc_matrices/hg19_chr21/",
        "10x_h5_v3=/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x_data/3.0.0/filtered_feature_bc_matrix.h5",
        "h5ad=/Users/warm/华大智造/TCGA/gdc/scanpy/src/scanpy/datasets/10x_pbmc68k_reduced.h5ad",
        "txt=/Users/warm/华大智造/TCGA/gdc/scanpy/src/scanpy/datasets/krumsiek11.txt",
        "zarr=/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x-10k-subset.zarr/",
        ""
    ])
    
    # 新创建的数据
    print("\n📦 新创建的测试数据:")
    new_data_count = 0
    for name, path in results.items():
        if 'scanpy' not in str(path):
            new_data_count += 1
            print(f"  {new_data_count}. {name}: {path}")
            manifest.append(f"{name}={path}")
    
    # 保存清单文件
    manifest_file = output_dir / 'test_data_manifest.txt'
    with open(manifest_file, 'w') as f:
        f.write('\n'.join(manifest))
    
    print(f"\n✓ 数据清单已保存: {manifest_file}")
    
    # ========== 生成测试脚本 ==========
    print("\n生成测试脚本...")
    
    test_script = output_dir / 'run_all_format_tests.sh'
    
    with open(test_script, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("# 测试所有格式的读取\n\n")
        f.write("echo '测试所有单细胞数据格式'\n")
        f.write("echo '=' '=' '=' '=' '=' '=' '=' '=' '=' '=' '=' '=' '=' '=' '=' '=' '=' '='\n\n")
        
        test_cases = [
            ("10X MTX v3 (压缩)", "/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x_data/3.0.0/filtered_feature_bc_matrix/", ""),
            ("10X MTX v2 (未压缩)", "/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x_data/1.2.0/filtered_gene_bc_matrices/hg19_chr21/", "--no-compressed"),
            ("10X H5 v3", "/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x_data/3.0.0/filtered_feature_bc_matrix.h5", ""),
            ("H5AD", "/Users/warm/华大智造/TCGA/gdc/scanpy/src/scanpy/datasets/10x_pbmc68k_reduced.h5ad", ""),
            ("Loom", str(results.get('loom', '')), "") if 'loom' in results else None,
            ("Zarr", "/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x-10k-subset.zarr/", ""),
            ("CSV", str(results.get('csv', '')), "--transpose") if 'csv' in results else None,
            ("TSV", str(results.get('tsv', '')), "--delimiter '\\t' --transpose") if 'tsv' in results else None,
            ("Excel", str(results.get('excel', '')), "--sheet expression_matrix") if 'excel' in results else None,
            ("MTX", str(results.get('mtx', '')), "") if 'mtx' in results else None,
            ("UMI-tools", str(results.get('umi_tools', '')), "") if 'umi_tools' in results else None,
        ]
        
        for i, test_case in enumerate(test_cases, 1):
            if test_case is None:
                continue
            
            name, path, args = test_case
            if not path or path == '':
                continue
                
            f.write(f"echo '{i}. 测试{name}'\n")
            f.write(f"python ../universal_scrna_reader.py \\\n")
            f.write(f"  '{path}' \\\n")
            f.write(f"  -o test_{name.replace(' ', '_').replace('(', '').replace(')', '').lower()}.h5ad")
            if args:
                f.write(f" \\\n  {args}")
            f.write("\n\n")
        
        f.write("echo '所有测试完成！'\n")
    
    # 设置执行权限
    os.chmod(test_script, 0o755)
    print(f"✓ 测试脚本已生成: {test_script}")
    
    # ========== 总结 ==========
    print("\n" + "=" * 70)
    print("总结")
    print("=" * 70)
    print(f"\n测试数据目录: {output_dir.absolute()}/")
    print(f"数据清单文件: {manifest_file}")
    print(f"测试脚本: {test_script}")
    
    print(f"\n创建的格式数量: {len(results)}")
    for name in results:
        print(f"  ✓ {name}")
    
    print("\n下一步:")
    print(f"  1. 查看数据清单: cat {manifest_file}")
    print(f"  2. 运行测试脚本: bash {test_script}")
    print(f"  3. 或手动测试: python universal_scrna_reader.py <数据路径> -o output.h5ad")
    
    return results


if __name__ == '__main__':
    results = create_all_test_data()
    
    print("\n" + "=" * 70)
    print("验证数据读取")
    print("=" * 70)
    
    # 验证几个关键格式
    print("\n快速验证...")
    
    # 验证10X MTX v3
    scanpy_10x_v3 = Path('/Users/warm/华大智造/TCGA/gdc/scanpy/tests/_data/10x_data/3.0.0/filtered_feature_bc_matrix')
    if scanpy_10x_v3.exists():
        adata = sc.read_10x_mtx(scanpy_10x_v3)
        print(f"✓ 10X MTX v3: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # 验证H5AD
    scanpy_h5ad = Path('/Users/warm/华大智造/TCGA/gdc/scanpy/src/scanpy/datasets/10x_pbmc68k_reduced.h5ad')
    if scanpy_h5ad.exists():
        adata = sc.read_h5ad(scanpy_h5ad)
        print(f"✓ H5AD: {adata.n_obs} cells × {adata.n_vars} genes")
    
    print("\n✓ 所有测试数据准备完成！")

