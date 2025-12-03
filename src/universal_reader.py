#!/usr/bin/env python3
"""
通用单细胞RNA-seq数据读取脚本
支持所有常见的单细胞数据格式，统一转换为H5AD格式

支持的格式：
- 10X Genomics (MTX格式和H5格式)
- H5AD (AnnData格式)
- Loom
- CSV/TSV
- Excel
- Zarr
- UMI-tools
- SOFT.GZ (GEO)
- RDS/H5Seurat (需要R转换)
"""

import os
import sys
from pathlib import Path
from typing import Optional, Literal
import warnings

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp


class UniversalScRNAReader:
    """通用单细胞RNA-seq数据读取器"""
    
    def __init__(self, verbose: bool = True):
        """
        初始化读取器
        
        参数:
            verbose: 是否打印详细信息
        """
        self.verbose = verbose
        self.supported_formats = {
            '10x_mtx': '10X Genomics MTX格式 (3个文件)',
            '10x_h5': '10X Genomics H5格式',
            'h5ad': 'AnnData H5AD格式',
            'loom': 'Loom格式',
            'csv': 'CSV格式',
            'tsv': 'TSV/TXT格式',
            'excel': 'Excel格式',
            'zarr': 'Zarr格式',
            'umi_tools': 'UMI-tools格式',
            'soft_gz': 'GEO SOFT.GZ格式',
            'mtx': '单个MTX文件',
            'hdf5': 'HDF5格式'
        }
    
    def detect_format(self, input_path: str) -> str:
        """
        自动检测数据格式
        
        参数:
            input_path: 输入文件或目录路径
            
        返回:
            检测到的格式类型
        """
        path = Path(input_path)
        
        # 检查是否是目录（可能是10X格式或Zarr）
        if path.is_dir():
            files = os.listdir(path)
            
            # 检查是否是10X MTX格式
            has_matrix = any(f in files for f in ['matrix.mtx', 'matrix.mtx.gz'])
            has_barcodes = any(f in files for f in ['barcodes.tsv', 'barcodes.tsv.gz'])
            has_features = any(f in files for f in ['features.tsv', 'features.tsv.gz', 'genes.tsv', 'genes.tsv.gz'])
            
            if has_matrix and has_barcodes and has_features:
                return '10x_mtx'
            
            # 检查是否是Zarr格式
            if '.zattrs' in files or '.zgroup' in files:
                return 'zarr'
            
            return 'unknown_directory'
        
        # 检查文件扩展名
        ext = ''.join(path.suffixes).lower()
        
        if ext in ['.h5ad']:
            return 'h5ad'
        elif ext in ['.h5', '.hdf5']:
            # 尝试判断是10X H5还是普通HDF5
            try:
                import h5py
                with h5py.File(path, 'r') as f:
                    if 'matrix' in f:
                        return '10x_h5'
                    else:
                        return 'hdf5'
            except:
                return 'hdf5'
        elif ext in ['.loom']:
            return 'loom'
        elif ext in ['.csv', '.csv.gz', '.csv.bz2']:
            return 'csv'
        elif ext in ['.tsv', '.tsv.gz', '.tsv.bz2', '.txt', '.txt.gz', '.tab', '.tab.gz']:
            return 'tsv'
        elif ext in ['.xlsx', '.xls']:
            return 'excel'
        elif ext in ['.mtx', '.mtx.gz']:
            return 'mtx'
        elif ext == '.soft.gz':
            return 'soft_gz'
        elif ext in ['.rds']:
            return 'rds'
        elif ext in ['.h5seurat']:
            return 'h5seurat'
        else:
            return 'unknown'
    
    def read_10x_mtx(
        self,
        path: str,
        var_names: Literal['gene_symbols', 'gene_ids'] = 'gene_symbols',
        make_unique: bool = True,
        cache: bool = True,
        gex_only: bool = True,
        compressed: bool = True
    ) -> ad.AnnData:
        """
        读取10X MTX格式
        
        参数:
            path: 包含matrix.mtx, barcodes.tsv, features.tsv的目录
            var_names: 使用基因符号或基因ID
            make_unique: 是否使基因名唯一
            cache: 是否使用缓存
            gex_only: 是否只保留基因表达数据
            compressed: 是否期望文件被压缩（Cell Ranger v3+: True, STARsolo: False）
        """
        if self.verbose:
            print(f"读取10X MTX格式: {path}")
            print(f"  - var_names: {var_names}")
            print(f"  - compressed: {compressed}")
        
        adata = sc.read_10x_mtx(
            path,
            var_names=var_names,
            make_unique=make_unique,
            cache=cache,
            gex_only=gex_only,
            compressed=compressed
        )
        
        if self.verbose:
            print(f"✓ 读取成功: {adata.n_obs} cells × {adata.n_vars} genes")
        
        return adata
    
    def read_10x_h5(
        self,
        path: str,
        genome: Optional[str] = None,
        gex_only: bool = True
    ) -> ad.AnnData:
        """
        读取10X H5格式
        
        参数:
            path: H5文件路径
            genome: 基因组名称（多基因组时需要指定）
            gex_only: 是否只保留基因表达数据
        """
        if self.verbose:
            print(f"读取10X H5格式: {path}")
            if genome:
                print(f"  - genome: {genome}")
        
        adata = sc.read_10x_h5(path, genome=genome, gex_only=gex_only)
        
        if self.verbose:
            print(f"✓ 读取成功: {adata.n_obs} cells × {adata.n_vars} genes")
        
        return adata
    
    def read_h5ad(self, path: str, backed: Optional[str] = None) -> ad.AnnData:
        """
        读取H5AD格式
        
        参数:
            path: H5AD文件路径
            backed: 'r'(只读) 或 'r+'(读写) 启用backed模式，None为完全加载
        """
        if self.verbose:
            print(f"读取H5AD格式: {path}")
            if backed:
                print(f"  - backed模式: {backed}")
        
        adata = sc.read_h5ad(path, backed=backed)
        
        if self.verbose:
            if backed:
                print(f"✓ 读取成功 (backed模式): {adata.n_obs} cells × {adata.n_vars} genes")
            else:
                print(f"✓ 读取成功: {adata.n_obs} cells × {adata.n_vars} genes")
        
        return adata
    
    def read_loom(self, path: str) -> ad.AnnData:
        """读取Loom格式"""
        if self.verbose:
            print(f"读取Loom格式: {path}")
        
        adata = sc.read_loom(path)
        
        if self.verbose:
            print(f"✓ 读取成功: {adata.n_obs} cells × {adata.n_vars} genes")
        
        return adata
    
    def read_csv(
        self,
        path: str,
        first_column_names: bool = True,
        transpose: bool = False
    ) -> ad.AnnData:
        """
        读取CSV格式
        
        参数:
            path: CSV文件路径
            first_column_names: 第一列是否是行名
            transpose: 是否转置（如果行是细胞，列是基因）
        """
        if self.verbose:
            print(f"读取CSV格式: {path}")
        
        adata = sc.read_csv(path, first_column_names=first_column_names)
        
        if transpose:
            if self.verbose:
                print("  - 转置矩阵")
            adata = adata.T
        
        if self.verbose:
            print(f"✓ 读取成功: {adata.n_obs} cells × {adata.n_vars} genes")
        
        return adata
    
    def read_text(
        self,
        path: str,
        delimiter: str = '\t',
        first_column_names: bool = True,
        transpose: bool = False
    ) -> ad.AnnData:
        """
        读取文本格式（TSV/TXT等）
        
        参数:
            path: 文件路径
            delimiter: 分隔符
            first_column_names: 第一列是否是行名
            transpose: 是否转置
        """
        if self.verbose:
            print(f"读取文本格式: {path}")
            print(f"  - delimiter: '{delimiter}'")
        
        adata = sc.read_text(path, delimiter=delimiter, first_column_names=first_column_names)
        
        if transpose:
            if self.verbose:
                print("  - 转置矩阵")
            adata = adata.T
        
        if self.verbose:
            print(f"✓ 读取成功: {adata.n_obs} cells × {adata.n_vars} genes")
        
        return adata
    
    def read_excel(
        self,
        path: str,
        sheet: str = 'Sheet1',
        transpose: bool = False
    ) -> ad.AnnData:
        """
        读取Excel格式
        
        参数:
            path: Excel文件路径
            sheet: Sheet名称
            transpose: 是否转置
        """
        if self.verbose:
            print(f"读取Excel格式: {path}")
            print(f"  - sheet: {sheet}")
        
        adata = sc.read_excel(path, sheet=sheet)
        
        if transpose:
            if self.verbose:
                print("  - 转置矩阵")
            adata = adata.T
        
        if self.verbose:
            print(f"✓ 读取成功: {adata.n_obs} cells × {adata.n_vars} genes")
        
        return adata
    
    def read_zarr(self, path: str) -> ad.AnnData:
        """读取Zarr格式"""
        if self.verbose:
            print(f"读取Zarr格式: {path}")
        
        adata = sc.read_zarr(path)
        
        if self.verbose:
            print(f"✓ 读取成功: {adata.n_obs} cells × {adata.n_vars} genes")
        
        return adata
    
    def read_umi_tools(self, path: str) -> ad.AnnData:
        """读取UMI-tools格式"""
        if self.verbose:
            print(f"读取UMI-tools格式: {path}")
        
        adata = sc.read_umi_tools(path)
        
        if self.verbose:
            print(f"✓ 读取成功: {adata.n_obs} cells × {adata.n_vars} genes")
        
        return adata
    
    def read_mtx(self, path: str) -> ad.AnnData:
        """读取单个MTX文件"""
        if self.verbose:
            print(f"读取MTX格式: {path}")
            print("  ⚠️  注意：单个MTX文件没有基因名和细胞条形码信息")
        
        adata = sc.read_mtx(path).T  # 转置
        
        if self.verbose:
            print(f"✓ 读取成功: {adata.n_obs} cells × {adata.n_vars} genes")
            print("  建议添加基因名和细胞条形码")
        
        return adata
    
    def read_hdf5(self, path: str, key: str = 'data') -> ad.AnnData:
        """
        读取HDF5格式
        
        参数:
            path: HDF5文件路径
            key: HDF5中的key
        """
        if self.verbose:
            print(f"读取HDF5格式: {path}")
            print(f"  - key: {key}")
        
        adata = sc.read_hdf(path, key=key)
        
        if self.verbose:
            print(f"✓ 读取成功: {adata.n_obs} cells × {adata.n_vars} genes")
        
        return adata
    
    def read_soft_gz(self, path: str) -> ad.AnnData:
        """读取GEO SOFT.GZ格式"""
        if self.verbose:
            print(f"读取SOFT.GZ格式: {path}")
        
        adata = sc.read(path, ext='soft.gz')
        
        if self.verbose:
            print(f"✓ 读取成功: {adata.n_obs} samples × {adata.n_vars} genes")
        
        return adata
    
    def read_auto(
        self,
        input_path: str,
        **kwargs
    ) -> ad.AnnData:
        """
        自动检测格式并读取
        
        参数:
            input_path: 输入文件或目录路径
            **kwargs: 传递给特定读取函数的参数
            
        返回:
            AnnData对象
        """
        # 检测格式
        format_type = self.detect_format(input_path)
        
        if self.verbose:
            print("=" * 70)
            print("通用单细胞数据读取器")
            print("=" * 70)
            print(f"\n检测到格式: {self.supported_formats.get(format_type, format_type)}")
            print(f"输入路径: {input_path}\n")
        
        # 根据格式调用相应的读取函数
        if format_type == '10x_mtx':
            adata = self.read_10x_mtx(input_path, **kwargs)
        elif format_type == '10x_h5':
            adata = self.read_10x_h5(input_path, **kwargs)
        elif format_type == 'h5ad':
            adata = self.read_h5ad(input_path, **kwargs)
        elif format_type == 'loom':
            adata = self.read_loom(input_path)
        elif format_type == 'csv':
            adata = self.read_csv(input_path, **kwargs)
        elif format_type == 'tsv':
            adata = self.read_text(input_path, **kwargs)
        elif format_type == 'excel':
            adata = self.read_excel(input_path, **kwargs)
        elif format_type == 'zarr':
            adata = self.read_zarr(input_path)
        elif format_type == 'mtx':
            adata = self.read_mtx(input_path)
        elif format_type == 'hdf5':
            adata = self.read_hdf5(input_path, **kwargs)
        elif format_type == 'soft_gz':
            adata = self.read_soft_gz(input_path)
        elif format_type == 'rds':
            raise ValueError(
                "RDS格式需要在R中转换为H5AD格式。\n"
                "请在R中运行:\n"
                "  library(SeuratDisk)\n"
                "  Convert('data.rds', dest='h5ad')"
            )
        elif format_type == 'h5seurat':
            raise ValueError(
                "H5Seurat格式需要在R中转换为H5AD格式。\n"
                "请在R中运行:\n"
                "  library(SeuratDisk)\n"
                "  Convert('data.h5seurat', dest='h5ad')"
            )
        else:
            raise ValueError(
                f"未知或不支持的格式: {format_type}\n"
                f"支持的格式: {list(self.supported_formats.keys())}"
            )
        
        return adata
    
    def validate_adata(self, adata: ad.AnnData) -> dict:
        """
        验证和检查AnnData对象
        
        返回:
            包含统计信息的字典
        """
        stats = {
            'n_obs': adata.n_obs,
            'n_vars': adata.n_vars,
            'is_sparse': sp.issparse(adata.X),
            'has_raw': adata.raw is not None,
            'obs_keys': list(adata.obs.columns),
            'var_keys': list(adata.var.columns),
            'obsm_keys': list(adata.obsm.keys()) if adata.obsm else [],
            'uns_keys': list(adata.uns.keys()) if adata.uns else [],
            'layers_keys': list(adata.layers.keys()) if adata.layers else []
        }
        
        if self.verbose:
            print("\n" + "=" * 70)
            print("数据验证和统计")
            print("=" * 70)
            print(f"\n维度信息:")
            print(f"  - 细胞数: {stats['n_obs']:,}")
            print(f"  - 基因数: {stats['n_vars']:,}")
            print(f"  - 是否稀疏矩阵: {stats['is_sparse']}")
            print(f"  - 是否有raw数据: {stats['has_raw']}")
            
            if stats['obs_keys']:
                print(f"\n细胞元数据列 (obs): {len(stats['obs_keys'])} 列")
                print(f"  {', '.join(stats['obs_keys'][:10])}")
                if len(stats['obs_keys']) > 10:
                    print(f"  ... 还有 {len(stats['obs_keys']) - 10} 列")
            
            if stats['var_keys']:
                print(f"\n基因元数据列 (var): {len(stats['var_keys'])} 列")
                print(f"  {', '.join(stats['var_keys'][:10])}")
                if len(stats['var_keys']) > 10:
                    print(f"  ... 还有 {len(stats['var_keys']) - 10} 列")
            
            if stats['obsm_keys']:
                print(f"\n多维数组 (obsm): {stats['obsm_keys']}")
            
            if stats['layers_keys']:
                print(f"\n数据层 (layers): {stats['layers_keys']}")
        
        return stats
    
    def standardize_adata(
        self,
        adata: ad.AnnData,
        ensure_sparse: bool = True,
        ensure_unique_names: bool = True
    ) -> ad.AnnData:
        """
        标准化AnnData对象
        
        参数:
            adata: 输入的AnnData对象
            ensure_sparse: 确保X是稀疏矩阵
            ensure_unique_names: 确保基因名唯一
        """
        if self.verbose:
            print("\n" + "=" * 70)
            print("标准化AnnData对象")
            print("=" * 70)
        
        # 确保使用稀疏矩阵
        if ensure_sparse and not sp.issparse(adata.X):
            if self.verbose:
                print("  - 转换为稀疏矩阵")
            adata.X = sp.csr_matrix(adata.X)
        
        # 确保基因名唯一
        if ensure_unique_names:
            n_vars_before = adata.n_vars
            adata.var_names_make_unique()
            n_duplicates = n_vars_before - len(set(adata.var_names))
            if n_duplicates > 0 and self.verbose:
                print(f"  - 处理了 {n_duplicates} 个重复基因名")
        
        # 确保obs_names和var_names是字符串
        if self.verbose:
            print("  - 确保索引为字符串类型")
        adata.obs_names = adata.obs_names.astype(str)
        adata.var_names = adata.var_names.astype(str)
        
        if self.verbose:
            print("✓ 标准化完成")
        
        return adata
    
    def save_h5ad(
        self,
        adata: ad.AnnData,
        output_path: str,
        compression: str = 'gzip'
    ) -> None:
        """
        保存为H5AD格式
        
        参数:
            adata: AnnData对象
            output_path: 输出文件路径
            compression: 压缩方式 ('gzip', 'lzf', None)
        """
        if self.verbose:
            print("\n" + "=" * 70)
            print("保存H5AD文件")
            print("=" * 70)
            print(f"输出路径: {output_path}")
            print(f"压缩方式: {compression}")
        
        adata.write_h5ad(output_path, compression=compression)
        
        # 检查文件大小
        file_size = os.path.getsize(output_path)
        file_size_mb = file_size / (1024 * 1024)
        
        if self.verbose:
            print(f"✓ 保存成功")
            print(f"  文件大小: {file_size_mb:.2f} MB")
    
    def process_single_sample(
        self,
        input_path: str,
        output_path: Optional[str] = None,
        sample_id: Optional[str] = None,
        **read_kwargs
    ) -> ad.AnnData:
        """
        处理单个样本：自动读取、验证、标准化、保存
        
        参数:
            input_path: 输入文件或目录路径
            output_path: 输出H5AD文件路径（可选）
            sample_id: 样本ID（可选）
            **read_kwargs: 传递给读取函数的参数
            
        返回:
            标准化后的AnnData对象
        """
        # 读取数据
        adata = self.read_auto(input_path, **read_kwargs)
        
        # 添加样本ID
        if sample_id:
            adata.obs['sample_id'] = sample_id
        
        # 验证数据
        stats = self.validate_adata(adata)
        
        # 标准化
        adata = self.standardize_adata(adata)
        
        # 保存（如果指定了输出路径）
        if output_path:
            self.save_h5ad(adata, output_path)
        
        if self.verbose:
            print("\n" + "=" * 70)
            print("处理完成")
            print("=" * 70)
        
        return adata


def main():
    """主函数：演示各种格式的读取"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='通用单细胞RNA-seq数据读取工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 自动检测格式
  python universal_scrna_reader.py input_data/ -o output.h5ad
  
  # 10X MTX格式
  python universal_scrna_reader.py filtered_feature_bc_matrix/ -o output.h5ad
  
  # 10X H5格式
  python universal_scrna_reader.py filtered_feature_bc_matrix.h5 -o output.h5ad
  
  # STARsolo输出（未压缩）
  python universal_scrna_reader.py Solo.out/Gene/filtered/ -o output.h5ad --no-compressed
  
  # CSV格式
  python universal_scrna_reader.py expression.csv -o output.h5ad --transpose
  
  # 指定样本ID
  python universal_scrna_reader.py data/ -o output.h5ad --sample-id sample1
        """
    )
    
    parser.add_argument('input', help='输入文件或目录路径')
    parser.add_argument('-o', '--output', help='输出H5AD文件路径', default=None)
    parser.add_argument('--sample-id', help='样本ID', default=None)
    parser.add_argument('--format', help='强制指定格式（可选）', default=None)
    parser.add_argument('--var-names', help='10X格式：使用gene_symbols或gene_ids',
                       choices=['gene_symbols', 'gene_ids'], default='gene_symbols')
    parser.add_argument('--genome', help='10X H5格式：基因组名称', default=None)
    parser.add_argument('--no-gex-only', help='10X格式：保留所有特征类型',
                       action='store_true')
    parser.add_argument('--no-compressed', help='10X MTX格式：文件未压缩（如STARsolo）',
                       action='store_true')
    parser.add_argument('--delimiter', help='文本格式：分隔符', default='\t')
    parser.add_argument('--transpose', help='文本格式：转置矩阵', action='store_true')
    parser.add_argument('--sheet', help='Excel格式：sheet名称', default='Sheet1')
    parser.add_argument('--backed', help='H5AD格式：使用backed模式',
                       choices=['r', 'r+'], default=None)
    parser.add_argument('--no-cache', help='10X格式：不使用缓存', action='store_true')
    parser.add_argument('--quiet', help='静默模式', action='store_true')
    
    args = parser.parse_args()
    
    # 创建读取器
    reader = UniversalScRNAReader(verbose=not args.quiet)
    
    # 准备读取参数
    read_kwargs = {}
    
    # 10X相关参数
    if args.var_names:
        read_kwargs['var_names'] = args.var_names
    if args.genome:
        read_kwargs['genome'] = args.genome
    if args.no_gex_only:
        read_kwargs['gex_only'] = False
    if args.no_compressed:
        read_kwargs['compressed'] = False
    if not args.no_cache:
        read_kwargs['cache'] = True
    
    # 文本格式参数
    if args.delimiter:
        read_kwargs['delimiter'] = args.delimiter
    if args.transpose:
        read_kwargs['transpose'] = args.transpose
    
    # Excel参数
    if args.sheet:
        read_kwargs['sheet'] = args.sheet
    
    # H5AD参数
    if args.backed:
        read_kwargs['backed'] = args.backed
    
    try:
        # 读取并处理数据
        adata = reader.process_single_sample(
            args.input,
            output_path=args.output,
            sample_id=args.sample_id,
            **read_kwargs
        )
        
        # 打印最终信息
        if not args.quiet:
            print(f"\n最终AnnData对象:")
            print(adata)
        
        return 0
    
    except Exception as e:
        print(f"\n❌ 错误: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())

