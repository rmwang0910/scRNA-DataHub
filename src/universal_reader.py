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
from typing import Optional, Literal, Union
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
            # 尝试判断是否是10X H5格式
            try:
                import h5py
                with h5py.File(path, 'r') as f:
                    keys = list(f.keys())
                    # 10X H5 v3: 有'matrix'键
                    if 'matrix' in keys:
                        return '10x_h5'
                    # 10X H5 v2: 第一层是genome名称，下面有'matrix', 'barcodes', 'gene_names'等
                    elif keys and len(keys) > 0:
                        first_key = keys[0]
                        if isinstance(f[first_key], h5py.Group):
                            subkeys = list(f[first_key].keys())
                            # 检查是否有10X特征性的键
                            if any(k in subkeys for k in ['matrix', 'barcodes', 'gene_names', 'genes']):
                                return '10x_h5'
                    # 不支持通用HDF5格式
                    return 'unknown'
            except:
                return 'unknown'
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
        
        # 检查features.tsv文件列数（处理DNB C4格式）
        from pathlib import Path
        import gzip
        
        path_obj = Path(path)
        suffix = '.gz' if compressed else ''
        
        # 检查features.tsv或genes.tsv
        features_file = path_obj / f'features.tsv{suffix}'
        genes_file = path_obj / f'genes.tsv{suffix}'
        
        feature_file_to_check = features_file if features_file.exists() else genes_file
        
        if feature_file_to_check.exists():
            # 读取第一行检查列数
            try:
                if str(feature_file_to_check).endswith('.gz'):
                    with gzip.open(feature_file_to_check, 'rt') as f:
                        first_line = f.readline().strip()
                else:
                    with open(feature_file_to_check, 'r') as f:
                        first_line = f.readline().strip()
                
                n_cols = len(first_line.split('\t'))
                
                if self.verbose:
                    print(f"  - 特征文件列数: {n_cols}")
                
                # 如果只有1列（DNB C4格式），需要手动处理
                if n_cols == 1:
                    if self.verbose:
                        print("  ℹ️  检测到DNB C4格式（features.tsv只有1列）")
                        print("  - 使用自定义读取方法...")
                    
                    # 手动读取DNB C4格式
                    adata = self._read_dnb_c4_mtx(path, var_names, make_unique, compressed)
                    
                    if self.verbose:
                        print(f"✓ 读取成功: {adata.n_obs} cells × {adata.n_vars} genes")
                    
                    return adata
                    
            except Exception as e:
                if self.verbose:
                    print(f"  ⚠️  特征文件检查失败: {e}")
        
        # 标准10X格式处理
        # 检查scanpy版本，判断是否支持compressed参数
        import inspect
        sig = inspect.signature(sc.read_10x_mtx)
        supports_compressed = 'compressed' in sig.parameters
        
        # 根据版本决定是否传递compressed参数
        if supports_compressed:
            adata = sc.read_10x_mtx(
                path,
                var_names=var_names,
                make_unique=make_unique,
                cache=cache,
                gex_only=gex_only,
                compressed=compressed
            )
        else:
            # 旧版本scanpy不支持compressed参数
            if self.verbose and not compressed:
                print("  ⚠️  当前scanpy版本不支持compressed参数，将尝试自动检测")
            
            adata = sc.read_10x_mtx(
                path,
                var_names=var_names,
                make_unique=make_unique,
                cache=cache,
                gex_only=gex_only
            )
        
        if self.verbose:
            print(f"✓ 读取成功: {adata.n_obs} cells × {adata.n_vars} genes")
        
        return adata
    
    def _read_dnb_c4_mtx(
        self,
        path: str,
        var_names: str = 'gene_symbols',
        make_unique: bool = True,
        compressed: bool = True
    ) -> ad.AnnData:
        """
        读取DNB C4格式的MTX数据
        DNB C4的features.tsv只有1列（基因名），需要特殊处理
        
        参数:
            path: 数据目录
            var_names: 变量名类型（对DNB C4只有gene_symbols可用）
            make_unique: 是否使基因名唯一
            compressed: 是否压缩
        """
        from pathlib import Path
        from scipy.io import mmread
        import gzip
        
        path_obj = Path(path)
        suffix = '.gz' if compressed else ''
        
        # 读取barcodes
        barcodes_file = path_obj / f'barcodes.tsv{suffix}'
        if compressed:
            with gzip.open(barcodes_file, 'rt') as f:
                barcodes = [line.strip() for line in f]
        else:
            with open(barcodes_file, 'r') as f:
                barcodes = [line.strip() for line in f]
        
        # 读取features（只有1列：基因名）
        features_file = path_obj / f'features.tsv{suffix}'
        if not features_file.exists():
            features_file = path_obj / f'genes.tsv{suffix}'
        
        if compressed:
            with gzip.open(features_file, 'rt') as f:
                gene_names = [line.strip() for line in f]
        else:
            with open(features_file, 'r') as f:
                gene_names = [line.strip() for line in f]
        
        # 读取matrix
        matrix_file = path_obj / f'matrix.mtx{suffix}'
        if compressed:
            # 需要先解压
            import tempfile
            import shutil
            with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.mtx') as tmp:
                with gzip.open(matrix_file, 'rb') as f_in:
                    shutil.copyfileobj(f_in, tmp)
                tmp_path = tmp.name
            
            matrix = mmread(tmp_path).T.tocsr()  # 转置为 cells × genes
            Path(tmp_path).unlink()  # 删除临时文件
        else:
            matrix = mmread(matrix_file).T.tocsr()
        
        # 创建AnnData对象
        adata = ad.AnnData(X=matrix)
        adata.obs_names = barcodes
        adata.var_names = gene_names
        
        # 添加gene_ids列（与gene_names相同）
        adata.var['gene_ids'] = gene_names
        adata.var['feature_types'] = ['Gene Expression'] * len(gene_names)
        
        # 处理重复基因名
        if make_unique:
            import anndata.utils
            adata.var_names = anndata.utils.make_index_unique(pd.Index(adata.var_names))
        
        if self.verbose:
            print(f"  - DNB C4格式：features.tsv只有1列（基因名）")
            print(f"  - 自动添加gene_ids和feature_types列")
        
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
            delimiter: 分隔符（None表示自动检测）
            first_column_names: 第一列是否是行名
            transpose: 是否转置
        """
        # 特殊文件处理：krumsiek11是scanpy内置数据集
        if 'krumsiek11' in str(path).lower():
            if self.verbose:
                print(f"读取scanpy内置数据集: krumsiek11")
            try:
                adata = sc.datasets.krumsiek11()
                if self.verbose:
                    print(f"✓ 读取成功: {adata.n_obs} cells × {adata.n_vars} genes")
                return adata
            except Exception as e:
                if self.verbose:
                    print(f"  ⚠️  使用内置数据集失败: {e}")
                    print(f"  - 尝试使用pandas读取多空格分隔文件")
        
        # 智能分隔符检测
        if delimiter == '\t' or delimiter is None:
            # 尝试检测实际的分隔符
            try:
                import gzip
                opener = gzip.open if path.endswith('.gz') else open
                with opener(path, 'rt') as f:
                    # 跳过注释行
                    line = f.readline()
                    while line.startswith('#'):
                        line = f.readline()
                    
                    # 检测分隔符
                    if '\t' in line:
                        detected_delimiter = '\t'
                    elif ',' in line:
                        detected_delimiter = ','
                    elif '  ' in line:
                        # 多个空格，需要用pandas处理
                        detected_delimiter = 'whitespace'
                    else:
                        detected_delimiter = delimiter
                    
                    if detected_delimiter != delimiter and self.verbose:
                        print(f"读取文本格式: {path}")
                        print(f"  - 检测到分隔符: '{detected_delimiter}' (而不是 '{delimiter}')")
                    delimiter = detected_delimiter
            except Exception as e:
                if self.verbose:
                    print(f"  ⚠️  分隔符检测失败，使用默认值: {e}")
        
        if self.verbose:
            print(f"读取文本格式: {path}")
            print(f"  - delimiter: '{delimiter}'")
        
        # 对于多空格分隔的文件，使用pandas读取
        if delimiter == 'whitespace' or '  ' in delimiter:
            if self.verbose:
                print(f"  - 检测到多空格分隔，使用pandas读取")
            try:
                df = pd.read_csv(path, sep=r'\s+', comment='#', index_col=0)
                adata = ad.AnnData(X=df.T.values, obs=pd.DataFrame(index=df.columns), 
                                  var=pd.DataFrame(index=df.index))
                if self.verbose:
                    print(f"✓ 读取成功: {adata.n_obs} cells × {adata.n_vars} genes")
                return adata
            except Exception as e:
                if self.verbose:
                    print(f"  ⚠️  pandas读取失败: {e}")
                # 继续尝试scanpy的方法
                delimiter = '\t'
        
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
        sheet: Union[str, int] = 0,
        transpose: bool = False
    ) -> ad.AnnData:
        """
        读取Excel格式
        
        参数:
            path: Excel文件路径
            sheet: Sheet名称或索引（默认0表示第一个sheet）
            transpose: 是否转置
        """
        if self.verbose:
            print(f"读取Excel格式: {path}")
            print(f"  - sheet: {sheet}")
        
        # 如果sheet是字符串数字，转换为整数
        if isinstance(sheet, str) and sheet.isdigit():
            sheet = int(sheet)
        
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
        
        # scanpy没有read_zarr，使用anndata的read_zarr
        try:
            # 首先尝试使用scanpy（如果有）
            adata = sc.read_zarr(path)
        except AttributeError:
            # 如果scanpy没有，使用anndata
            if self.verbose:
                print("  - 使用 anndata.read_zarr()")
            adata = ad.read_zarr(path)
        
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
        
        # 根据格式调用相应的读取函数，过滤不适用的参数
        if format_type == '10x_mtx':
            # 10X MTX格式支持的参数
            valid_params = ['var_names', 'make_unique', 'cache', 'cache_compression', 
                          'gex_only', 'prefix', 'compressed']
            filtered_kwargs = {k: v for k, v in kwargs.items() if k in valid_params}
            adata = self.read_10x_mtx(input_path, **filtered_kwargs)
        elif format_type == '10x_h5':
            # 10X H5格式支持的参数
            valid_params = ['genome', 'gex_only', 'backup_url']
            filtered_kwargs = {k: v for k, v in kwargs.items() if k in valid_params}
            adata = self.read_10x_h5(input_path, **filtered_kwargs)
        elif format_type == 'h5ad':
            # H5AD格式支持的参数
            valid_params = ['backed', 'chunk_size']
            filtered_kwargs = {k: v for k, v in kwargs.items() if k in valid_params}
            adata = self.read_h5ad(input_path, **filtered_kwargs)
        elif format_type == 'loom':
            adata = self.read_loom(input_path)
        elif format_type == 'csv':
            # CSV格式支持的参数
            valid_params = ['first_column_names', 'transpose']
            filtered_kwargs = {k: v for k, v in kwargs.items() if k in valid_params}
            adata = self.read_csv(input_path, **filtered_kwargs)
        elif format_type == 'tsv':
            # 文本格式支持的参数
            valid_params = ['delimiter', 'first_column_names', 'transpose']
            filtered_kwargs = {k: v for k, v in kwargs.items() if k in valid_params}
            adata = self.read_text(input_path, **filtered_kwargs)
        elif format_type == 'excel':
            # Excel格式支持的参数
            valid_params = ['sheet', 'transpose']
            filtered_kwargs = {k: v for k, v in kwargs.items() if k in valid_params}
            adata = self.read_excel(input_path, **filtered_kwargs)
        elif format_type == 'zarr':
            adata = self.read_zarr(input_path)
        elif format_type == 'mtx':
            adata = self.read_mtx(input_path)
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
            
            # 检查数据类型，如果是object类型，先转换为float
            if adata.X.dtype == object:
                if self.verbose:
                    print("  - 检测到object类型，转换为float64")
                try:
                    adata.X = adata.X.astype(np.float64)
                except (ValueError, TypeError) as e:
                    if self.verbose:
                        print(f"  ⚠️  转换失败: {e}")
                        print("  - 尝试逐元素转换")
                    # 尝试使用pandas转换
                    import pandas as pd
                    adata.X = pd.DataFrame(adata.X).apply(pd.to_numeric, errors='coerce').fillna(0).values.astype(np.float64)
            
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
        
        # 清理uns中的非字符串键（H5AD格式要求）
        if hasattr(adata, 'uns') and adata.uns:
            cleaned_uns = {}
            for key, value in adata.uns.items():
                # 如果值是字典，检查其键是否为整数
                if isinstance(value, dict):
                    cleaned_dict = {}
                    for k, v in value.items():
                        # 将整数键转换为字符串
                        str_key = str(k) if not isinstance(k, str) else k
                        cleaned_dict[str_key] = v
                    cleaned_uns[key] = cleaned_dict
                else:
                    cleaned_uns[key] = value
            adata.uns = cleaned_uns
            if self.verbose and len(cleaned_uns) != len(adata.uns):
                print(f"  - 清理了uns中的非字符串键")
        
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

