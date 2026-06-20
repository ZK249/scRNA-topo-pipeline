"""
数据 I/O 封装。屏蔽底层格式细节，统一接口。
"""

import scanpy as sc
import anndata
from pathlib import Path
from .config import RAW_DIR, PROCESSED_DIR


def read_10x_filtered(matrix_dir: Path = None) -> sc.AnnData:
    """
    读取 10x Genomics 过滤后的标准输出。
    
    Args:
        matrix_dir: 包含 matrix.mtx, genes.tsv, barcodes.tsv 的目录
                   默认使用 config 中配置的 RAW_DIR
    
    Returns:
        AnnData 对象，原始 UMI 计数
    """
    if matrix_dir is None:
        matrix_dir = RAW_DIR
    
    if not matrix_dir.exists():
        raise FileNotFoundError(
            f"数据目录不存在: {matrix_dir}\n"
            f"请确认已下载并解压 pbmc3k_filtered_gene_bc_matrices.tar.gz 到 data/raw/"
        )
    
    adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", cache=True)
    adata.var_names_make_unique()  # 处理重复基因名
    return adata


def save_h5ad(adata: sc.AnnData, filename: str, dir_path: Path = None) -> Path:
    """保存 AnnData 为 h5ad 格式，返回保存路径。"""
    if dir_path is None:
        dir_path = PROCESSED_DIR
    dir_path.mkdir(parents=True, exist_ok=True)
    out_path = dir_path / filename
    adata.write_h5ad(out_path)
    return out_path


def load_h5ad(filename: str, dir_path: Path = None) -> sc.AnnData:
    """加载 h5ad 文件。"""
    if dir_path is None:
        dir_path = PROCESSED_DIR
    return sc.read_h5ad(dir_path / filename)