"""
聚类与降维模块。标准单细胞分析流程 + 拓扑增强验证。
"""

import scanpy as sc
from . import config


def run_neighbors(adata, n_pcs=None, n_neighbors=None):
    """构建 KNN 图，基于 PCA 嵌入。"""
    if n_pcs is None:
        n_pcs = config.N_PCS_CLUSTERING
    if n_neighbors is None:
        n_neighbors = config.NEIGHBORS_N
    
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    return adata


def run_umap(adata, min_dist=0.5):
    """UMAP 非线性降维。"""
    sc.tl.umap(adata, min_dist=min_dist)
    return adata


def run_leiden(adata, resolution=None, key_added="leiden"):
    """Leiden 图聚类。"""
    if resolution is None:
        resolution = config.LEIDEN_RESOLUTION
    sc.tl.leiden(adata, resolution=resolution, key_added=key_added)
    return adata