"""
轨迹推断模块。使用 PAGA 构建细胞类型图 + Diffusion Pseudotime 推断分化方向。
"""

import scanpy as sc
import numpy as np
from . import config


def run_paga(adata, groups="cell_type"):
    """
    PAGA：基于分区图的抽象，连接不同细胞类型，揭示过渡关系。
    """
    sc.tl.paga(adata, groups=groups)
    return adata


def run_diffmap(adata, n_comps=15):
    """
    Diffusion Map：捕捉非线性流形结构，比 PCA 更适合轨迹数据。
    """
    sc.tl.diffmap(adata, n_comps=n_comps)
    return adata


def run_dpt(adata, root_group="T cells", root_key="cell_type"):
    """
    Diffusion Pseudotime：从根节点计算每个细胞的"伪时间"。
    
    Args:
        root_group: 作为分化起点的细胞类型
        root_key: 存储细胞类型标签的列名
    """
    # 找到 root_group 中 Diffusion Map 坐标的质心作为根
    root_mask = adata.obs[root_key] == root_group
    if not root_mask.any():
        raise ValueError(f"Root group '{root_group}' not found in adata.obs['{root_key}']")
    
    # 计算该细胞群在 Diffusion Map 中的平均坐标
    root_coords = adata.obsm["X_diffmap"][root_mask]
    root_mean = root_coords.mean(axis=0)
    
    # 找到最接近质心的细胞作为根
    distances = np.linalg.norm(root_coords - root_mean, axis=1)
    root_idx = np.where(root_mask)[0][np.argmin(distances)]
    
    adata.uns["iroot"] = root_idx
    sc.tl.dpt(adata)
    return adata