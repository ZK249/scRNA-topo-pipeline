"""
细胞类型自动注释模块。
基于差异表达基因（DEG）和已知 marker 库进行注释。
"""

import json
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path


def load_marker_db(marker_path=None):
    """加载 marker 基因库。"""
    if marker_path is None:
        marker_path = Path(__file__).resolve().parent.parent / "markers" / "pbmc_markers.json"
    with open(marker_path, "r") as f:
        return json.load(f)


def rank_cluster_genes(adata, group_key="leiden", method="wilcoxon"):
    """
    计算每个簇的差异表达基因（DEG）。
    使用 scanpy 内置的 rank_genes_groups。
    """
    # 使用 raw 数据（包含全部基因）进行 DEG 分析
    sc.tl.rank_genes_groups(
        adata, 
        groupby=group_key, 
        method=method, 
        use_raw=True,
        key_added="rank_genes"
    )
    return adata


def get_top_markers(adata, cluster_id, n_genes=10, key="rank_genes"):
    """获取指定簇的 top DEG。"""
    return sc.get.rank_genes_groups_df(adata, group=cluster_id, key=key).head(n_genes)


def score_cell_types(adata, marker_db):
    """
    基于 marker 基因库对每个簇进行评分。
    计算每个簇中 marker 基因的平均表达量。
    
    Returns:
        DataFrame: 簇 × 细胞类型 的得分矩阵
    """
    # 使用 raw 数据（log1p 标准化后的全部基因）
    expr = pd.DataFrame(
        adata.raw.X.toarray() if hasattr(adata.raw.X, "toarray") else adata.raw.X,
        index=adata.obs_names,
        columns=adata.raw.var_names
    )
    
    scores = {}
    clusters = adata.obs["leiden"].cat.categories
    
    for cell_type, info in marker_db.items():
        markers = info["markers"]
        # 只保留存在于数据中的基因
        valid_markers = [m for m in markers if m in expr.columns]
        if not valid_markers:
            continue
        
        # 计算每个细胞中该细胞类型 marker 的平均表达
        cell_scores = expr[valid_markers].mean(axis=1)
        
        # 按簇聚合
        cluster_scores = adata.obs[["leiden"]].copy()
        cluster_scores["score"] = cell_scores
        mean_by_cluster = cluster_scores.groupby("leiden")["score"].mean()
        scores[cell_type] = mean_by_cluster
    
    score_df = pd.DataFrame(scores).fillna(0)
    return score_df


def annotate_clusters(score_df, top_n=1):
    """
    根据得分矩阵为每个簇分配最可能的细胞类型。
    
    Returns:
        dict: {cluster_id: cell_type}
    """
    annotations = {}
    for cluster in score_df.index:
        top_types = score_df.loc[cluster].sort_values(ascending=False).head(top_n)
        annotations[cluster] = top_types.index[0]
    return annotations


def add_annotation_to_adata(adata, annotations, score_df=None):
    """
    将注释结果写入 adata.obs。
    """
    adata.obs["cell_type"] = adata.obs["leiden"].map(annotations)
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
    
    if score_df is not None:
        # 保存得分到 uns
        adata.uns["cell_type_scores"] = score_df
    return adata