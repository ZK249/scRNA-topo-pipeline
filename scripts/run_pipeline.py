"""
scRNA-topo-pipeline: 一键运行完整分析流程。

Usage:
    python scripts/run_pipeline.py

输出：
    - data/processed/pbmc3k_*.h5ad
    - results/figures/*.png
    - results/tables/*.csv
"""

import sys
from pathlib import Path

# 将项目根目录加入 Python 路径
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

sc.settings.verbosity = 2  # 减少日志输出
sc.settings.set_figure_params(dpi=100, facecolor="white")

from src import config
from src.io import read_10x_filtered, save_h5ad
from src.clustering import run_neighbors, run_umap, run_leiden
from src.topology import compute_persistence, plot_persistence_diagrams, extract_topological_stats
from src.annotation import (
    load_marker_db,
    rank_cluster_genes,
    score_cell_types,
    annotate_clusters,
    add_annotation_to_adata,
)
from src.trajectory import run_paga, run_diffmap, run_dpt


def phase1_qc_preprocess():
    """Phase 1: 质控与预处理"""
    print("\n" + "=" * 50)
    print("Phase 1: QC & Preprocessing")
    print("=" * 50)

    adata = read_10x_filtered()
    print(f"原始数据: {adata.n_obs} cells × {adata.n_vars} genes")

    # 质控指标
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # 过滤
    adata = adata[adata.obs.n_genes_by_counts > config.QC_MIN_GENES, :]
    adata = adata[adata.obs.n_genes_by_counts < config.QC_MAX_GENES, :]
    adata = adata[adata.obs.pct_counts_mt < config.QC_MAX_PCT_MT, :]
    sc.pp.filter_genes(adata, min_cells=config.QC_MIN_CELLS)
    print(f"过滤后: {adata.n_obs} cells × {adata.n_vars} genes")

    # 标准化
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=config.N_TOP_GENES, subset=True)
    adata.raw = adata

    # PCA
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=config.N_PCS)

    save_h5ad(adata, "pbmc3k_qc_normalized.h5ad")
    print("Phase 1 完成，数据已保存")
    return adata


def phase2_clustering_topology(adata):
    """Phase 2: 聚类 + 拓扑验证"""
    print("\n" + "=" * 50)
    print("Phase 2: Clustering & Topological Validation")
    print("=" * 50)

    # 聚类
    run_neighbors(adata, n_pcs=config.N_PCS_CLUSTERING, n_neighbors=config.NEIGHBORS_N)
    run_umap(adata)
    run_leiden(adata, resolution=config.LEIDEN_RESOLUTION)

    n_clusters = adata.obs["leiden"].nunique()
    print(f"Leiden 聚类: {n_clusters} clusters")

    # 全量拓扑
    pca_topo = adata.obsm["X_pca"][:, :config.TOPO_PCS]
    print(f"计算全量持久同调: {pca_topo.shape}")
    result_full = compute_persistence(pca_topo, maxdim=config.TOPO_MAXDIM)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    plot_persistence_diagrams(result_full["dgms"], ax=axes, title="Full PBMC Data")
    plt.tight_layout()
    plt.savefig(config.FIGURES_DIR / "02_persistence_full.png", dpi=150, bbox_inches="tight")
    plt.close()

    # 各簇拓扑
    cluster_ids = adata.obs["leiden"].cat.categories
    cluster_stats = []
    fig, axes = plt.subplots(len(cluster_ids), 2, figsize=(12, 4.5 * len(cluster_ids)))

    for i, cid in enumerate(cluster_ids):
        mask = adata.obs["leiden"] == cid
        coords = adata.obsm["X_pca"][mask, :config.TOPO_PCS]
        res = compute_persistence(coords, maxdim=config.TOPO_MAXDIM)
        stats = extract_topological_stats(res["dgms"])
        stats["cluster"] = cid
        stats["n_cells"] = mask.sum()
        cluster_stats.append(stats)
        plot_persistence_diagrams(res["dgms"], ax=axes[i], title=f"Cluster {cid} (n={mask.sum()})")

    plt.tight_layout()
    plt.savefig(config.FIGURES_DIR / "02_persistence_by_cluster.png", dpi=150, bbox_inches="tight")
    plt.close()

    df_stats = pd.DataFrame(cluster_stats)
    df_stats.to_csv(config.TABLES_DIR / "02_cluster_topology_stats.csv", index=False)

    # UMAP 图
    fig, ax = plt.subplots(figsize=(8, 7))
    sc.pl.umap(adata, color="leiden", ax=ax, show=False, legend_loc="on data",
               palette="tab20", title=f"Leiden Clustering (n={n_clusters})")
    plt.tight_layout()
    plt.savefig(config.FIGURES_DIR / "02_umap_leiden.png", dpi=150, bbox_inches="tight")
    plt.close()

    save_h5ad(adata, "pbmc3k_clustered.h5ad")
    print("Phase 2 完成")
    return adata, df_stats


def phase3_annotation(adata, df_stats):
    """Phase 3: 细胞类型注释"""
    print("\n" + "=" * 50)
    print("Phase 3: Cell Type Annotation")
    print("=" * 50)

    rank_cluster_genes(adata, group_key="leiden", method="wilcoxon")

    marker_db = load_marker_db()
    score_df = score_cell_types(adata, marker_db)
    annotations = annotate_clusters(score_df, top_n=1)

    print("自动注释结果:")
    for k, v in annotations.items():
        print(f"  Cluster {k}: {v}")

    add_annotation_to_adata(adata, annotations, score_df)

    # 保存注释表
    annot_df = pd.DataFrame({
        "leiden": adata.obs["leiden"].cat.categories,
        "cell_type": [annotations[c] for c in adata.obs["leiden"].cat.categories],
        "n_cells": adata.obs["leiden"].value_counts().sort_index().values,
    })
    annot_df.to_csv(config.TABLES_DIR / "03_cluster_annotation.csv", index=False)

    # UMAP
    fig, ax = plt.subplots(figsize=(8, 7))
    sc.pl.umap(adata, color="cell_type", ax=ax, show=False, legend_loc="on data", title="Cell Type Annotation")
    plt.tight_layout()
    plt.savefig(config.FIGURES_DIR / "03_umap_cell_type.png", dpi=150, bbox_inches="tight")
    plt.close()

    save_h5ad(adata, "pbmc3k_annotated.h5ad")
    print("Phase 3 完成")
    return adata


def phase4_trajectory(adata, df_stats):
    """Phase 4: 轨迹推断"""
    print("\n" + "=" * 50)
    print("Phase 4: Trajectory Inference")
    print("=" * 50)

    # PAGA
    run_paga(adata, groups="cell_type")
    fig, ax = plt.subplots(figsize=(8, 8))
    sc.pl.paga(adata, color="cell_type", node_size_scale=3, edge_width_scale=1, ax=ax, show=False,
               title="PAGA: Cell Type Connectivity")
    plt.tight_layout()
    plt.savefig(config.FIGURES_DIR / "04_paga_graph.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Diffusion Map + DPT
    run_diffmap(adata, n_comps=15)
    run_dpt(adata, root_group="T cells", root_key="cell_type")

    print("各细胞类型平均伪时间:")
    dpt_by_type = adata.obs.groupby("cell_type")["dpt_pseudotime"].agg(["mean", "std", "count"])
    print(dpt_by_type.round(3))

    # UMAP + DPT
    fig, ax = plt.subplots(figsize=(8, 7))
    sc.pl.umap(adata, color="dpt_pseudotime", ax=ax, show=False, title="Diffusion Pseudotime on UMAP", color_map="viridis")
    plt.tight_layout()
    plt.savefig(config.FIGURES_DIR / "04_umap_dpt.png", dpi=150, bbox_inches="tight")
    plt.close()

    # PAGA + DPT
    fig, ax = plt.subplots(figsize=(8, 8))
    sc.pl.paga(adata, color="dpt_pseudotime", node_size_scale=3, edge_width_scale=1, ax=ax, show=False,
               title="PAGA colored by Pseudotime")
    plt.tight_layout()
    plt.savefig(config.FIGURES_DIR / "04_paga_dpt.png", dpi=150, bbox_inches="tight")
    plt.close()

    # Violin
    fig, ax = plt.subplots(figsize=(10, 5))
    dpt_df = adata.obs[["cell_type", "dpt_pseudotime"]].copy()
    sns.violinplot(data=dpt_df, x="cell_type", y="dpt_pseudotime", palette="tab10", inner="box", ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_title("Pseudotime Distribution by Cell Type")
    plt.tight_layout()
    plt.savefig(config.FIGURES_DIR / "04_dpt_violin.png", dpi=150, bbox_inches="tight")
    plt.close()

    # 拓扑 × 轨迹联合
    df_stats["cell_type"] = df_stats["cluster"].astype(str).map(
        dict(zip(adata.obs["leiden"].cat.categories.astype(str),
                 adata.obs.groupby("leiden")["cell_type"].first()))
    )
    merged = df_stats.merge(dpt_by_type.reset_index(), on="cell_type", how="left")

    fig, ax = plt.subplots(figsize=(8, 6))
    for ct in merged["cell_type"]:
        row = merged[merged["cell_type"] == ct].iloc[0]
        ax.scatter(row["mean"], row["mean_persistence_h0"], s=row["count"] * 2, alpha=0.7, label=ct)
        ax.annotate(ct, (row["mean"], row["mean_persistence_h0"]), fontsize=9, ha="center", va="bottom")
    ax.set_xlabel("Mean Diffusion Pseudotime")
    ax.set_ylabel("H0 Mean Persistence")
    ax.set_title("Topological Independence vs Developmental Distance")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(config.FIGURES_DIR / "04_topo_vs_dpt.png", dpi=150, bbox_inches="tight")
    plt.close()

    save_h5ad(adata, "pbmc3k_final.h5ad")
    print("Phase 4 完成，最终数据已保存")
    return adata


def main():
    print("=" * 50)
    print("scRNA-topo-pipeline: 一键运行")
    print("=" * 50)

    # Phase 1
    adata = phase1_qc_preprocess()

    # Phase 2
    adata, df_stats = phase2_clustering_topology(adata)

    # Phase 3
    adata = phase3_annotation(adata, df_stats)

    # Phase 4
    adata = phase4_trajectory(adata, df_stats)

    print("\n" + "=" * 50)
    print("全部完成！")
    print(f"结果目录: {config.RESULTS_DIR}")
    print("=" * 50)


if __name__ == "__main__":
    main()