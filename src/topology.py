"""
拓扑数据分析模块。对单细胞 PCA 嵌入计算持久同调，
验证聚类结构的拓扑稳定性。
"""

import numpy as np
import matplotlib.pyplot as plt
from ripser import ripser


def compute_persistence(point_cloud, maxdim=1):
    """
    计算点云的持久同调。
    
    Args:
        point_cloud: (n_points, n_dims) 数组
        maxdim: 最高维度（0=连通分支, 1=环）
    
    Returns:
        dict, ripser 输出，包含 'dgms'（持久图列表）
    """
    result = ripser(point_cloud, maxdim=maxdim, metric="euclidean")
    return result


def plot_persistence_diagrams(diagrams, ax=None, title="Persistence Diagram"):
    """
    绘制持久图（H0 和 H1）。
    
    Args:
        diagrams: list of arrays, ripser['dgms']
        ax: matplotlib axes (1x2)
    """
    if ax is None:
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # H0: 连通分支
    d0 = diagrams[0]
    # 过滤无穷远点（death=inf）
    d0_finite = d0[d0[:, 1] != np.inf]
    birth0, death0 = d0_finite[:, 0], d0_finite[:, 1]
    persistence0 = death0 - birth0
    
    ax[0].scatter(birth0, death0, s=20, alpha=0.6, c="tab:blue")
    ax[0].plot([0, max(death0.max(), birth0.max())], 
               [0, max(death0.max(), birth0.max())], 
               "k--", alpha=0.3, label="diagonal")
    ax[0].set_xlabel("Birth")
    ax[0].set_ylabel("Death")
    ax[0].set_title(f"H0 (Connected Components)\nN={len(d0_finite)}, "
                    f"Mean Persistence={persistence0.mean():.3f}")
    ax[0].legend()
    
    # H1: 环
    if len(diagrams) > 1 and len(diagrams[1]) > 0:
        d1 = diagrams[1]
        birth1, death1 = d1[:, 0], d1[:, 1]
        persistence1 = death1 - birth1
        ax[1].scatter(birth1, death1, s=20, alpha=0.6, c="tab:red")
        ax[1].plot([0, max(death1.max(), birth1.max())], 
                   [0, max(death1.max(), birth1.max())], 
                   "k--", alpha=0.3, label="diagonal")
        ax[1].set_xlabel("Birth")
        ax[1].set_ylabel("Death")
        ax[1].set_title(f"H1 (Loops)\nN={len(d1)}, "
                        f"Mean Persistence={persistence1.mean():.3f}")
        ax[1].legend()
    else:
        ax[1].text(0.5, 0.5, "No H1 features detected", 
                   ha="center", va="center", transform=ax[1].transAxes)
        ax[1].set_title("H1 (Loops)")
    
    plt.suptitle(title)
    return ax


def extract_topological_stats(diagrams):
    """
    从持久图中提取统计特征。
    
    Returns:
        dict: n_h0, n_h1, mean_persistence_h0, mean_persistence_h1
    """
    stats = {"n_h0": 0, "n_h1": 0, "mean_persistence_h0": 0.0, "mean_persistence_h1": 0.0}
    
    # H0
    d0 = diagrams[0]
    d0_finite = d0[d0[:, 1] != np.inf]
    stats["n_h0"] = len(d0_finite)
    if len(d0_finite) > 0:
        stats["mean_persistence_h0"] = (d0_finite[:, 1] - d0_finite[:, 0]).mean()
    
    # H1
    if len(diagrams) > 1 and len(diagrams[1]) > 0:
        d1 = diagrams[1]
        stats["n_h1"] = len(d1)
        stats["mean_persistence_h1"] = (d1[:, 1] - d1[:, 0]).mean()
    
    return stats