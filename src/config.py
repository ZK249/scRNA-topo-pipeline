"""
全局配置参数。所有路径、阈值、分析参数集中管理。
修改这里即可适配其他数据集。
"""

from pathlib import Path

# ========== 路径配置 ==========
PROJECT_ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw" / "pbmc3k_filtered_gene_bc_matrices" / "hg19"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
TABLES_DIR = RESULTS_DIR / "tables"

# 确保目录存在
for d in [PROCESSED_DIR, FIGURES_DIR, TABLES_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ========== 质控阈值 ==========
# 基于 10x PBMC 3k 经典参数，可根据实际数据分布调整
QC_MIN_GENES = 200          # 每个细胞至少检测到的基因数
QC_MAX_GENES = 2500         # 上限，排除潜在双细胞
QC_MIN_CELLS = 3            # 每个基因至少在多少个细胞中表达
QC_MAX_PCT_MT = 5           # 线粒体基因比例上限（%）

# ========== 分析参数 ==========
N_TOP_GENES = 2000            # 高变基因数
N_PCS = 50                    # PCA 主成分数
LEIDEN_RESOLUTION = 0.5       # 聚类分辨率（越小簇越少）
NEIGHBORS_N = 15              # KNN 邻居数

# ========== 聚类与拓扑参数 ==========
NEIGHBORS_N = 15               # KNN 邻居数
N_PCS_CLUSTERING = 30          # 聚类时使用的 PCA 维度
LEIDEN_RESOLUTION = 0.5        # 聚类分辨率

# 拓扑分析
TOPO_MAXDIM = 1                # 计算 0 维(连通分支) + 1 维(环)
TOPO_PCS = 10                  # 持久同调使用的 PCA 维度（太高 ripser 会慢）