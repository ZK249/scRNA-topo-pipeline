# scRNA-topo-pipeline

单细胞转录组数据分析 Pipeline：整合拓扑数据分析（TDA）与标准生物信息学流程，实现细胞亚群识别、拓扑验证、自动注释与发育轨迹推断。

---

## 项目结构

```
scRNA-topo-pipeline/
├── data/                    # 数据目录（gitignore）
│   ├── raw/                 # 10x PBMC 3k 原始矩阵
│   └── processed/           # 清洗后的 AnnData (.h5ad)
├── src/                     # 模块化源码
│   ├── config.py            # 全局参数
│   ├── io.py                # 数据读写
│   ├── qc.py                # 质控
│   ├── preprocess.py        # 标准化、降维
│   ├── clustering.py        # Leiden 聚类
│   ├── topology.py          # 持久同调（Ripser）
│   ├── annotation.py        # 细胞类型自动注释
│   └── trajectory.py        # PAGA + Diffusion Pseudotime
├── notebooks/               # 交互式分析（按阶段）
│   ├── 01_qc_and_preprocess.ipynb
│   ├── 02_clustering_and_topology.ipynb
│   ├── 03_cell_annotation.ipynb
│   └── 04_trajectory_inference.ipynb
├── scripts/
│   └── run_pipeline.py      # 一键运行完整流程
├── markers/
│   └── pbmc_markers.json    # PBMC 细胞类型标记基因库
├── results/                 # 输出结果（gitignore）
│   ├── figures/
│   └── tables/
└── requirements.txt
```

---

## 快速开始

### 1. 环境安装

```bash
# 创建 conda 环境（推荐）
conda create -n scrna python=3.10
conda activate scrna

# 安装依赖
pip install -r requirements.txt
```

### 2. 数据下载

```bash
mkdir -p data/raw
cd data/raw

# 下载 10x PBMC 3k 过滤后的表达矩阵
wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz

cd ../..
```

### 3. 一键运行

```bash
python scripts/run_pipeline.py
```

或分阶段运行 Jupyter Notebook：

```bash
jupyter notebook notebooks/
```

---

## 分析流程

### Phase 1: 质控与预处理
- 读取 10x 标准输出（`matrix.mtx` + `genes.tsv` + `barcodes.tsv`）
- 质控指标：基因数、UMI 数、线粒体比例
- 过滤：200 < n_genes < 2500, pct_mt < 5%
- 标准化：Total UMI 归一化 → log1p → 高变基因选择（top 2000）
- PCA 降维（50 主成分）

### Phase 2: 聚类与拓扑验证
- KNN 图构建（n_neighbors=15, n_pcs=30）
- UMAP 非线性降维
- **Leiden 图聚类**（resolution=0.5）
- **核心创新**：对 PCA 嵌入（前 10 维）计算 **Rips 复形持久同调**
  - 验证聚类结构的拓扑稳定性
  - Platelets 的 H0 平均持久性高达 **11.58**，是其他簇的 3-6 倍

### Phase 3: 细胞类型注释
- Wilcoxon rank-sum 检验计算差异表达基因（DEG）
- 基于 `CellMarker` 知识库构建 PBMC marker 基因映射
- 自动评分 + 类型分配
- 结果：6 个簇 → T cells / NK cells / B cells / Monocytes / Dendritic cells / Platelets

### Phase 4: 发育轨迹推断
- **PAGA**：构建细胞类型间的过渡图
- **Diffusion Map**：捕捉非线性分化流形
- **Diffusion Pseudotime**：以 T cells 为根推断分化方向
- **拓扑 × 轨迹联合分析**：验证拓扑独立性与发育距离的正相关性

---

## 核心结果

| 细胞类型 | 细胞数 | 拓扑持久性 (H0) | 伪时间 | 生物学解读 |
|----------|--------|----------------|--------|-----------|
| T cells | 1181 | 1.93 | 0.020 (根) | 高度异质性，最松散 |
| NK cells | 427 | 3.12 | 0.046 | 与 T 细胞共享淋巴系起源 |
| B cells | 342 | 2.39 | 0.272 | 独立淋巴系 |
| Monocytes | 639 | 2.70 | **0.890** | 髓系，伪时间最远 |
| Dendritic cells | 36 | 3.99 | 0.719 | 髓系末端，稀有亚型 |
| **Platelets** | **13** | **11.58** | **0.768** | **无核，拓扑完全隔离** |

---

## 关键产出

| 文件 | 说明 |
|------|------|
| `results/figures/01_qc_before_filter.png` | 质控分布图 |
| `results/figures/02_umap_leiden.png` | UMAP 聚类图 |
| `results/figures/02_persistence_full.png` | 全量数据持久图 |
| `results/figures/02_persistence_by_cluster.png` | 各簇持久图对比 |
| `results/figures/03_umap_cell_type.png` | 细胞类型注释 UMAP |
| `results/figures/04_paga_graph.png` | PAGA 细胞类型连接图 |
| `results/figures/04_umap_dpt.png` | 伪时间 UMAP 映射 |
| `results/figures/04_topo_vs_dpt.png` | **拓扑持久性 vs 发育距离** |
| `results/tables/03_cluster_annotation.csv` | 簇注释表 |
| `results/tables/02_cluster_topology_stats.csv` | 拓扑统计表 |

---

## 技术栈

- **Scanpy** / **AnnData**: 单细胞数据结构与标准分析
- **Leiden** / **UMAP**: 图聚类与降维
- **Ripser**: 快速持久同调计算
- **PAGA** / **Diffusion Map**: 轨迹推断
- **Seaborn** / **Matplotlib**: 可视化

---

## License

MIT
