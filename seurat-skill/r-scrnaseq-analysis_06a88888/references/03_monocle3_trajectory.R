# =============================================================================
# 03_monocle3_trajectory.R
# Monocle3 轨迹推断与拟时序（pseudotime）分析模板
# 官方文档: https://cole-trapnell-lab.github.io/monocle3/docs/introduction/
# =============================================================================

# ---- 0. 安装 ----------------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install(c("BiocGenerics","DelayedArray","DelayedMatrixStats",
#                        "limma","lme4","lattice","SingleCellExperiment",
#                        "SummarizedExperiment"))
# devtools::install_github("cole-trapnell-lab/leidenbase")
# devtools::install_github("cole-trapnell-lab/monocle3")
# # Seurat 互转：
# remotes::install_github("satijalab/seurat-wrappers")

suppressPackageStartupMessages({
  library(monocle3)
  library(Seurat)
  library(SeuratWrappers)
  library(ggplot2)
  library(dplyr)
})
set.seed(42)

# ---- 0.1 项目目录与 figures/（所有图片必须保存到这里）------------------------
# proj_dir <- "path/to/project"
# fig_dir  <- file.path(proj_dir, "figures")
# dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
# 后续所有 ggsave / pdf 都用 file.path(fig_dir, "xxx.pdf")

# ---- 1. 从 Seurat 对象构建 cell_data_set（CDS） -----------------------------
seurat <- readRDS("seurat_processed.rds")

# 方案 A：SeuratWrappers 直接转（推荐，最省事）
cds <- as.cell_data_set(seurat)
# 若 rowData 缺少 gene_short_name 列：
rowData(cds)$gene_short_name <- rownames(cds)

# 方案 B：手动构建（更可控）
# exprs <- GetAssayData(seurat, layer = "counts")
# cell_meta <- seurat@meta.data
# gene_meta <- data.frame(gene_short_name = rownames(exprs),
#                         row.names       = rownames(exprs))
# cds <- new_cell_data_set(
#   expression_data = exprs,
#   cell_metadata   = cell_meta,
#   gene_metadata   = gene_meta
# )

# ---- 2. 预处理 + 降维 + 聚类 ------------------------------------------------
# num_dim 根据数据复杂度调整，常见 30–100
cds <- preprocess_cds(cds, num_dim = 50, method = "PCA", verbose = TRUE)
plot_pc_variance_explained(cds)

cds <- reduce_dimension(
  cds,
  reduction_method = "UMAP",
  preprocess_method = "PCA",
  umap.n_neighbors  = 15L,
  umap.min_dist     = 0.1
)

# 注意：Monocle3 的 resolution 单位是 1e-4 量级，和 Seurat 不同
cds <- cluster_cells(
  cds,
  reduction_method = "UMAP",
  resolution       = 1e-3,
  k                = 20
)

# 用 Seurat 中已注释的细胞类型覆盖，方便后续选起点
colData(cds)$celltype <- seurat$celltype[colnames(cds)]

p1 <- plot_cells(cds, color_cells_by = "cluster",   label_groups_by_cluster = FALSE,
                 label_leaves = FALSE, label_branch_points = FALSE)
p2 <- plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster = FALSE,
                 label_leaves = FALSE, label_branch_points = FALSE,
                 group_label_size = 4)
print(p1 | p2)

# ---- 3. 学习轨迹主图 --------------------------------------------------------
cds <- learn_graph(
  cds,
  use_partition        = TRUE,
  close_loop           = FALSE,
  learn_graph_control  = list(ncenter = 400, minimal_branch_len = 10)
)

plot_cells(cds,
           color_cells_by      = "celltype",
           label_groups_by_cluster = FALSE,
           label_leaves        = TRUE,
           label_branch_points = TRUE,
           label_roots         = TRUE,
           graph_label_size    = 3,
           group_label_size    = 4)

# ---- 4. 选择根节点并计算 pseudotime -----------------------------------------
# 关键：根节点必须由已知的“起始状态”决定（如干/祖细胞、最早时间点）
# 方式 A：按已知祖细胞 cluster 选根（推荐）
root_type <- "Stem / Progenitor"   # ← 改成你数据里的起始细胞类型名
root_cells <-  colData(cds) %>% as.data.frame() %>%
  filter(celltype == root_type) %>% rownames()
cds <- order_cells(cds, root_cells = root_cells)

# 方式 B：按 principal graph node 选根
# cds <- order_cells(cds, root_pr_nodes = "Y_12")

# 方式 C：交互式选根（RStudio 里会提示点图）
# cds <- order_cells(cds)

# ⚠️ 画图颜色按 §0.1 与用户确认；连续变量（pseudotime）直接给渐变，例如：
#   scale_color_viridis_c(option = "plasma", na.value = "grey90")
#   scale_color_gradient2(low = "#3C5488", mid = "white", high = "#E64B35",
#                         midpoint = median(pseudotime, na.rm=TRUE))
p_pseudo <- plot_cells(cds,
                       color_cells_by      = "pseudotime",
                       label_cell_groups   = FALSE,
                       label_leaves        = TRUE,
                       label_branch_points = TRUE,
                       label_roots         = TRUE,
                       graph_label_size    = 3) +
  scale_color_viridis_c(option = "plasma", na.value = "grey90") +
  ggtitle("Pseudotime")
print(p_pseudo)

# ---- 5. 沿轨迹的基因模块与差异基因 ------------------------------------------
# 5.1 graph_test：Moran's I 自相关，找随轨迹变化的基因
# ⚠️ 筛选前必须按 SKILL.md §0.3 与用户确认：
#   - 是否按 q_value、p_value、morans_I 筛选
#   - 各自阈值（常见组合见下方三选一）
#   - Top N 取多少基因进入下游热图/模块
qval_thresh <- 0.05    # ← 与用户确认
moran_thresh <- NULL   # ← 与用户确认；标准 0.1、严格 0.25，不筛设 NULL
top_n       <- 50      # ← 与用户确认

graph_test_res <- graph_test(
  cds,
  neighbor_graph = "principal_graph",
  cores          = 4
)

# 按用户确认的阈值筛选
traj_de <- graph_test_res
if (!is.null(qval_thresh))  traj_de <- subset(traj_de, q_value  < qval_thresh)
if (!is.null(moran_thresh)) traj_de <- subset(traj_de, morans_I > moran_thresh)

# 取最显著的前 N 个基因
traj_genes <- traj_de %>%
  arrange(q_value, desc(morans_I)) %>%
  head(top_n) %>%
  pull(gene_short_name)

# ⚠️ 画图颜色按 §0.1 与用户确认；pseudotime 连续变量示例：
# pseudo_cols <- viridis::viridis(100, option = "plasma")
plot_genes_in_pseudotime(
  cds[traj_genes[1:6],],
  color_cells_by = "celltype",
  min_expr       = 0.5,
  ncol           = 3   # ← 多图分面列数，按 §0.1 与用户确认
)

# 基因表达在轨迹 UMAP 上的位置
plot_cells(cds, genes = traj_genes[1:4],
           show_trajectory_graph = TRUE,
           label_cell_groups     = FALSE)

# 5.2 回归分析 fit_models：同时控制批次
# ⚠️ 模型公式与系数筛选阈值按 §0.3 与用户确认：
#   - 除 pseudotime 外是否加入 batch / sample / cell_cycle
#   - coefficient_table 用哪个 term、q_value 阈值多少
# fits <- fit_models(cds[rowData(cds)$gene_short_name %in% traj_genes,],
#                    model_formula_str = "~ pseudotime + batch",
#                    expression_family = "quasipoisson")
# coef_tab <- coefficient_table(fits) %>%
#   filter(term == "pseudotime", q_value < 0.05)
# head(coef_tab)

# 5.3 模块聚类：找共调控基因模块
genemod <- find_gene_modules(
  cds[traj_genes,],
  resolution = 1e-2,
  cores      = 4
)
cell_group <- tibble::tibble(
  cell   = rownames(colData(cds)),
  group  = colData(cds)$celltype
)
agg_mat <- aggregate_gene_expression(cds, genemod, cell_group)
# ⚠️ 热图颜色向量按 §0.1 与用户确认：
# pheatmap::pheatmap(agg_mat, scale = "column", clustering_method = "ward.D2",
#                    color = colorRampPalette(c("#3C5488","white","#E64B35"))(100))
pheatmap::pheatmap(agg_mat, scale = "column", clustering_method = "ward.D2")

# ---- 6. 分支分析（branch-dependent genes）-----------------------------------
# 在含分叉的轨迹上，交互选分支：
# branches <- choose_graph_segments(cds)
# 然后对该分支内细胞做 graph_test
# cds_sub <- choose_cells(cds)
# branch_test <- graph_test(cds_sub, neighbor_graph = "principal_graph", cores=4)

# 若需要 Monocle2 风格的 BEAM 分支检验，可参考 03_monocle2_BEAM.R
# （Monocle3 已用 graph_test + 选分支替代大部分 BEAM 用法）

# ---- 7. 保存 ----------------------------------------------------------------
saveRDS(cds, "cds_monocle3.rds")
ggsave("monocle3_pseudotime.pdf", p_pseudo, width = 7, height = 6)

sessionInfo()
