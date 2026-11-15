# =============================================================================
# 01_seurat_qc_clustering.R
# Seurat v5/v6 兼容的标准 scRNA-seq 分析流程模板
# 官方文档: https://satijalab.org/seurat/articles/pbmc3k_tutorial
# DoubletFinder: https://github.com/chris-mcginnis-ucsf/DoubletFinder
#
# “帮我完成一个单细胞转录组基本分析”请使用
# references/05_seurat_v6_basic_pipeline.R，而不是本通用模板。
# =============================================================================

# ---- 0. 安装依赖（首次运行时取消注释）---------------------------------------
# install.packages("pak")
# pak::pkg_install(c("Seurat", "SeuratObject", "patchwork", "ggplot2",
#                    "dplyr", "hdf5r", "sctransform", "harmony"))
# BiocManager::install(c("SingleR", "celldex"))
# # DoubletFinder（需要时安装）：
# remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", force = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

set.seed(42)
options(future.globals.maxSize = 8 * 1024^3)

# ---- 0.1 项目目录与 figures/ 文件夹（强制）----------------------------------
# ⚠️ proj_dir 必须替换为用户指定的项目文件夹或原始数据所在文件夹。
#    所有图片保存到 fig_dir，数据文件保存到 proj_dir，不要混放。
proj_dir <- "path/to/project"                      # ← 改成用户的项目路径
fig_dir  <- file.path(proj_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
cat("图片将保存到:", fig_dir, "\n")

# ---- 1. 读入数据 ------------------------------------------------------------
# 10x CellRanger 输出目录（包含 matrix.mtx/barcodes.tsv/features.tsv）
data_dir <- file.path(proj_dir, "filtered_feature_bc_matrix")  # 按实际路径调整
counts   <- Read10X(data.dir = data_dir)

# 若是单个 .h5 文件：
# counts <- Read10X_h5(file.path(proj_dir, "filtered_feature_bc_matrix.h5"))

seurat <- CreateSeuratObject(
  counts       = counts,
  project      = "scRNA_demo",
  min.cells    = 3,    # 基因至少在多少个细胞中表达才保留
  min.features = 200   # 细胞至少检测到多少个基因才保留
)

# ---- 2. 计算 QC 指标 --------------------------------------------------------
# ⚠️ 物种由用户确认；人用 "^MT-" / "^RP[SL]"，小鼠用 "^mt-" / "^Rp[sl]"
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
seurat[["percent.rb"]] <- PercentageFeatureSet(seurat, pattern = "^RP[SL]")

# ---- 2.1 先出 QC 分布图，保存到 figures/，然后询问阈值（强制）---------------
# ⚠️ 在看到这些图之前，不要直接 subset！
#    把图片路径给用户看，然后询问：
#      - nFeature_RNA 下限/上限（根据小提琴长尾和散点图定）
#      - nCount_RNA 下限/上限
#      - percent.mt 上限（典型 15%，高代谢组织可放宽到 20%）
#      - 是否需要按 percent.rb 过滤
#    用户确认后再填到下面的 subset 条件里。
p_vln <- VlnPlot(seurat,
                 features = c("nFeature_RNA", "nCount_RNA",
                              "percent.mt", "percent.rb"),
                 ncol = 4, pt.size = 0.1)
ggsave(file.path(fig_dir, "QC_violin.pdf"), p_vln,
       width = 12, height = 5, device = cairo_pdf)

p_sc1 <- FeatureScatter(seurat, feature1 = "nCount_RNA",
                        feature2 = "percent.mt", pt.size = 0.3)
p_sc2 <- FeatureScatter(seurat, feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA", pt.size = 0.3)
ggsave(file.path(fig_dir, "QC_scatter_nCount_vs_mt.pdf"), p_sc1,
       width = 6, height = 5, device = cairo_pdf)
ggsave(file.path(fig_dir, "QC_scatter_nCount_vs_nFeature.pdf"), p_sc2,
       width = 6, height = 5, device = cairo_pdf)

cat("QC 图已保存到", fig_dir, "，请查看后确认阈值。\n")

# 下面阈值是占位示例，必须替换为用户看图后确认的值：
nfeat_lower <- 200     # ← 与用户确认
nfeat_upper <- 6000    # ← 与用户确认
ncount_lower <- 500    # ← 与用户确认
ncount_upper <- Inf    # ← 与用户确认（不需要上限设 Inf）
mt_upper    <- 15      # ← 与用户确认（百分比）

seurat <- subset(
  seurat,
  subset = nFeature_RNA > nfeat_lower & nFeature_RNA < nfeat_upper &
           nCount_RNA   > ncount_lower & nCount_RNA   < ncount_upper &
           percent.mt   < mt_upper
)

# ---- 3. 标准化（默认 LogNormalize；用户明确要求才用 SCTransform）------------
# 方案 A（默认）：LogNormalize
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize",
                        scale.factor = 1e4)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(seurat), 20)
p_hvg <- VariableFeaturePlot(seurat)
p_hvg <- LabelPoints(plot = p_hvg, points = top20, repel = TRUE)
ggsave(file.path(fig_dir, "HVG_variablefeatures.pdf"), p_hvg,
       width = 7, height = 6, device = cairo_pdf)

# 如需回归细胞周期 / 线粒体比例：
# s.genes   <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# seurat    <- CellCycleScoring(seurat, s.features=s.genes, g2m.features=g2m.genes)
# seurat    <- ScaleData(seurat, vars.to.regress=c("percent.mt","S.Score","G2M.Score"))
seurat <- ScaleData(seurat)

# 方案 B（仅在用户明确要求时启用）：SCTransform (vst.flavor="v2")
# seurat <- SCTransform(seurat, vst.flavor = "v2",
#                       vars.to.regress = "percent.mt", verbose = FALSE)
# 后续 PCA/聚类在 SCT assay 上进行。
use_sct <- FALSE   # ← 用户明确要 SCT 时改成 TRUE

# ---- 4. 线性降维：PCA -------------------------------------------------------
seurat <- RunPCA(seurat, features = VariableFeatures(seurat), npcs = 50)
p_elbow <- ElbowPlot(seurat, ndims = 50)
ggsave(file.path(fig_dir, "PCA_elbow.pdf"), p_elbow,
       width = 6, height = 4, device = cairo_pdf)

dims_use <- 1:30  # 根据 ElbowPlot 调整，常见 15–40；与用户确认

# ---- 4.5 DoubletFinder 双细胞筛查（询问用户后启用）--------------------------
# ⚠️ 先询问用户是否需要做 DoubletFinder。注意：
#    - 适用于单个样本，多样本 merge/整合后不要直接跑。
#    - 下面代码默认 use_sct=FALSE（LogNormalize）；若用 SCT，把 sct 参数改 TRUE。
run_doubletfinder <- FALSE   # ← 与用户确认（TRUE/FALSE）

if (run_doubletfinder) {
  library(DoubletFinder)

  # 4.5.1 pK 参数扫描
  sweep.res <- paramSweep(seurat, PCs = dims_use, sct = use_sct)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  p_bcmvn <- ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
    geom_point() + geom_line() + theme_bw()
  ggsave(file.path(fig_dir, "DoubletFinder_BCmvn.pdf"), p_bcmvn,
         width = 7, height = 4, device = cairo_pdf)

  # 取 BCmvn 最大的 pK（也可以让用户看图后指定）
  pK_val <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  cat("自动选择 pK =", pK_val, "\n")

  # 4.5.2 估计预期 doublet 数
  # ⚠️ doublet_rate 根据 10x 上样细胞数查官方表格，典型 0.4%–8%，
  #    教程默认 7.5%。必须与用户确认。
  doublet_rate <- 0.075   # ← 与用户确认
  nExp_poi <- round(doublet_rate * ncol(seurat))

  # 若已有注释/聚类结果，可用 modelHomotypic 校正同型 doublet：
  # annotations <- seurat$seurat_clusters
  # homotypic.prop <- modelHomotypic(annotations)
  # nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  nExp_poi.adj <- nExp_poi   # 无注释时不校正

  # 4.5.3 跑两次：Poisson 估计 + 同型校正
  seurat <- doubletFinder(seurat, PCs = dims_use, pN = 0.25,
                          pK = pK_val, nExp = nExp_poi,
                          reuse.pANN = FALSE, sct = use_sct)
  pann_col <- grep("^pANN", colnames(seurat@meta.data), value = TRUE)[1]
  seurat <- doubletFinder(seurat, PCs = dims_use, pN = 0.25,
                          pK = pK_val, nExp = nExp_poi.adj,
                          reuse.pANN = pann_col, sct = use_sct)

  df_class_col <- grep("^DF.classifications", colnames(seurat@meta.data),
                       value = TRUE)[1]
  cat("DoubletFinder 分类列:", df_class_col, "\n")
  print(table(seurat[[df_class_col]]))

  # 先跑 UMAP 看 doublet 分布（颜色按 §0.1 与用户确认，此处用默认）
  seurat <- RunUMAP(seurat, dims = dims_use, reduction = "pca")
  p_df <- DimPlot(seurat, group.by = df_class_col, raster = FALSE) +
    ggtitle("DoubletFinder")
  ggsave(file.path(fig_dir, "DoubletFinder_UMAP.pdf"), p_df,
         width = 6, height = 5, device = cairo_pdf)

  # 过滤 doublet
  n_before <- ncol(seurat)
  seurat <- subset(seurat, cells = colnames(seurat)[seurat[[df_class_col]][,1] == "Singlet"])
  cat("Doublet 过滤：", n_before, "->", ncol(seurat), "cells\n")
}

# ---- 5. 细胞聚类与非线性降维 ------------------------------------------------
seurat <- FindNeighbors(seurat, dims = dims_use)
seurat <- FindClusters(seurat, resolution = 0.6)  # 0.4–1.2 间调整
seurat <- RunUMAP(seurat, dims = dims_use, reduction = "pca")
# seurat <- RunTSNE(seurat, dims = dims_use)

# ⚠️ 画图前必须按 SKILL.md §0.1 与用户确认配色、分面、尺寸等；
#    颜色直接给向量，例如：
# cluster_colors <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F",
#                     "#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85")
p_umap <- DimPlot(seurat, label = TRUE, repel = TRUE) +
  ggtitle(sprintf("UMAP | %d cells, %d clusters",
                  ncol(seurat), length(levels(Idents(seurat)))))
ggsave(file.path(fig_dir, "UMAP_clusters.pdf"), p_umap,
       width = 7, height = 6, device = cairo_pdf)

# ---- 6. 簇间差异基因（marker） ----------------------------------------------
# ⚠️ 调用 FindAllMarkers / FindMarkers 前必须按 SKILL.md §0.2 与用户确认：
#   - only.pos (TRUE/FALSE)
#   - logfc.threshold（默认 0.25）
#   - 结果筛选：p_val_adj / p_val 阈值（默认 0.05），是否再按 avg_log2FC 过滤
only_pos    <- TRUE     # ← 与用户确认
lfc_thresh  <- 0.25    # ← 与用户确认
padj_thresh <- 0.05    # ← 与用户确认

markers_all <- FindAllMarkers(
  seurat,
  only.pos        = only_pos,
  min.pct         = 0.25,
  logfc.threshold = lfc_thresh,
  test.use        = "wilcox"
)
# markers_all <- subset(markers_all, p_val_adj < padj_thresh)
write.csv(markers_all, file.path(proj_dir, "cluster_markers.csv"),
          row.names = FALSE)

# 每个 cluster 取 top N（N 按 §0.1 与用户确认）
top_n <- 10
top_markers <- markers_all %>% group_by(cluster) %>%
  slice_max(n = top_n, order_by = avg_log2FC)

# ⚠️ 热图颜色向量按 §0.1 确认，例如：
# heatmap_cols <- circlize::colorRamp2(c(-2,0,2), c("#3C5488","white","#E64B35"))
p_hm <- DoHeatmap(seurat, features = top_markers$gene) + NoLegend()
ggsave(file.path(fig_dir, "marker_heatmap_top10.pdf"), p_hm,
       width = 12, height = 10, device = cairo_pdf)

# 两 cluster 之间差异
de_0_vs_1 <- FindMarkers(seurat, ident.1 = "0", ident.2 = "1",
                         only.pos = only_pos,
                         min.pct = 0.25, logfc.threshold = lfc_thresh)

# ---- 7. 细胞类型注释 --------------------------------------------------------
# 方案 A：基于已知 marker 手工注释（最可靠）
# new_ids <- c("0" = "...", "1" = "...")  # 按实际 top markers 填
# seurat <- RenameIdents(seurat, new_ids)
# seurat$celltype <- Idents(seurat)
#
# canonical_markers <- c(...)
# p_dot <- DotPlot(seurat, features = canonical_markers) + RotatedAxis()
# ggsave(file.path(fig_dir, "dotplot_canonical_markers.pdf"), p_dot,
#        width = 8, height = 5, device = cairo_pdf)
# p_fp <- FeaturePlot(seurat, features = c("...","..."))
# ggsave(file.path(fig_dir, "featureplot_markers.pdf"), p_fp,
#        width = 10, height = 8, device = cairo_pdf)

# 方案 B：SingleR 自动注释
# ref <- celldex::HumanPrimaryCellAtlasData()
# pred <- SingleR::SingleR(test = GetAssayData(seurat, layer="data"),
#                          ref = ref, labels = ref$label.main)
# seurat$singler <- pred$labels
#
# 人源数据的 Seurat v6 自动注释见 references/06_azimuth_scslide_v6.R：
# seurat <- NormalizeData(seurat, assay = "RNA")
# seurat <- AzimuthAPI::CloudAzimuth(seurat)

# ---- 8. 多样本整合（可选，参考） --------------------------------------------
# obj.list <- SplitObject(seurat, split.by = "orig.ident")
# obj.list <- lapply(obj.list, SCTransform, vst.flavor = "v2")
# features <- SelectIntegrationFeatures(obj.list, nfeatures = 3000)
# obj.list <- PrepSCTIntegration(obj.list, anchor.features = features)
# anchors  <- FindIntegrationAnchors(obj.list, normalization.method = "SCT",
#                                    anchor.features = features)
# combined <- IntegrateData(anchors, normalization.method = "SCT")
# # Seurat v5/v6 推荐：IntegrateLayers(method=CCAIntegration|HarmonyIntegration)

# ---- 9. 保存 ----------------------------------------------------------------
saveRDS(seurat, file.path(proj_dir, "seurat_processed.rds"))
# 如用户要求写回 RData：
# save(seurat, file = file.path(proj_dir, "basal_analysis.RData"))

sessionInfo()
