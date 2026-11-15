# =============================================================================
# 04_cellchat_communication.R
# CellChat 细胞间通讯分析模板
# 来源: 新格元《单细胞数据科学实战》CellChat 章节 + 官方文档
#       https://github.com/sqjin/CellChat
# =============================================================================

# ---- 0. 安装 ----------------------------------------------------------------
# if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes")
# remotes::install_github("sqjin/CellChat")
# install.packages(c("patchwork","ggplot2","dplyr","Matrix","NMF"))

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
})
options(stringsAsFactors = FALSE)

# ---- 0.1 项目目录与 figures/（所有图片必须保存到这里）------------------------
# proj_dir <- "path/to/project"
# fig_dir  <- file.path(proj_dir, "figures")
# dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
# CellChat 的 pdf/ggsave 都用 file.path(fig_dir, "xxx.pdf")

# ⚠️ CellChat 大量出图（netVisual_circle/aggregate/chord/bubble/individual、
#    netAnalysis_contribution）。写实际代码前必须按 SKILL.md §0.1 与用户确认：
#      - 展示多少条通路 / 多少对 LR
#      - 多图 layout（hierarchy/circle/chord）、ncol/nrow/patchwork 分面
#      - 配色：直接给颜色向量，例如
#          group_cols <- c("#E64B35","#4DBBD5","#00A087","#3C5488",...)
#          并通过 netVisual_* 的 color.use/group.cell.colors 等参数传入
#      - PDF width/height（chord 图建议正方形 7x7–9x9）

# ---- 1. 从 Seurat 构造 CellChat 对象 ----------------------------------------
seurat <- readRDS("seurat_processed.rds")

# CellChat 要求 normalized data matrix（Seurat 的 "data" layer）
# v5：LayerData；旧版：GetAssayData(seurat, slot="data")
data.input <- LayerData(seurat, layer = "data")

meta <- data.frame(
  labels  = seurat$celltype,          # 细胞类型注释
  row.names = colnames(seurat)
)

# 按条件筛选（如只看某个 group）：
# use_cells <- rownames(meta)[meta$labels %in% c("T cells","B cells","Monocytes")]
# data.input <- data.input[, use_cells]
# meta       <- meta[use_cells, , drop=FALSE]

cellchat <- createCellChat(object = data.input,
                           meta   = meta,
                           group.by = "labels")

# ---- 2. 载入配受体数据库 ----------------------------------------------------
CellChatDB <- CellChatDB.human     # 小鼠数据用 CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

# 可只取部分分类："Secreted Signaling" / "ECM-Receptor" / "Cell-Cell Contact"
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
# 也可直接用全库：CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

# ---- 3. 预处理 + 计算通讯概率 ------------------------------------------------
cellchat <- subsetData(cellchat, features = NULL)
future::plan("multisession", workers = 4)   # 并行，根据机器改

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# PPI 投影（人/鼠选对应的 PPI.human / PPI.mouse）
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# 把通讯细胞数过少的对去掉
cellchat <- filterCommunication(cellchat, min.cells = 10)

# 细胞水平的通讯结果（数据框）
df.net <- subsetCommunication(cellchat)
head(df.net)

# 通路水平的通讯结果
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")

future::plan("sequential")

# ---- 4. 网络可视化 ----------------------------------------------------------
groupSize <- as.numeric(table(cellchat@idents))

# 4.1 全局圆圈图：相互作用数量 / 总权重
par(mfrow = c(1,3), xpd = TRUE)
netVisual_circle(cellchat@net$count,
                 vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight,
                 vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, title.name = "Interaction weights/strength")
netVisual_circle(cellchat@net$weight,
                 vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, title.name = "Targets: cDC2",
                 targets.use = "cDC2")

# 4.2 信号通路层级图 / 圆圈图 / 弦图
pathways.show <- cellchat@netP$pathways
head(pathways.show, 20)

# Hierarchy
netVisual_aggregate(cellchat, signaling = pathways.show[1],
                    layout = "hierarchy",
                    vertex.receiver = seq_along(levels(cellchat@idents)))
# Circle
netVisual_aggregate(cellchat, signaling = pathways.show[1], layout = "circle")
# Chord
netVisual_aggregate(cellchat, signaling = pathways.show[1], layout = "chord")

# 4.3 每个细胞类型作为 sender 的独立圆圈图（循环）
# mat <- cellchat@net$weight
# par(mfrow = c(3,4), xpd = TRUE)
# for (i in seq_len(nrow(mat))) {
#   mat2 <- matrix(0, nrow=nrow(mat), ncol=ncol(mat), dimnames=dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight=groupSize, weight.scale=TRUE,
#                    edge.weight.max=max(mat), title.name=rownames(mat)[i])
# }

# ---- 5. 配受体对水平可视化 --------------------------------------------------
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show[1])
LR.show <- pairLR[1, ]   # 取该通路第一对显著配受体

# Hierarchy / circle / chord
netVisual_individual(cellchat, signaling = pathways.show,
                     pairLR.use = LR.show, layout = "hierarchy")
netVisual_individual(cellchat, signaling = pathways.show,
                     pairLR.use = LR.show, layout = "circle")
netVisual_individual(cellchat, signaling = pathways.show,
                     pairLR.use = LR.show, layout = "chord")

# 贡献条形图
p1 <- netAnalysis_contribution(cellchat, signaling = pathways.show[1])
print(p1)

# ---- 6. 气泡图：指定来源 / 接收者 / 通路 ------------------------------------
levels(cellchat@idents)   # 先看分组顺序，再用下标或名字

b1 <- netVisual_bubble(cellchat,
                       sources.use = 4,
                       targets.use = c(5:11),
                       remove.isolate = FALSE)

b2 <- netVisual_bubble(cellchat,
                       sources.use = "CD4+ T cells",
                       targets.use = c("B cells","Monocytes"),
                       signaling  = c("CCL","CXCL","MIF"),
                       remove.isolate = FALSE)
print(b1); print(b2)

# ---- 7. 多组对比（如 control vs treatment） ---------------------------------
# 若有多个条件，分别为每个条件构建 CellChat 对象后：
# cellchat_list <- list(CTRL = cellchat_ctrl, STIM = cellchat_stim)
# cellchat_merge <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
# compareInteractions(cellchat_merge, show.legend = FALSE)
# netVisual_diffInteraction(cellchat_merge, weight.scale = TRUE)
# rankNet(cellchat_merge, mode = "comparison", stacked = TRUE)
# 详见 CellChat 多数据集分析 vignette。

# ---- 8. 保存结果 ------------------------------------------------------------
saveRDS(cellchat, "cellchat_result.rds")
write.csv(df.net,  "cellchat_LR_pairs.csv",  row.names = FALSE)
write.csv(df.netp, "cellchat_pathways.csv",  row.names = FALSE)

# 批量保存每个通路的 chord 图：
# dir.create("cellchat_chord", showWarnings = FALSE)
# for (pw in pathways.show) {
#   pdf(sprintf("cellchat_chord/%s.pdf", pw), width=7, height=7)
#   try(netVisual_aggregate(cellchat, signaling=pw, layout="chord"))
#   dev.off()
# }

sessionInfo()
