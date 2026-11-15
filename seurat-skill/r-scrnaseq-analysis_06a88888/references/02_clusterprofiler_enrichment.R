# =============================================================================
# 02_clusterprofiler_enrichment.R
# clusterProfiler 4.x + MSigDB/fgsea：GO/KEGG ORA 与 GSEA 模板
# 官方文档: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
# =============================================================================

# ---- 0. 安装 ----------------------------------------------------------------
# if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "enrichplot", "DOSE",
#                        "org.Hs.eg.db", "org.Mm.eg.db",
#                        "ReactomePA", "msigdbr", "fgsea"))
# install.packages(c("ggplot2", "dplyr", "tidyr"))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
  library(ReactomePA)
  library(msigdbr)
  library(fgsea)
  library(org.Hs.eg.db)   # 人；小鼠换成 org.Mm.eg.db
  library(ggplot2)
  library(dplyr)
})

# 物种参数（小鼠改：OrgDb=org.Mm.eg.db, kegg="mmu", species="Mus musculus"）
ORG_DB       <- org.Hs.eg.db
KEGG_ORG     <- "hsa"
MSIGDB_SPP   <- "Homo sapiens"

# ⚠️ 本模板涉及大量富集可视化（dotplot/barplot/cnetplot/emapplot/gseaplot2/
#    ridgeplot/heatplot）。写实际分析代码前，必须按 SKILL.md §0.1 与用户确认：
#      - showCategory（每张图展示前多少通路/基因集）
#      - 多图分面（ncol/nrow/facet_grid 变量）
#      - 配色（直接给颜色向量，不要只报调色板名）
#      - PDF 设备 width/height（dotplot 通路多时 height 取 7-10）
#    所有 ggsave 必须保存到项目目录下的 figures/ 子目录：
#      proj_dir <- "path/to/project"
#      fig_dir  <- file.path(proj_dir, "figures"); dir.create(fig_dir, recursive=TRUE)
#      ggsave(file.path(fig_dir, "xxx.pdf"), ...)
#    另外，输入的 deg 来自 FindMarkers/FindAllMarkers，其 only.pos /
#    logfc.threshold / p_val_adj 阈值必须按 §0.2 与用户确认后再跑差异。

# ---- 1. 从 Seurat FindAllMarkers/FindMarkers 结果读入差异基因 ----------------
# markers <- read.csv("cluster_markers.csv")          # FindAllMarkers 输出
# deg <- subset(markers, cluster == "3")               # 取某一簇
# deg <- read.csv("condition_vs_control_deg.csv")       # 或 FindMarkers 输出

# 为便于演示，这里构造一个示例 deg（实际使用时替换为真实数据）
# deg 必须包含：gene (SYMBOL)、avg_log2FC、p_val_adj（或 pvalue）
set.seed(1)
deg <- data.frame(
  gene     = c("TP53","BRCA1","MYC","EGFR","STAT3","IL6","TNF","VEGFA",
               "CD4","CD8A","MS4A1","CD14","NCAM1","EPCAM","KRT19","LGR5"),
  avg_log2FC = c(2.1, 1.8, -1.5, 2.4, 1.2, 0.9, -0.8, 1.5,
                 -2.0, -1.7, 0.6, 1.1, -0.5, 1.9, 2.2, -1.1),
  p_val_adj = c(1e-8, 2e-6, 3e-5, 1e-10, 4e-4, 5e-3, 2e-2, 8e-5,
                1e-7, 3e-6, 2e-3, 6e-4, 3e-2, 5e-7, 2e-9, 7e-3)
)

# ---- 2. SYMBOL -> ENTREZID 转换 ---------------------------------------------
gene_df <- bitr(deg$gene, fromType = "SYMBOL", toType = "ENTREZID",
                OrgDb = ORG_DB)
cat("输入基因:", length(deg$gene), " | 成功映射:", nrow(gene_df), "\n")

deg <- merge(deg, gene_df, by.x = "gene", by.y = "SYMBOL")

# ---- 3. ORA：过表达分析 ----------------------------------------------------
# 3.1 取显著上调基因（典型：|log2FC|>1, p.adj<0.05；用户可放宽）
sig_up <- deg %>% filter(avg_log2FC > 1 & p_val_adj < 0.05) %>% pull(ENTREZID)
sig_dn <- deg %>% filter(avg_log2FC < -1 & p_val_adj < 0.05) %>% pull(ENTREZID)

# 3.2 GO 富集（BP/MF/CC/ALL）
ego <- enrichGO(
  gene          = sig_up,
  OrgDb         = ORG_DB,
  keyType       = "ENTREZID",
  ont           = "BP",                # "BP" / "MF" / "CC" / "ALL"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE                 # 自动把 ENTREZID 转成 SYMBOL
)
head(ego, 20)
dotplot(ego, showCategory = 20) + ggtitle("GO-BP ORA")
barplot(ego, showCategory = 20)
cnetplot(ego, showCategory = 10, foldChange = setNames(deg$avg_log2FC, deg$SYMBOL))
emapplot(pairwise_termsim(ego), showCategory = 20)

# 3.3 KEGG 富集
ekk <- enrichKEGG(
  gene          = sig_up,
  organism      = KEGG_ORG,
  keyType       = "kegg",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)
ekk <- setReadable(ekk, OrgDb = ORG_DB, keyType = "ENTREZID")
dotplot(ekk, showCategory = 20) + ggtitle("KEGG ORA")
# 通路图（需联网，pathview）
# pathview::pathview(gene.data = sig_up, pathway.id = "hsa04110", species = KEGG_ORG)

# 3.4 Reactome 富集
er <- enrichPathway(gene = sig_up, organism = "human", pvalueCutoff = 0.05,
                    readable = TRUE)
dotplot(er, showCategory = 20)

# 3.5 多组/多簇比较：compareCluster
# markers$ENTREZID <- mapIds(ORG_DB, markers$gene, "ENTREZID", "SYMBOL")
# formula 形式：ENTREZID ~ cluster
# ck <- compareCluster(ENTREZID ~ cluster, data = markers, fun = "enrichGO",
#                      OrgDb = ORG_DB, ont = "BP", pAdjustMethod = "BH")
# dotplot(ck, x = ~Cluster, showCategory = 10)

# ---- 4. GSEA：基因集富集分析 -----------------------------------------------
# 4.1 构造排序基因列表：必须包含所有检测基因，按统计量（log2FC 或 signed p）递减
geneList <- sort(setNames(deg$avg_log2FC, deg$ENTREZID), decreasing = TRUE)
# 去掉 NA 和重复
geneList <- geneList[!is.na(geneList) & !duplicated(names(geneList))]
stopifnot(!is.unsorted(rev(geneList)))  # 必须从大到小

# 4.2 基于 GO 的 GSEA
gsea_go <- gseGO(
  geneList      = geneList,
  OrgDb         = ORG_DB,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  eps           = 0
)
gsea_go <- setReadable(gsea_go, ORG_DB, "ENTREZID")
gseaplot2(gsea_go, geneSetID = 1:3, pvalue_table = TRUE)
ridgeplot(gsea_go, showCategory = 20)

# 4.3 基于 KEGG 的 GSEA
gsea_kegg <- gseKEGG(
  geneList      = geneList,
  organism      = KEGG_ORG,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  eps           = 0
)
gsea_kegg <- setReadable(gsea_kegg, ORG_DB, "ENTREZID")
gseaplot2(gsea_kegg, geneSetID = 1:3)

# ---- 5. MSigDB 基因集 + fgsea / clusterProfiler::GSEA ----------------------
# 5.1 取 Hallmark (H) 集合；其他常用：C2(CP:KEGG/REACTOME), C5(GO), C7(immune)
msig_h <- msigdbr(species = MSIGDB_SPP, collection = "H")
term2gene_h <- msig_h %>% select(gs_name, gene_symbol)

# clusterProfiler 风格（输入 SYMBOL 即可，无需 ENTREZ）
geneList_sym <- sort(setNames(deg$avg_log2FC, deg$gene), decreasing = TRUE)
gsea_h <- GSEA(
  geneList      = geneList_sym,
  TERM2GENE     = term2gene_h,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  eps           = 0
)
gseaplot2(gsea_h, geneSetID = 1:5)
dotplot(gsea_h, split = ".sign", font.size = 9) + facet_grid(~.sign)

# 5.2 fgsea 快速版（同样输入 MSigDB 列表 + 排序向量）
pathways_h <- split(msig_h$gene_symbol, msig_h$gs_name)
fgsea_res <- fgsea(
  pathways = pathways_h,
  stats    = geneList_sym,
  minSize  = 15,
  maxSize  = 500,
  eps      = 0
)
data.table::fwrite(fgsea_res, "fgsea_hallmark.csv")
plotEnrichment(pathways_h[["HALLMARK_INFLAMMATORY_RESPONSE"]], geneList_sym)

# ---- 6. 结果保存 ------------------------------------------------------------
write.csv(as.data.frame(ego),      "enrichGO_BP.csv",      row.names = FALSE)
write.csv(as.data.frame(ekk),      "enrichKEGG.csv",       row.names = FALSE)
write.csv(as.data.frame(gsea_go),  "gseaGO_BP.csv",        row.names = FALSE)
write.csv(as.data.frame(gsea_kegg),"gseaKEGG.csv",         row.names = FALSE)
saveRDS(list(ora_go=ego, ora_kegg=ekk, gsea_go=gsea_go, gsea_kegg=gsea_kegg),
        "enrichment_results.rds")

sessionInfo()
