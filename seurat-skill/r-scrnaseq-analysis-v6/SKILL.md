---
name: r-scrnaseq-analysis-v6
description: 使用 R 语言完成单细胞转录组测序（scRNA-seq）数据分析。覆盖 Seurat v6 标准流程、SCTransform、Leiden 聚类、AzimuthAPI 自动注释、scSLIDE 样本水平分析、clusterProfiler GO/KEGG/Reactome ORA 与 GSEA、Monocle3/Monocle2 轨迹分析和 CellChat 细胞通讯。当用户要求写 R 代码做单细胞数据分析、Seurat 流程、细胞类型注释、差异表达、富集、轨迹或细胞通讯，或说“帮我完成一个单细胞转录组基本分析”时触发。
---

# R 单细胞转录组测序分析工具

帮助用户用 R 完成单细胞 RNA-seq（scRNA-seq）的完整分析：从 10x 数据读入、质控、降维聚类、细胞类型注释，到差异基因功能富集、GSEA、拟时序轨迹、细胞通讯。输出可直接运行的 R 脚本，代码注释清晰，关键参数给出选择依据。

## 核心工具栈（官方来源）

- **Seurat v6**：质控、SCTransform、整合、PCA/UMAP、Leiden 聚类和可视化；v6 核心流程与 v5 向后兼容。文档 https://satijalab.org/seurat/ ，版本说明 https://satijalab.org/seurat/articles/announcements ，changelog https://satijalab.org/seurat/news/index.html ，命令速查 https://satijalab.org/seurat/articles/essential_commands
- **AzimuthAPI（Seurat v6 的 Pan-human Azimuth）**：云端 `CloudAzimuth()` 或本地 `ANNotate()`，返回层级标签、置信度和 128 维 `azimuth_embed`。文档 https://satijalab.org/seurat/articles/intro_azimuthapi_vignette
- **scSLIDE（Seurat v6）**：从多样本单细胞数据构造 landmark-by-sample 密度矩阵，用于样本聚类、样本轨迹和轨迹相关差异表达。文档 https://satijalab.org/seurat/articles/scSLIDE_vignette
- **clusterProfiler 4.x（Bioconductor）**：ORA 与 GSEA，支持 GO/KEGG/Reactome 等数千物种，多组比较与 LLM 解读。文档 https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html ，正式方法见 Wu et al. 2021 The Innovation；最新多组学方案见 Xu et al. 2024 Nat Protoc。
- **GSEA / MSigDB**：基因集富集分析方法与基因集数据库（H/C2/C5 等集合，2026.1 已更新到 Ensembl 115）。https://www.gsea-msigdb.org/gsea/index.jsp 。R 中用 `clusterProfiler::GSEA()` 或 `fgsea::fgsea()` 跑，基因集用 `msigdbr` 包获取。
- **Monocle3**：UMAP 流形学习轨迹、聚类、图自相关差异分析、拟时序。https://cole-trapnell-lab.github.io/monocle3/
- **Monocle 2**：DDRTree 反向图嵌入、BEAM 分支依赖表达分析、Census 计数转换。https://cole-trapnell-lab.github.io/monocle-release/
- **CellChat**：细胞间通讯网络推断（配体-受体）。https://github.com/sqjin/CellChat

## 固定触发：单细胞转录组基本分析（Seurat v6 自动档）

当用户明确说“帮我完成一个单细胞转录组基本分析”或语义等价地要求按本 skill 的“基本分析自动档”执行时，采用 `references/05_seurat_v6_basic_pipeline.md`，直接完成全流程。该短语视为用户已经授权下列固定参数，不再按 §0.1 询问图片尺寸、分面、配色，也不询问 PC 数、聚类分辨率或是否使用 SCTransform。

固定流程：

1. **读取与识别输入**：优先定位用户给出的 10x `filtered_feature_bc_matrix` 目录、`Read10X_h5()` 可读的 `.h5` 文件；若用户提供已有 Seurat RDS，可直接使用。若项目内有多个候选输入且无法唯一判断，先列出候选并让用户指定，不能猜。
2. **判断物种**：根据基因名同时检查人/鼠线粒体与核糖体基因写法及常见管家基因写法。人类使用 `^MT-`、`^RP[SL]`，小鼠使用 `^mt-`、`^Rp[sl]`。证据不足或冲突时停止并报告，不要猜物种。
3. **基础 QC**：计算并保存 `nFeature_RNA`、`nCount_RNA`、`percent.mt`、`percent.rb` 图。自动采用 `nFeature_RNA >= 200`、`nCount_RNA >= 500`；`percent.mt` 上限为 `min(30, max(10, median + 3*MAD))`，`percent.rb` 上限为 `max(5, median + 3*MAD)`，并把实际阈值写入清单。若自动阈值会移除超过 30% 细胞，保留 QC 图和原始对象后停止，提示先诊断数据，不要静默放宽阈值。
4. **SCTransform 归一化与回归**：固定 `SCTransform(..., assay="RNA", new.assay.name="SCT", vst.flavor="v2", vars.to.regress=c("percent.mt","percent.rb"))`，在 SCT assay 上继续。`vars_to_regress` 是本 skill 的固定触发档参数名；Seurat 实际参数名是 `vars.to.regress`。
5. **PCA 与 ndims**：在 SCT assay 上 `RunPCA()`，保存肘图；从 PC 标准差曲线的最大曲率自动识别转折点，`ndims` 取该点并按可用 PC 数截断。把识别结果和曲线图同时交付。
6. **UMAP 与 Leiden**：`FindNeighbors(dims=1:ndims)`，`RunUMAP(dims=1:ndims)`，然后 `FindClusters(algorithm=4, resolution=0.5, cluster.name="leiden_res0.5")`。Leiden 固定使用 `resolution=0.5`，不要擅自改成 Louvain。
7. **Azimuth 自动注释**：先在 RNA assay 上执行 `NormalizeData()`，再调用 `AzimuthAPI::CloudAzimuth()`；云端不可用时，仅当本地 `panhumanpy` 环境已配置才尝试 `AzimuthAPI::ANNotate()`。保存 `azimuth_medium`、`azimuth_fine`/`azimuth_fine_filtered` 和 `final_level_confidence` 图。**Pan-human Azimuth 只能注释人源数据；小鼠或其他物种不得套用该模型。** 非人物种应明确说明限制，并改用手工 marker、SingleR 或用户指定的物种适配参考，其他已完成步骤和结果不得丢弃。
8. **统一图片规格**：本自动档产生的所有图片（QC、HVG、肘图、UMAP、Leiden、Azimuth embedding/注释/置信度/QC heatmap）全部保存为 `300 dpi`、LZW 压缩的 `.tif`，输出到 `figures_tif/`。本档不使用 PDF/PNG 作为最终图。
9. **交付**：保存处理后的 Seurat RDS、`analysis_manifest.tsv`（输入、物种、QC 阈值、ndims、resolution、Azimuth 模式和软件版本）和 `sessionInfo.txt`，并核对每个预期图层至少有一个 TIFF。

详细可执行模板见 `references/05_seurat_v6_basic_pipeline.md`。如果用户在同一句话中另外指定了参数，以用户参数覆盖本档的对应默认值。

## 工作流程

### 0. 写代码前必须确认的参数（强制）

**只要本轮任务会画出图、跑 FindMarkers/FindAllMarkers、或跑拟时序轨迹变化基因，在写任何 R 代码之前，必须先就以下参数向用户询问并拿到答复。** 不要用默认值"先跑再说"，也不要把参数选择藏在代码注释里让用户自己改。用户若明确说"你定/按常规/用默认"，才可以使用下面的推荐默认值，并在回复里注明用了哪个默认。唯一例外是上方“固定触发：单细胞转录组基本分析（Seurat v6 自动档）”，该短语已经预授权固定参数。

#### 0.1 出图参数（任何 DimPlot/FeaturePlot/VlnPlot/DotPlot/DoHeatmap/dotplot/barplot/cnetplot/emapplot/heatplot/plot_cells 等）

逐条问清：

1. **展示前 N 项**：图上要展示多少个 cluster / 通路 / 基因？（如 dotplot 的 `showCategory`、DoHeatmap 的 top N 基因、DotPlot 选哪些 marker）
2. **多图分面方式**：
   - 一张图还是多张图？
   - 多张时用 `ncol`/`nrow` 排成几列几行？还是用 `patchwork` 的 `|` / `/` 拼接？
   - 分面变量是什么（按 cluster、按样本、按时间点、按通路分类…）？
3. **配色（直接给出颜色向量）**：
   - **必须直接给具体的颜色向量**，例如 `c("#E64B35","#4DBBD5","#00A087","#3C5488",...)` 或 RColorBrewer 向量 `RColorBrewer::brewer.pal(9,"Set1")`，并在询问时把候选颜色**列出来让用户挑或确认**，而不是只说"用 Set1""用 npg 调色板"让用户自己去查。
   - 连续变量（pseudotime、表达量、p 值）用哪个渐变：如 `scale_color_viridis_c(option="plasma")`、`scale_color_gradient2(low="#3C5488",mid="white",high="#E64B35")`，同样直接写明颜色。
   - 询问时可以给 2–3 套**已展开的颜色向量**作为选项（离散型/连续型各一套），让用户选一套或自定。
4. **图形设备与尺寸**：PDF/PNG、`width`、`height`、`dpi`（默认 PDF、width=8、height=6，用户不特别说明可直接用，但要在交付时点一句）。
5. **其他展示项**：是否显示标签（`label=TRUE`）、是否 repel、图例位置、标题文字、是否按特定变量分组着色（`group.by`）、点大小/透明度等。用户没提可以用合理默认，但 cluster 标签、标题、分组着色这三项在图里属于关键信息，仍要确认。

#### 0.2 FindMarkers / FindAllMarkers 差异分析参数

在写差异分析代码前，逐条问：

1. **`only.pos`**：是否只保留上调基因（`only.pos=TRUE`），还是上下调都要（`FALSE`）？默认推荐：做 marker 鉴定 / ORA 富集时用 `TRUE`；做完整差异对比、GSEA 时用 `FALSE`。
2. **`logfc.threshold`**：log2FC 的最低限？默认推荐 `0.25`（Seurat 默认）；严格筛选可用 `0.5`，宽松筛选可用 `0.1`。
3. **结果筛选阈值**：拿到差异结果后，是否按 `p_val_adj`（校正后 p 值，推荐）或 `p_val` 过滤？阈值是多少（默认 `0.05`）？是否同时按 `avg_log2FC` 再做一次过滤？
4. （可选）`min.pct`（默认 0.25）、`test.use`（默认 `"wilcox"`，大样本可考虑 `"MAST"`/`"DESeq2"`）这两个用户通常不关心，但如果检测基因数或细胞数差异很大要主动提一句。

询问时把推荐默认值写在括号里，用户回复"默认"即可整组采用。

#### 0.3 拟时序轨迹变化基因参数（graph_test / fit_models / BEAM）

在写轨迹差异基因代码前，逐条问：

1. **筛选依据组合**：用哪些指标筛选轨迹变化基因？Monocle3 `graph_test` 返回 `morans_I`、`p_value`、`q_value`、`status`。常见组合：
   - 宽松：`q_value < 0.05`
   - 标准：`q_value < 0.05 & morans_I > 0.1`
   - 严格：`q_value < 0.01 & morans_I > 0.25`
2. **取多少个基因进入下游**：Top N（如 50/100/200）用于热图、模块聚类、可视化？
3. **若用 `fit_models` 回归**：模型公式里除了 `pseudotime` 还要不要校正 `batch`/样本/细胞周期？系数表筛选阈值（`q_value`、`term`）怎么定？
4. **若用 Monocle2 BEAM**：`branch_point` 选第几个分叉？`qval` 阈值？

#### 0.4 询问方式

- 把上面 0.1–0.3 中**本轮真正涉及**的参数整合成一段简短清单问用户，不要一次把所有类别都抛出去（没画图就别问配色，没跑轨迹就别问 Moran's I）。
- 用户已经在需求里明确给过的参数（如"resolution=0.3""top20"）不要重复问。
- 询问时给出推荐默认值和一句话理由，方便用户直接说"按默认"。
- 拿到答复后再写代码；代码里的参数必须和用户答复一致，不要在代码里偷偷改。

### 1. 先明确分析目标与输入
确认：物种（人/小鼠，决定注释数据库 `org.Hs.eg.db` / `org.Mm.eg.db`）、数据来源（10x CellRanger 输出的 filtered_feature_bc_matrix、h5 还是现成 Seurat/rds）、分析目的（标准聚类？多组整合？轨迹？通讯？）。不要在用户没给数据路径时编造路径，用占位符并提示替换。

**项目目录约定（强制）**：
- 在用户指定的项目文件夹（或原始数据所在文件夹）下新建 `figures/` 子目录（`dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)`）。
- **脚本运行中产生的所有图片**（QC 小提琴图/散点图、UMAP/tSNE、marker 热图/点图、富集 dotplot/cnetplot、拟时序图、CellChat 图等）**全部保存到 `figures/` 目录**，统一用 `ggsave(file.path(fig_dir, "xxx.pdf"), ...)`，不要散落在工作目录。
- 固定触发的“基本分析自动档”改用 `figures_tif/`，所有图片均为 `300 dpi` `.tif`；不要在该档同时输出 PDF/PNG 终稿。
- 数据文件（RData/RDS/CSV）保存在项目根目录，和 `figures/` 平级；不要把图片和数据混在一起。
- 脚本开头定义 `proj_dir` 和 `fig_dir <- file.path(proj_dir, "figures")`，后续所有保存路径都基于这两个变量拼接。

### 2. 标准 Seurat v6 流程（向后兼容 v5）
参考模板 `references/01_seurat_qc_clustering.md`。关键步骤顺序固定：

1. `Read10X()` / `Read10X_h5()` 读入 → `CreateSeuratObject(min.cells=3, min.features=200)`

2. **QC 前必须先出图、再问阈值（强制）**：
   - 计算 `percent.mt`（人 `^MT-`，小鼠 `^mt-`）和 `percent.rb`（人 `^RP[SL]`，小鼠 `^Rp[sl]`）。
   - **必须把 `nFeature_RNA`、`nCount_RNA`、`percent.mt`、`percent.rb` 四个指标的分布图保存到 `figures/` 目录**（`VlnPlot` 四联小提琴图 + `nCount_RNA` vs `percent.mt`、`nCount_RNA` vs `nFeature_RNA` 两张散点图），并把图片路径呈现给用户。
   - 然后**询问用户质控阈值**：`nFeature` 上下限、`nCount` 上下限、`percent.mt` 上限、是否过滤 `percent.rb` 异常值。根据图中分布给出建议值（如看到小提琴长尾在 6000 附近可建议 `<6000`），但最终阈值由用户确认。
   - 用户确认后才执行 `subset()` 过滤。禁止跳过出图直接套默认阈值。

3. **双细胞筛查（询问用户是否需要，不默认开启）**：
   - 询问用户是否使用 [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) 进行双细胞（doublet）预测。说明：DoubletFinder 适用于单个样本（或同一样本拆 multiple lanes），**不要在多样本 merge/整合后的数据上直接跑**（会产生不可能的跨样本人造 doublet）。
   - 用户同意后，流程：
     a. 先用标准流程（NormalizeData + FindVariableFeatures + ScaleData + RunPCA）预处理，**注意：DoubletFinder 的 `paramSweep` 有 `sct` 参数，用 SCTransform 时设 `sct=TRUE`，用 LogNormalize 时设 `sct=FALSE`**。
     b. `sweep.res <- paramSweep(seu, PCs = 1:dims, sct = FALSE)` → `sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)` → `bcmvn <- find.pK(sweep.stats)`；画 BCmvn 图帮用户选 pK（取 BCmvn 最大值对应 pK）。
     c. 估计双细胞数：`nExp_poi <- round(doublet_rate * ncol(seu))`，doublet_rate 根据 10x 上样量查官方表格（典型 0.4%–8%，教程默认 7.5%），需向用户确认。
     d. 若已有聚类注释，用 `modelHomotypic(annotations)` 估计同型 doublet 比例并校正：`nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))`。
     e. 跑两次 `doubletFinder(seu, PCs=1:dims, pN=0.25, pK=..., nExp=nExp_poi, sct=FALSE)` 和带 `reuse.pANN` 的校正版本；最终用校正后的 DF.class 列过滤 `subset(seu, DF.class_xxx == "Singlet")`。
   - 把 pK-BCmvn 图、doublet UMAP 着色图保存到 `figures/`。
   - DoubletFinder 2023 年后已兼容 Seurat v5，函数名去掉 `_v3` 后缀（直接 `paramSweep`/`doubletFinder`）；安装：`remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force=TRUE)`。

4. **标准化（常规任务默认 LogNormalize；固定触发的基本分析档强制 SCTransform）**：
   - **默认使用 LogNormalize**：`NormalizeData(normalization.method="LogNormalize", scale.factor=1e4)` → `FindVariableFeatures(selection.method="vst", nfeatures=2000)` → `ScaleData`。
   - **只有用户明确要求 SCTransform，或触发“基本分析自动档”时才用**：`SCTransform(seu, vst.flavor="v2", vars.to.regress=...)`，后续 PCA/聚类在 SCT assay 上进行。
   - 不要主动推荐或默认切换到 SCTransform；脚本中把 SCT 方案写成注释掉的备选块。

5. 降维：`RunPCA` → `ElbowPlot`（保存到 figures/）选 PC 数 → `FindNeighbors(dims=1:N)` → `FindClusters(resolution=...)` → `RunUMAP`/`RunTSNE`
6. 聚类与注释：`DimPlot`（保存到 figures/）、`FindAllMarkers`、`FindMarkers`——**调用前必须按 §0.2 询问 `only.pos`、`logfc.threshold` 和 p/q 值筛选**；依据经典 marker 手动重命名 cluster（`RenameIdents`），或用 SingleR/CellTypist/AzimuthAPI 做自动注释
7. 可视化：`DimPlot(label=TRUE)`、`FeaturePlot`、`VlnPlot`、`DotPlot`、`DoHeatmap`。**画图前必须按 §0.1 确认展示数量、分面、配色（直接给颜色向量）等参数；所有图保存到 figures/。**

### 3. 多数据集/多样本整合
- Seurat v6 继续首选 `IntegrateLayers`：`IntegrateLayers(object, method=CCAIntegration|RPCAIntegration|HarmonyIntegration)`，然后在 integrated assay 上 `FindNeighbors/FindClusters/RunUMAP`。v4 的 `FindIntegrationAnchors + IntegrateData` 仍兼容，但不要和新流程混用。
- 大数据用 BPCells 磁盘后端 + sketch 抽样（见 Seurat v6 的 sketch 与 scSLIDE vignette）。
- Harmony 需安装 `harmony` 包，用 `IntegrateLayers(method=HarmonyIntegration)` 或在 `RunHarmony` 后直接用 harmony 维度做聚类。
- 整合前后都要展示 `DimPlot(group.by="orig.ident")` 评估批次混合效果。

### 4. 差异基因功能富集（clusterProfiler）
参考模板 `references/02_clusterprofiler_enrichment.md`。

- **输入 ID 类型**：单细胞 DE 结果默认是 SYMBOL（基因名）。clusterProfiler 的 `enrichGO/gseGO` 要求 ENTREZID，用 `bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)` 转换；转换会丢失一部分基因（别名/已撤回），这是正常现象，报告转换前后数量。
- **ORA（过表达分析，用 FindAllMarkers 上调基因）**：`enrichGO(gene, OrgDb, keyType="ENTREZID", ont="BP|MF|CC|ALL", pAdjustMethod="BH", qvalueCutoff=0.05, readable=TRUE)`；KEGG 用 `enrichKEGG(gene, organism="hsa"|"mmu", pvalueCutoff=0.05)`；Reactome 用 `ReactomePA::enrichPathway`。
- **GSEA（需要排好序的基因列表）**：用 `FindMarkers` 返回的 `avg_log2FC` 排序（不要用 p 值），`geneList <- sort(setNames(deg$avg_log2FC, deg$entrez), decreasing=TRUE)`，然后 `gseGO(geneList, OrgDb, ont="BP", pvalueCutoff=0.05)` 或 `gseKEGG`。
- **可视化**用配套的 `enrichplot`：`dotplot`（点图）、`barplot`、`cnetplot`（基因-术语网络）、`emapplot`（术语重叠）、`heatplot`、`gseaplot2`（GSEA  running enrichment score）、`ridgeplot`。配合 `clusterProfiler::pairwise_termsim` 做语义相似性简化。
- 多组对比用 `compareCluster(geneList~group, fun="enrichGO")` 统一画对比图。
- clusterProfiler 4.0 已集成 LLM 自动解读（需配置），默认不主动使用，除非用户要求。

### 5. GSEA + MSigDB 基因集
- MSigDB 基因集在 R 里用 `msigdbr::msigdbr(species="Homo sapiens", collection="H"|"C2"|"C5"|...)` 获取（当前 2026.1 对应 Ensembl 115，msigdbr 版本要匹配）。
- 经典用法：取 Hallmark（H）集合 → `split(x=genes$gene_symbol, f=genes$gs_name)` → `fgsea::fgsea(pathways, geneList, minSize=15, maxSize=500, eps=0)` 快速 GSEA，或 `clusterProfiler::GSEA(geneList, TERM2GENE=...)`。
- 单样本富集用 `GSVA::gsva` 或 `singscore`，或通过 GenePattern 跑 ssGSEA。
- 运行 GSEA 时强调：基因列表必须按排序统计量（log2FC 或 signed -log10 pvalue）从大到小排好且含全部检测基因，不要只输入差异基因。

### 6. Monocle3 轨迹与拟时序
参考模板 `references/03_monocle3_trajectory.md`。核心调用链（来自官方 docs）：
1. 构建 `cell_data_set`：`new_cell_data_set(expression_matrix, cell_metadata, gene_metadata)`；从 Seurat 转换用 `seurat_wrappers::as.cell_data_set` 或手工取 `GetAssayData(seurat, layer="counts")` + `seurat@meta.data`。
2. `preprocess_cds(cds, num_dim=50, method="PCA")` → `reduce_dimension(cds, reduction_method="UMAP", preprocess_method="PCA")`
3. `cluster_cells(cds, resolution=1e-4 到 1e-2)`（注意：Monocle3 的 resolution 是聚类图的分辨率，和 Seurat 不同，常用 `1e-4`/`1e-3` 这种小值）
4. `learn_graph(cds, use_partition=TRUE, learn_graph_control=list(ncenter=...))` 学习主图
5. `order_cells(cds, root_cells=...)` 或 `order_cells(cds, root_pr_nodes=...)`：必须指定起点。起点选择依据：已知的祖细胞/干/起始状态 cluster，或用 Seurat 里高表达的干性 marker；用户没给方向时主动提示需要根节点信息。
6. 可视化：`plot_cells(cds, color_cells_by="pseudotime", label_branch_points=TRUE, label_leaves=TRUE, label_roots=TRUE)`；按基因 `plot_cells(cds, genes=c(...), show_trajectory_graph=FALSE)`。
7. 差异分析：`graph_test(cds, neighbor_graph="principal_graph", cores=4)` 找沿轨迹变化的基因；**筛选前必须按 §0.3 确认 q 值、Moran's I 阈值与 Top N**。`top_markers(cds, group_cells_by="cluster")`；回归分析用 `fit_models(cds, model_formula_str="~pseudotime + batch")` + `coefficient_table`。
8. 分支分析：`choose_graph_segments(cds)` 交互式选分支，或用 `choose_cells`/`subset` 提取分支后再 `graph_test`。

**Monocle 2（备用）**：当需要 BEAM 分支依赖分析、DDRTree 二维轨迹或与旧代码兼容时使用。关键函数：`importCDS(seurat_cds)` 或 `newCellDataSet`、`estimateSizeFactors/dispersion`、`reduceDimension(method="DDRTree", max_components=2)`、`orderCells`、`BEAM(HSMM, branch_point=1, relative_expr=TRUE)`、`plot_genes_branched_heatmap/plot_genes_branched_pseudotime`。Monocle 2 要求 exprs 为相对表达值（FPKM/TPM 或 Census 计数），用相对 counts 时考虑 Census 转换。

### 7. CellChat 细胞通讯
参考模板 `references/04_cellchat_communication.md`。核心流程：
1. `createCellChat(object=data.input, meta=meta, group.by="labels")`，`data.input` 用 normalized data matrix（Seurat 中 `LayerData(seurat, layer="data")` 或 `GetAssayData(seurat, slot="data")`）。
2. `CellChatDB <- CellChatDB.human`（鼠用 `.mouse`），可 `subsetDB(search="Secreted Signaling"|"ECM-Receptor"|"Cell-Cell Contact")`，`cellchat@DB <- CellChatDB.use`。
3. `subsetData` → `identifyOverExpressedGenes` → `identifyOverExpressedInteractions` → `projectData(PPI.human)` → `computeCommunProb(raw.use=TRUE)` → `filterCommunication(min.cells=10)` → `computeCommunProbPathway` → `aggregateNet`。
4. 可视化：`netVisual_circle`、`netVisual_aggregate(signaling=..., layout="hierarchy"|"circle"|"chord")`、`netVisual_bubble(sources.use=, targets.use=, signaling=)`、`netAnalysis_contribution`、`netVisual_individual(pairLR.use=)`。
5. 结果提取：`subsetCommunication(cellchat)` 返回数据框，可按 sources/targets/signaling 筛选并 `write.csv`。

### 8. Seurat v6 增量：AzimuthAPI 自动注释

参考 `references/06_azimuth_scslide_v6.md`。V6 vignette 展示的 Pan-human Azimuth 入口是 `AzimuthAPI`，核心输出是层级神经网络标签和 128 维 embedding；旧式 `RunAzimuth()` 组织参考映射在确有匹配参考时仍可用，但不要把两条工作流或输出列名混为一谈。

1. 安装并载入：`remotes::install_github("satijalab/AzimuthAPI")`，`library(AzimuthAPI)`。输入必须是已完成 `NormalizeData()` 的人源 Seurat 对象。
2. 云端：`obj <- CloudAzimuth(obj)`。这是最简单的入口，但会走 AWS API；有数据外送限制时不要使用云端。
3. 本地：创建 `python=3.9` 的 conda 环境，`pip install panhumanpy`，用 `reticulate::use_condaenv("pan-human-azimuth", required=TRUE)` 指定解释器，再运行 `obj <- ANNotate(obj)`。首次运行会下载并缓存模型权重。
4. 输出字段：`final_level_labels`、`final_level_confidence`、`full_hierarchical_labels`、`full_consistent_hierarchy`、`azimuth_broad`、`azimuth_medium`、`azimuth_fine`，以及 reduction `azimuth_embed`。层级标签用 `|` 分隔。
5. 可视化：用 `azimuth_medium` 看主要细胞类型，用 `PrepLabel(obj, "azimuth_fine", "azimuth_fine_filtered", cutoff=20)` 合并小标签，用 `final_level_confidence` 评估可靠性；需要时用 `RunUMAP(reduction="azimuth_embed", dims=1:128, reduction.name="azimuth_umap")`。
6. 低置信度细胞不要静默删除。可增加 `azimuth_confidence_pass <- final_level_confidence > 0.5`，把低置信度细胞仅显示为 `Unassigned_low_confidence`，并保留原始预测列。
7. `make_azimuth_QC_heatmaps()` 返回按 broad label 分组的图列表，逐图保存，不能只打印到设备。所有这些图也要按用户指定的设备、尺寸和 dpi 保存。

**物种边界**：`AzimuthAPI` 是 Pan-human 模型，只适用于人源数据。小鼠、其他物种或未确定物种的数据不能套用；应使用匹配物种的人工 marker、SingleR/celldex 参考或用户指定的其他参考模型。

### 9. Seurat v6 增量：scSLIDE 样本水平分析

参考 `references/06_azimuth_scslide_v6.md`。`scSLIDE` 是 **Single-cell Sample-Level Inference using Density Estimation**，用于多样本数据，不是空间转录组中的 Slide-seq。它把每个样本表示为“landmark 细胞周围的密度向量”，再在样本水平做聚类、轨迹和差异分析。

输入至少需要：一个多样本 Seurat 对象、样本/供体列、用于监督 PLS 的表型或多变量响应列、细胞类型列或可用于构造细胞类型嵌入的 reduction。只有单样本、没有重复样本或没有可用响应变量时，不适用 scSLIDE。

推荐先构造两类 embedding：
- 细胞类型 embedding：Pan-human Azimuth 的 `azimuth_embed`，或 PCA/其他整合与 reference mapping reduction。
- 表型 embedding：`RunPLS()` 对响应变量做监督降维。

最简封装流程：
```r
object <- PrepareSampleObject(
  object = object,
  assay = "RNA",
  add.hvg = TRUE,
  group.by.CorTest = "azimuth_label",
  Y = "Disease",
  sketch.training = TRUE,
  group.by.Sketch = "Donor",
  ncells.per.group = 1000,
  ncells.landmark = 2000,
  ncomp = 20,
  pls.function = "cppls",
  k.nn = 5,
  name.reduction.1 = "azimuth_embed",
  dims.reduction.1 = 1:128,
  fix.wnn.weights = c(0.5, 0.5)
)

sample_obj <- GenerateSampleObject(
  object = object,
  sketch.assay = "LANDMARK",
  nn.name = "weighted.nn",
  group.by = "Donor",
  rename.group.by = "azimuth_label",
  k.nn = 5,
  normalization.method = "ChiSquared",
  add.meta.data = TRUE,
  return.seurat = TRUE
)
```

需要逐步控制时按以下链条：
1. `SketchDataByGroup()` 按样本均匀抽样训练集；样本细胞数已平衡时可省略。
2. `SketchData(method="LeverageScore")` 选覆盖异质性的 landmark 细胞。
3. `QuickCorTest()` 可补选与响应变量相关的基因；HVG 已足够时省略。
4. `RunPLS(pls.function="plsr"|"cppls"|"spls")` 得到表型 embedding；`ProjectCellEmbeddings()` 投影到全部细胞。
5. `FindmmNN()` 合并 cell-type 与 phenotype 两类 embedding；只用单一模态时用 `FindsmNN()`。
6. `GenerateSampleObject(normalization.method="ChiSquared")` 生成 landmark-by-sample 密度矩阵并返回样本水平 Seurat 对象。
7. `RunDiffusionMap()` 对样本水平矩阵降维；可用 `princurve::principal_curve()` 把前几个 DC 拟合为联合伪轨迹。
8. 轨迹相关差异表达使用 pseudobulk：`AggregateExpression()` → 合并样本轨迹 metadata → 按细胞类型 `subset()` → `TrajDETest(traj.var=..., latent.vars=...)`；热图用 `PlotAnnotatedHeatmap()`。

大规模数据优先 BPCells on-disk assay。V6 的函数 `setThreads(n)`/`getThreads()` 可为 `NormalizeData`、`RunPCA`、`FindNeighbors`、`FindClusters`、`RunUMAP` 及部分整合函数启用多线程；启用后随机数和浮点误差可能使 UMAP/簇边界轻微变化，应固定 seed 并记录线程数。

scSLIDE 常见边界：
- 不要把 landmark 数设得远少于细胞类型数，否则稀有群体会在密度矩阵中消失。
- `group.by` 必须真的是生物学样本/供体，不能使用技术 barcode lane 冒充重复样本。
- 表型 `Y` 不可缺失；多分类或连续变量应先明确编码，并检查每个水平是否有足够样本。
- 样本水平结论不能由单样本推断；样本数太少时优先描述性展示，不强行做组间显著性检验。
- 患者轨迹与细胞轨迹是两种不同分析；不要拿 Monocle3 的细胞 pseudotime 替代 scSLIDE 的 sample trajectory。

### 10. 输出与交付
- 每次输出一个完整可运行的 R 脚本（`.R` 文件），分段加注释；如用户要求 Notebook 风格则输出 `.Rmd` 或 Quarto `.qmd`。
- 脚本开头给出 `sessionInfo` 占位与 `install.packages/BiocManager::install` 安装代码块（注释掉，让用户按需执行）。
- 关键参数（PC 数、resolution、过滤阈值、qvalueCutoff、根节点选择）必须说明为什么这么选，而不是直接扔参数。
- 图表用 `ggplot2`/Seurat 自带；常规任务可用 `ggsave("xxx.pdf", width=8, height=6)`；基本分析自动档必须使用 `device="tiff", dpi=300, compression="lzw"`。
- 引用工具时标注官方文档链接，在结论附近以注释或正文形式给出。

## 边界与常见坑

- **写代码前先问参数（最高优先级）**：凡是要画图、跑 `FindMarkers`/`FindAllMarkers`、跑拟时序 `graph_test`/`fit_models`/`BEAM`，必须先按 §0.1/§0.2/§0.3 与用户确认参数；禁止"先跑后改"。配色必须直接给颜色向量，不要只报调色板名字让用户自己查。
- **固定触发例外**：“帮我完成一个单细胞转录组基本分析”已经预授权上方自动档参数，不再重复询问；只有输入不唯一、物种无法可靠判断或数据导致异常 QC 损失时才停下确认。
- **Seurat 对象版本**：v5/v6 用 `LayerData` 取 counts/data；旧代码的 `GetAssayData(object, slot="counts")` 仍可用但会提示替代。整合函数从 v4 的 `IntegrateData` 改为 `IntegrateLayers`，不要混用导致 assay 结构错误。
- **Azimuth 物种限制**：`AzimuthAPI` 的 Pan-human 模型只能用于人；不得把 `CloudAzimuth()`/`ANNotate()` 产生的标签套到小鼠或其他物种。
- **scSLIDE 不是空间分析**：它分析多样本细胞密度和样本轨迹，不是 Slide-seq 空间坐标；不要用 `Load10X_Spatial()`/`SpatialPlot()` 实现 scSLIDE。
- **SCTransform 参数名**：用户说 `vars_to_regress` 时，R 代码应写 `vars.to.regress`；回归变量名必须实际存在于 Seurat metadata。
- **Leiden 参数**：`FindClusters(algorithm=4)` 才是 Leiden；自动档固定 `resolution=0.5`，而 `algorithm=1` 是 Louvain。
- **物种前缀**：人线粒体 `^MT-`、核糖体 `^RP[SL]`；小鼠 `^mt-`、`^Rp[sl]`；KEGG organism `hsa` vs `mmu`；org.db 不要装错。
- **ID 转换丢失**：`bitr` 会因别名/废弃 ID 丢基因，属于正常；富集结果里基因数少于输入时不要报错，报告转换率。
- **GSEA 排序**：必须包含所有检测到的基因、按统计量连续排序且无 NA；只传差异基因会被 GSEA 拒绝或结果失真。
- **Monocle3 根节点**：不指定 root 就 `order_cells` 会报错或得到无意义 pseudotime；务必先确认起点 cluster。
- **Monocle3 resolution 单位**：和 Seurat 不同，常用 `1e-4`、`5e-4`、`1e-3`，不要用 0.5。
- **大数据内存**：>10 万细胞提示用 BPCells/sketch；不要让用户一次性把数百万细胞读进 dense matrix。
- **多线程可复现性**：`setThreads()` 可能改变随机数顺序和浮点结果；记录线程数和 seed，并在比较结果时保持配置一致。
- **不要编造函数或参数**：不确定的包版本/函数签名以官方文档为准；clusterProfiler、Seurat、Monocle3 都在快速迭代，优先用当前稳定版（Seurat v6、clusterProfiler 4.x、Monocle3 beta）。
- 本技能只产出 R 代码与分析方案，不在沙箱里实际运行 R（运行环境无 R/Bioconductor）；用户拿到脚本在本地/RStudio 运行。

## 参考文件

- `references/01_seurat_qc_clustering.md` — Seurat v5/v6 兼容标准流程模板
- `references/02_clusterprofiler_enrichment.md` — clusterProfiler ORA + GSEA 模板
- `references/03_monocle3_trajectory.md` — Monocle3 拟时序模板
- `references/04_cellchat_communication.md` — CellChat 细胞通讯模板
- `references/05_seurat_v6_basic_pipeline.md` — 基本分析自动档：SCTransform/肘图/UMAP/Leiden 0.5/Azimuth/300 dpi TIFF
- `references/06_azimuth_scslide_v6.md` — Seurat v6 AzimuthAPI 与 scSLIDE 参考流程

代码量较大或用户要求完整 pipeline 时，把对应模板读出来，在其基础上根据用户数据路径、物种、分组变量做替换与增删，而不是从零生成。
