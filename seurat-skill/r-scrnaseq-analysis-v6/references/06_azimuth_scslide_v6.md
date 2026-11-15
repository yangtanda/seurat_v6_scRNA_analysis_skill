# =============================================================================
# 06_azimuth_scslide_v6.R
# Seurat v6 additions:
#   1. AzimuthAPI Pan-human cell-type annotation
#   2. scSLIDE sample-level analysis
#
# Official vignettes:
#   https://satijalab.org/seurat/articles/intro_azimuthapi_vignette
#   https://satijalab.org/seurat/articles/scSLIDE_vignette
#
# Run the Azimuth section before scSLIDE when an Azimuth cell-type embedding is
# available. scSLIDE requires multiple biological samples and a phenotype or
# other response variable; it is not a spatial Slide-seq workflow.
# =============================================================================

# ---- 0. Packages ------------------------------------------------------------
# install.packages(c(
#   "Seurat", "princurve", "ggplot2", "ggridges", "BPCells",
#   "reticulate", "patchwork", "igraph"
# ))
# remotes::install_github("satijalab/AzimuthAPI")
# BiocManager::install("destiny")

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

# Use the project-specific path and a 300 dpi TIFF output directory when the
# user requests TIFF figures. Otherwise follow the active request.
proj_dir <- "path/to/project"
fig_dir <- file.path(proj_dir, "figures_tif")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_tif <- function(plot, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path(fig_dir, filename),
    plot = plot,
    device = "tiff",
    dpi = 300,
    width = width,
    height = height,
    units = "in",
    compression = "lzw",
    bg = "white"
  )
}

# =============================================================================
# Part A. AzimuthAPI: Pan-human hierarchical annotation
# =============================================================================

# AzimuthAPI is Pan-human only. Stop rather than applying a human model to
# mouse or another species.
species <- "human"
stopifnot(identical(species, "human"))

object <- readRDS(file.path(proj_dir, "seurat_processed.rds"))

# AzimuthAPI expects normalized RNA data.
object <- NormalizeData(object, assay = "RNA", verbose = FALSE)
original_assay <- DefaultAssay(object)
DefaultAssay(object) <- "RNA"

if (!requireNamespace("AzimuthAPI", quietly = TRUE)) {
  stop("Install AzimuthAPI with remotes::install_github('satijalab/AzimuthAPI').")
}

# Option 1: cloud API. This sends the analysis request to AWS.
object_ann <- AzimuthAPI::CloudAzimuth(object)
DefaultAssay(object_ann) <- if (original_assay %in% Assays(object_ann)) {
  original_assay
} else {
  "RNA"
}

# Option 2: local panhumanpy environment.
# conda create -n pan-human-azimuth python=3.9
# conda activate pan-human-azimuth
# pip install panhumanpy
#
# reticulate::use_condaenv("pan-human-azimuth", required = TRUE)
# object_ann <- AzimuthAPI::ANNotate(object)

# Key output columns:
#   final_level_labels, final_level_confidence, full_hierarchical_labels,
#   full_consistent_hierarchy, azimuth_broad, azimuth_medium, azimuth_fine
azimuth_metadata <- c(
  "final_level_labels",
  "final_level_confidence",
  "full_hierarchical_labels",
  "full_consistent_hierarchy",
  "azimuth_broad",
  "azimuth_medium",
  "azimuth_fine"
)
stopifnot(all(azimuth_metadata %in% colnames(object_ann[[]])))
stopifnot("azimuth_embed" %in% Reductions(object_ann))

# Preserve confidence information instead of silently deleting low-confidence
# cells. Filtering can be requested separately after inspecting the histogram.
object_ann$azimuth_confidence_pass <- object_ann$final_level_confidence > 0.5

p_confidence_hist <- ggplot(
  object_ann[[]],
  aes(x = final_level_confidence)
) +
  geom_histogram(bins = 20, fill = "#4DBBD5", color = "white") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "#E64B35") +
  labs(
    x = "Azimuth final-level confidence",
    y = "Cell count",
    title = "Pan-human Azimuth confidence"
  ) +
  theme_bw()
save_tif(p_confidence_hist, "Azimuth_confidence_histogram.tif", width = 7, height = 5)

# Merge sparse labels for display while retaining the original predictions.
if ("PrepLabel" %in% getNamespaceExports("AzimuthAPI")) {
  object_ann <- AzimuthAPI::PrepLabel(
    object_ann,
    "azimuth_fine",
    "azimuth_fine_filtered",
    cutoff = 20
  )
}

azimuth_dims <- min(128L, ncol(Embeddings(object_ann[["azimuth_embed"]])))
object_ann <- RunUMAP(
  object_ann,
  reduction = "azimuth_embed",
  dims = 1:azimuth_dims,
  reduction.name = "azimuth_umap",
  seed.use = 42,
  verbose = FALSE
)

p_medium <- DimPlot(
  object_ann,
  reduction = "azimuth_umap",
  group.by = "azimuth_medium",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Azimuth medium labels")
save_tif(p_medium, "Azimuth_UMAP_medium.tif", width = 9, height = 6)

fine_column <- if ("azimuth_fine_filtered" %in% colnames(object_ann[[]])) {
  "azimuth_fine_filtered"
} else {
  "azimuth_fine"
}
p_fine <- DimPlot(
  object_ann,
  reduction = "azimuth_umap",
  group.by = fine_column,
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Azimuth fine labels")
save_tif(p_fine, "Azimuth_UMAP_fine.tif", width = 10, height = 6)

p_confidence <- FeaturePlot(
  object_ann,
  features = "final_level_confidence",
  reduction = "azimuth_umap",
  min.cutoff = 0,
  max.cutoff = 1
) +
  ggtitle("Azimuth confidence on azimuth_embed UMAP")
save_tif(p_confidence, "Azimuth_UMAP_confidence.tif", width = 8, height = 6)

if ("make_azimuth_QC_heatmaps" %in% getNamespaceExports("AzimuthAPI")) {
  az_heatmaps <- AzimuthAPI::make_azimuth_QC_heatmaps(object_ann)
  for (i in seq_along(az_heatmaps)) {
    panel_name <- names(az_heatmaps)[[i]]
    if (is.null(panel_name) || !nzchar(panel_name)) {
      panel_name <- paste0("panel_", i)
    }
    panel_name <- gsub("[^A-Za-z0-9_-]+", "_", panel_name)
    save_tif(
      az_heatmaps[[i]],
      paste0("Azimuth_QC_heatmap_", panel_name, ".tif"),
      width = 10,
      height = 8
    )
  }
}

saveRDS(object_ann, file.path(proj_dir, "seurat_azimuth_annotated.rds"))

# =============================================================================
# Part B. scSLIDE: sample-level density, clustering, and trajectory
# =============================================================================

# Required metadata:
#   Donor  - biological sample/donor identifier
#   Disease - response used by supervised PLS; replace with the actual phenotype
#   azimuth_label (or another cell-type column)
#
# Do not substitute a technical lane or one sample per group for biological
# replication. Sample-level inference requires multiple biological samples.

object <- object_ann
variable_metadata <- c("Donor", "Disease", "azimuth_label")
missing_metadata <- setdiff(variable_metadata, colnames(object[[]]))
if (length(missing_metadata) > 0L) {
  stop(
    "Missing scSLIDE metadata: ",
    paste(missing_metadata, collapse = ", ")
  )
}
object$Disease <- as.factor(object$Disease)

object <- NormalizeData(object, assay = "RNA", verbose = FALSE)
object <- FindVariableFeatures(object, nfeatures = 2000, verbose = FALSE)

# Seurat v6 supports multithreading for many core functions. Record the setting
# because multithreading can alter random-number order and floating-point results.
setThreads(8)

# Wrapper workflow from the official vignette.
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
  dims.reduction.1 = 1:azimuth_dims,
  fix.wnn.weights = c(0.5, 0.5),
  verbose = TRUE
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

# Center, but do not z-score, the sample-level density matrix before diffusion
# map analysis.
sample_obj <- ScaleData(
  sample_obj,
  features = rownames(sample_obj),
  do.scale = FALSE,
  do.center = TRUE,
  verbose = FALSE
)

sample_assay <- DefaultAssay(sample_obj)
sample_obj <- RunDiffusionMap(
  sample_obj,
  assay = sample_assay,
  layer = "scale.data",
  features = rownames(sample_obj),
  ncomp = 10,
  verbose = FALSE
)

p_samples <- DimPlot(
  sample_obj,
  reduction = "DiffMap",
  group.by = "Disease"
) +
  ggtitle("scSLIDE sample-level diffusion map")
save_tif(p_samples, "scSLIDE_sample_diffusion_map.tif", width = 8, height = 6)

# Store leading diffusion components as sample metadata for trajectory models.
sample_obj$traj_DC1 <- sample_obj[["DiffMap"]]@cell.embeddings[, 1]
sample_obj$traj_DC2 <- sample_obj[["DiffMap"]]@cell.embeddings[, 2]
if (ncol(sample_obj[["DiffMap"]]@cell.embeddings) >= 3L) {
  sample_obj$traj_DC3 <- sample_obj[["DiffMap"]]@cell.embeddings[, 3]
}

# Optional principal curve over the first two DCs.
princurve_res <- princurve::principal_curve(
  sample_obj[["DiffMap"]]@cell.embeddings[, 1:2],
  plot = FALSE,
  maxit = 20,
  trace = TRUE
)
sample_obj$traj_overall <- princurve_res$lambda

p_curve <- DimPlot(
  sample_obj,
  reduction = "DiffMap",
  group.by = "Disease"
) +
  geom_path(
    data = data.frame(
      x = princurve_res$s[princurve_res$ord, 1],
      y = princurve_res$s[princurve_res$ord, 2]
    ),
    aes(x = x, y = y),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  ggtitle("scSLIDE principal curve")
save_tif(p_curve, "scSLIDE_principal_curve.tif", width = 8, height = 6)

# ---- Optional trajectory-associated differential expression ----------------
# Replace "CD14 monocyte" with the cell type being tested. The sample-level
# trajectory metadata must be joined by Donor before TrajDETest().
bulk_obj <- AggregateExpression(
  object,
  assays = "RNA",
  group.by = c("Donor", "azimuth_label"),
  return.seurat = TRUE
)

# Attach the biological sample labels back to the pseudobulk columns. The
# aggregate column names follow the official vignette's Donor_celltype format.
bulk_meta <- object@meta.data[, c("Donor", "azimuth_label"), drop = FALSE]
bulk_meta <- bulk_meta[!duplicated(bulk_meta), , drop = FALSE]
rownames(bulk_meta) <- paste0(
  gsub("_", "-", bulk_meta$Donor),
  "_",
  gsub("_", "-", bulk_meta$azimuth_label)
)
bulk_obj <- AddMetaData(bulk_obj, metadata = bulk_meta)

bulk_counts <- LayerData(bulk_obj, assay = "RNA", layer = "counts")
bulk_obj$nCount_RNA <- Matrix::colSums(bulk_counts)
bulk_obj$nFeature_RNA <- Matrix::colSums(bulk_counts > 0)

target_cell_type <- "CD14 monocyte"
sub_bulk_obj <- subset(bulk_obj, subset = azimuth_label == target_cell_type)
if (ncol(sub_bulk_obj) == 0L) {
  stop("target_cell_type was not found in the pseudobulk object.")
}

# The pseudobulk column names encode Donor_celltype; the scSLIDE vignette first
# attaches the matching Donor metadata and then merges sample trajectory values.
bulk_meta <- data.frame(
  Donor = sub_bulk_obj$Donor,
  row.names = colnames(sub_bulk_obj)
)
sample_trajectory <- sample_obj[[]][, c("Donor", "traj_DC1")]
if ("traj_DC2" %in% colnames(sample_obj[[]])) {
  sample_trajectory$traj_DC2 <- sample_obj$traj_DC2
}
if ("traj_DC3" %in% colnames(sample_obj[[]])) {
  sample_trajectory$traj_DC3 <- sample_obj$traj_DC3
}
bulk_meta <- merge(
  bulk_meta,
  sample_trajectory,
  by = "Donor",
  all.x = TRUE,
  sort = FALSE
)
rownames(bulk_meta) <- colnames(sub_bulk_obj)
sub_bulk_obj <- AddMetaData(sub_bulk_obj, metadata = bulk_meta)
sub_bulk_obj$log_count <- log10(sub_bulk_obj$nCount_RNA)

de_dc1 <- TrajDETest(
  sub_bulk_obj,
  assay = "RNA",
  layer = "counts",
  traj.var = "traj_DC1",
  latent.vars = "log_count"
)
write.csv(
  de_dc1,
  file.path(proj_dir, paste0("scSLIDE_", target_cell_type, "_trajDC1_DE.csv")),
  row.names = FALSE
)

# For a joint trajectory, test traj_overall instead. Always verify that the
# trajectory column exists and that enough samples have non-missing values.

# Large-data branch:
# - Store counts with BPCells and keep the assay on disk.
# - Use the same PrepareSampleObject()/GenerateSampleObject() chain.
# - Filter samples with too few cells only when that exclusion is pre-specified.
# - Use name.reduction.1 = "azimuth_embed" after Pan-human Azimuth or a PCA
#   reduction when no cell-type-oriented embedding is available.

setThreads(1)
capture.output(sessionInfo(), file = file.path(proj_dir, "sessionInfo_v6_reference.txt"))
