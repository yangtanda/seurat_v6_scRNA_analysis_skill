# =============================================================================
# 05_seurat_v6_basic_pipeline.R
# Trigger:
#   "帮我完成一个单细胞转录组基本分析"
#
# Fixed workflow:
#   Seurat input -> species detection -> QC -> SCTransform -> PCA/elbow ->
#   UMAP -> Leiden (resolution = 0.5) -> AzimuthAPI -> 300 dpi TIFF
#
# This profile is pre-authorized by the trigger phrase. See SKILL.md before
# changing normalization, regression variables, dimensionality, or resolution.
# =============================================================================

# ---- 0. Dependencies --------------------------------------------------------
# install.packages(c(
#   "Seurat", "SeuratObject", "ggplot2", "dplyr", "patchwork",
#   "Matrix", "hdf5r", "igraph", "reticulate"
# ))
# remotes::install_github("satijalab/AzimuthAPI")
#
# AzimuthAPI has two modes:
# - CloudAzimuth(): simpler, but sends data to a cloud API.
# - ANNotate(): local Python environment with panhumanpy; use for restricted data.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

set.seed(42)
options(future.globals.maxSize = 8 * 1024^3)

# ---- 1. Project and output directories --------------------------------------
# Replace this before running. The agent should fill it from the user's request
# or the unambiguous project directory.
proj_dir <- "path/to/project"
stopifnot(dir.exists(proj_dir))
proj_dir <- normalizePath(proj_dir, winslash = "/", mustWork = TRUE)

fig_dir <- file.path(proj_dir, "figures_tif")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_tif <- function(plot, filename, width = 8, height = 6) {
  stopifnot(tools::file_ext(filename) == "tif")
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

# ---- 2. Locate input and read with Seurat -----------------------------------
# Set input_path explicitly when known. Otherwise this function searches for
# exactly one plausible 10x input. Multiple candidates stop the script instead
# of silently choosing the wrong sample.
input_path <- Sys.getenv("SCRNA_INPUT", unset = "")

find_10x_candidates <- function(root) {
  h5 <- list.files(
    root,
    pattern = "\\.h5$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  h5 <- h5[grepl("filtered|feature_bc_matrix", basename(h5), ignore.case = TRUE)]

  dirs <- list.dirs(root, recursive = TRUE, full.names = TRUE)
  is_10x_dir <- vapply(
    dirs,
    FUN.VALUE = logical(1),
    FUN = function(x) {
      mtx <- list.files(x, pattern = "^matrix\\.mtx(\\.gz)?$", full.names = FALSE)
      features <- list.files(x, pattern = "^(features|genes)\\.tsv(\\.gz)?$")
      barcodes <- list.files(x, pattern = "^barcodes\\.tsv(\\.gz)?$")
      length(mtx) == 1L && length(features) == 1L && length(barcodes) == 1L
    }
  )

  unique(c(h5, dirs[is_10x_dir]))
}

if (!nzchar(input_path)) {
  candidates <- find_10x_candidates(proj_dir)
  if (length(candidates) == 0L) {
    stop("No 10x matrix or filtered_feature_bc_matrix.h5 was found under proj_dir.")
  }
  if (length(candidates) > 1L) {
    stop(
      "Multiple candidate inputs were found; ask the user to choose one:\n",
      paste(candidates, collapse = "\n")
    )
  }
  input_path <- candidates[[1]]
}

input_path <- normalizePath(input_path, winslash = "/", mustWork = TRUE)

read_seurat_input <- function(path) {
  if (dir.exists(path)) {
    counts <- Read10X(data.dir = path)
  } else if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (!inherits(obj, "Seurat")) {
      stop("The supplied RDS is not a Seurat object.")
    }
    return(obj)
  } else if (grepl("\\.h5$", path, ignore.case = TRUE)) {
    counts <- Read10X_h5(filename = path)
  } else {
    stop("input_path must be a 10x directory, a 10x .h5 file, or a Seurat .rds.")
  }

  if (is.list(counts)) {
    preferred <- grep("Gene Expression", names(counts), value = TRUE)
    if (length(preferred) == 0L) {
      stop("Read10X returned multiple modalities; specify the gene-expression matrix.")
    }
    counts <- counts[[preferred[[1]]]]
  }

  CreateSeuratObject(
    counts = counts,
    project = "scRNA_basic",
    min.cells = 3,
    min.features = 200
  )
}

seurat <- read_seurat_input(input_path)

# ---- 3. Detect species from feature naming -----------------------------------
detect_species <- function(features) {
  features <- unique(as.character(features))

  human_house <- c("GAPDH", "ACTB", "MALAT1", "B2M", "TMSB10", "UBC")
  mouse_house <- c("Gapdh", "Actb", "Malat1", "B2m", "Tmsb10", "Ubc")

  human_score <- 2 * sum(grepl("^MT-", features)) +
    sum(grepl("^RP[SL]", features)) +
    2 * sum(features %in% human_house)

  mouse_score <- 2 * sum(grepl("^mt-", features)) +
    sum(grepl("^Rp[sl]", features)) +
    2 * sum(features %in% mouse_house)

  if (max(human_score, mouse_score) < 5L) {
    stop("Species could not be determined reliably from feature names.")
  }

  ratio <- (human_score + 1) / (mouse_score + 1)
  if (ratio >= 2) {
    return(list(
      species = "human",
      mt_pattern = "^MT-",
      rb_pattern = "^RP[SL]",
      human_score = human_score,
      mouse_score = mouse_score
    ))
  }
  if (ratio <= 0.5) {
    return(list(
      species = "mouse",
      mt_pattern = "^mt-",
      rb_pattern = "^Rp[sl]",
      human_score = human_score,
      mouse_score = mouse_score
    ))
  }

  stop(
    "Species signal is ambiguous (human score = ", human_score,
    ", mouse score = ", mouse_score, "). Ask the user to confirm species."
  )
}

species_info <- detect_species(rownames(seurat))
species <- species_info$species
seurat$species <- species

# ---- 4. Automatic QC with diagnostic TIFFs ----------------------------------
require_features <- c("nFeature_RNA", "nCount_RNA")
if (!all(require_features %in% colnames(seurat[[]]))) {
  stop("nFeature_RNA or nCount_RNA is missing from the Seurat metadata.")
}

seurat[["percent.mt"]] <- PercentageFeatureSet(
  seurat,
  pattern = species_info$mt_pattern
)
seurat[["percent.rb"]] <- PercentageFeatureSet(
  seurat,
  pattern = species_info$rb_pattern
)

qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb")
p_qc_before <- VlnPlot(
  seurat,
  features = qc_features,
  ncol = 4,
  pt.size = 0
)
save_tif(p_qc_before, "QC_violin_before_filtering.tif", width = 12, height = 5)

p_qc_mt <- FeatureScatter(
  seurat,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt",
  pt.size = 0.3
)
p_qc_feature <- FeatureScatter(
  seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  pt.size = 0.3
)
save_tif(p_qc_mt | p_qc_feature, "QC_scatter_before_filtering.tif", width = 11, height = 5)

safe_mad <- function(x) {
  value <- stats::mad(x, na.rm = TRUE)
  if (!is.finite(value) || value == 0) {
    value <- stats::IQR(x, na.rm = TRUE) / 1.349
  }
  if (!is.finite(value) || value == 0) 0 else value
}

mt_cut <- min(
  30,
  max(10, stats::median(seurat$percent.mt, na.rm = TRUE) +
    3 * safe_mad(seurat$percent.mt))
)
rb_cut <- max(
  5,
  stats::median(seurat$percent.rb, na.rm = TRUE) +
    3 * safe_mad(seurat$percent.rb)
)

keep_cells <- seurat$nFeature_RNA >= 200 &
  seurat$nCount_RNA >= 500 &
  seurat$percent.mt < mt_cut &
  seurat$percent.rb < rb_cut

removed_fraction <- 1 - mean(keep_cells)
if (removed_fraction > 0.30) {
  saveRDS(seurat, file.path(proj_dir, "seurat_before_automatic_QC.rds"))
  stop(
    "Automatic QC would remove ", round(100 * removed_fraction, 1),
    "% of cells. Inspect the saved QC TIFFs; do not silently relax thresholds."
  )
}

seurat <- subset(seurat, cells = colnames(seurat)[keep_cells])

p_qc_after <- VlnPlot(
  seurat,
  features = qc_features,
  ncol = 4,
  pt.size = 0
)
save_tif(p_qc_after, "QC_violin_after_filtering.tif", width = 12, height = 5)

# ---- 5. SCTransform with mitochondrial/ribosomal regression -----------------
seurat <- SCTransform(
  seurat,
  assay = "RNA",
  new.assay.name = "SCT",
  vst.flavor = "v2",
  vars.to.regress = c("percent.mt", "percent.rb"),
  variable.features.n = 3000,
  verbose = TRUE
)
DefaultAssay(seurat) <- "SCT"

top_hvg <- head(VariableFeatures(seurat), 20)
p_hvg <- VariableFeaturePlot(seurat)
p_hvg <- LabelPoints(
  plot = p_hvg,
  points = top_hvg,
  repel = TRUE
)
save_tif(p_hvg, "SCTransform_variable_features.tif", width = 8, height = 6)

# ---- 6. PCA and automatic elbow-based ndims --------------------------------
npcs_requested <- min(
  50L,
  ncol(seurat) - 1L,
  length(VariableFeatures(seurat)) - 1L
)
npcs_requested <- max(2L, as.integer(npcs_requested))

seurat <- RunPCA(
  seurat,
  assay = "SCT",
  features = VariableFeatures(seurat),
  npcs = npcs_requested,
  verbose = FALSE
)

infer_elbow <- function(pca_sd, min_dims = 10L, max_dims = 50L) {
  pca_sd <- as.numeric(pca_sd)
  max_dims <- min(max_dims, length(pca_sd))
  min_dims <- min(max(min_dims, 2L), max_dims)

  scaled_sd <- (pca_sd - min(pca_sd)) /
    max(max(pca_sd) - min(pca_sd), .Machine$double.eps)
  if (length(scaled_sd) < 4L) {
    return(as.integer(min_dims))
  }

  # Light smoothing prevents one noisy PC from defining the elbow.
  smoothed <- as.numeric(stats::filter(
    scaled_sd,
    filter = rep(1 / 3, 3),
    sides = 2
  ))
  smoothed[!is.finite(smoothed)] <- scaled_sd[!is.finite(smoothed)]

  curvature <- diff(diff(smoothed))
  eligible <- seq.int(min_dims - 1L, length(curvature))
  eligible <- eligible[
    eligible >= 1L & eligible <= length(curvature)
  ]
  if (length(eligible) == 0L) {
    return(as.integer(min_dims))
  }
  knee_from_curvature <- eligible[which.max(curvature[eligible])] + 1L

  as.integer(max(min_dims, min(knee_from_curvature, max_dims)))
}

pca_sd <- Stdev(seurat, reduction = "pca")
ndims <- infer_elbow(pca_sd, min_dims = 10L, max_dims = npcs_requested)

p_elbow <- ElbowPlot(
  seurat,
  reduction = "pca",
  ndims = npcs_requested,
  plot_type = "stdev"
) +
  geom_vline(xintercept = ndims, linetype = "dashed", color = "#E64B35") +
  labs(title = paste0("PCA elbow: selected ", ndims, " dimensions"))
save_tif(p_elbow, "PCA_elbow_selected_ndims.tif", width = 7, height = 5)

# ---- 7. UMAP and Leiden clustering at resolution 0.5 -----------------------
seurat <- FindNeighbors(
  seurat,
  reduction = "pca",
  dims = 1:ndims,
  graph.name = "sct_snn",
  verbose = FALSE
)

seurat <- RunUMAP(
  seurat,
  reduction = "pca",
  dims = 1:ndims,
  reduction.name = "umap",
  seed.use = 42,
  verbose = FALSE
)

seurat <- FindClusters(
  seurat,
  graph.name = "sct_snn",
  resolution = 0.5,
  algorithm = 4,
  leiden_method = "igraph",
  leiden_objective_function = "modularity",
  cluster.name = "leiden_res0.5",
  random.seed = 0,
  verbose = TRUE
)
Idents(seurat) <- "leiden_res0.5"

n_clusters <- length(levels(Idents(seurat)))
cluster_colors <- grDevices::hcl.colors(
  max(n_clusters, 2L),
  palette = "Dark 3"
)
p_umap_clusters <- DimPlot(
  seurat,
  reduction = "umap",
  group.by = "leiden_res0.5",
  cols = cluster_colors,
  label = TRUE,
  repel = TRUE
) +
  ggtitle(paste0("Leiden 0.5 | ", n_clusters, " clusters"))
save_tif(p_umap_clusters, "UMAP_Leiden_resolution_0.5.tif", width = 8, height = 6)

# Save the completed base analysis before a network-dependent annotation step.
saveRDS(seurat, file.path(proj_dir, "seurat_before_azimuth.rds"))

# ---- 8. AzimuthAPI annotation ------------------------------------------------
azimuth_mode <- "not_run"

if (species != "human") {
  writeLines(
    c(
      paste0("species\t", species),
      "azimuth_status\tskipped: AzimuthAPI is Pan-human only",
      paste0("ndims\t", ndims),
      "leiden_resolution\t0.5"
    ),
    file.path(proj_dir, "analysis_manifest.tsv")
  )
  stop(
    "AzimuthAPI is a Pan-human model and cannot annotate ", species,
    " data. Base analysis results were saved before this stop."
  )
}

if (!requireNamespace("AzimuthAPI", quietly = TRUE)) {
  stop(
    "AzimuthAPI is required for automatic annotation. ",
    "Install with remotes::install_github('satijalab/AzimuthAPI')."
  )
}

# AzimuthAPI expects normalized RNA data even though SCT remains the default
# assay used for the unsupervised UMAP and Leiden analysis above.
seurat <- NormalizeData(
  seurat,
  assay = "RNA",
  normalization.method = "LogNormalize",
  verbose = FALSE
)
DefaultAssay(seurat) <- "RNA"

annotate_with_azimuth <- function(object) {
  cloud_error <- NULL

  annotated <- tryCatch(
    AzimuthAPI::CloudAzimuth(object),
    error = function(e) {
      cloud_error <<- conditionMessage(e)
      NULL
    }
  )
  if (!is.null(annotated)) {
    return(list(object = annotated, mode = "CloudAzimuth"))
  }

  local_available <- requireNamespace("reticulate", quietly = TRUE) &&
    isTRUE(tryCatch(
      reticulate::py_module_available("panhumanpy"),
      error = function(e) FALSE
    ))

  if (local_available) {
    annotated <- AzimuthAPI::ANNotate(object)
    return(list(object = annotated, mode = "ANNotate"))
  }

  stop(
    "CloudAzimuth failed and local panhumanpy is unavailable. Cloud error: ",
    cloud_error
  )
}

azimuth_result <- annotate_with_azimuth(seurat)
seurat <- azimuth_result$object
azimuth_mode <- azimuth_result$mode
DefaultAssay(seurat) <- "SCT"

if (!"final_level_confidence" %in% colnames(seurat[[]])) {
  stop("Azimuth completed but final_level_confidence is missing.")
}
if (!"azimuth_embed" %in% Reductions(seurat)) {
  stop("Azimuth completed but azimuth_embed is missing.")
}

seurat$azimuth_confidence_pass <- seurat$final_level_confidence > 0.5

if ("azimuth_fine" %in% colnames(seurat[[]]) &&
    "PrepLabel" %in% getNamespaceExports("AzimuthAPI")) {
  seurat <- AzimuthAPI::PrepLabel(
    seurat,
    "azimuth_fine",
    "azimuth_fine_filtered",
    cutoff = 20
  )
}

fine_label <- if ("azimuth_fine_filtered" %in% colnames(seurat[[]])) {
  "azimuth_fine_filtered"
} else {
  "azimuth_fine"
}

if (fine_label %in% colnames(seurat[[]])) {
  seurat$azimuth_label_auto <- as.character(seurat[[fine_label]][, 1])
  seurat$azimuth_label_auto[!seurat$azimuth_confidence_pass] <-
    "Unassigned_low_confidence"
}

embedding_dims <- min(128L, ncol(Embeddings(seurat[["azimuth_embed"]])))
seurat <- RunUMAP(
  seurat,
  reduction = "azimuth_embed",
  dims = 1:embedding_dims,
  reduction.name = "azimuth_umap",
  seed.use = 42,
  verbose = FALSE
)

azimuth_cols <- grDevices::hcl.colors(
  max(length(unique(seurat$azimuth_medium)), 2L),
  palette = "Dark 3"
)
p_umap_medium <- DimPlot(
  seurat,
  reduction = "umap",
  group.by = "azimuth_medium",
  cols = azimuth_cols,
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Azimuth medium labels on unsupervised UMAP")
save_tif(p_umap_medium, "UMAP_Azimuth_medium_labels.tif", width = 9, height = 6)

if ("azimuth_label_auto" %in% colnames(seurat[[]])) {
  auto_cols <- grDevices::hcl.colors(
    max(length(unique(seurat$azimuth_label_auto)), 2L),
    palette = "Dark 3"
  )
  p_umap_fine <- DimPlot(
    seurat,
    reduction = "azimuth_umap",
    group.by = "azimuth_label_auto",
    cols = auto_cols,
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle("Azimuth fine labels with low-confidence cells retained")
  save_tif(p_umap_fine, "UMAP_Azimuth_fine_confidence_aware.tif", width = 10, height = 6)
}

p_confidence <- FeaturePlot(
  seurat,
  features = "final_level_confidence",
  reduction = "azimuth_umap",
  min.cutoff = 0,
  max.cutoff = 1
) +
  ggtitle("Azimuth final-level confidence")
save_tif(p_confidence, "UMAP_Azimuth_confidence.tif", width = 8, height = 6)

if ("make_azimuth_QC_heatmaps" %in% getNamespaceExports("AzimuthAPI")) {
  azimuth_heatmaps <- tryCatch(
    AzimuthAPI::make_azimuth_QC_heatmaps(seurat),
    error = function(e) {
      warning("Azimuth QC heatmaps failed: ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(azimuth_heatmaps)) {
    for (i in seq_along(azimuth_heatmaps)) {
      label <- names(azimuth_heatmaps)[[i]]
      label <- if (is.null(label) || !nzchar(label)) paste0("panel_", i) else label
      label <- gsub("[^A-Za-z0-9_-]+", "_", label)
      save_tif(
        azimuth_heatmaps[[i]],
        paste0("Azimuth_QC_heatmap_", label, ".tif"),
        width = 10,
        height = 8
      )
    }
  }
}

# ---- 9. Save object, manifest, and session information ----------------------
saveRDS(seurat, file.path(proj_dir, "seurat_basic_analysis.rds"))

manifest <- data.frame(
  key = c(
    "input_path", "species", "human_feature_score", "mouse_feature_score",
    "nFeature_RNA_min", "nCount_RNA_min", "percent_mt_max", "percent_rb_max",
    "sc_transform_flavor", "sc_transform_regression", "ndims", "leiden_algorithm",
    "leiden_resolution", "azimuth_mode", "tiff_dpi", "tiff_compression"
  ),
  value = c(
    input_path, species, species_info$human_score, species_info$mouse_score,
    200, 500, round(mt_cut, 3), round(rb_cut, 3),
    "v2", "percent.mt + percent.rb", ndims, "Leiden (algorithm=4)",
    0.5, azimuth_mode, 300, "lzw"
  ),
  stringsAsFactors = FALSE
)
write.table(
  manifest,
  file = file.path(proj_dir, "analysis_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

capture.output(sessionInfo(), file = file.path(proj_dir, "sessionInfo.txt"))

message("Basic Seurat v6 analysis complete.")
message("TIFF figures: ", fig_dir)
message("Seurat object: ", file.path(proj_dir, "seurat_basic_analysis.rds"))
