# ============================================================
# 05_CancerCells.R
# Annotation and cleaning of cancer cells
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(scCustomize)
  library(ggplot2)
})

set.seed(42)

# ============================================================
# Configuration
# ============================================================

config <- read.table(
  "config.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

cfg <- setNames(
  config$value,
  config$key
)

rsdir <- cfg[["rsdir"]]
outdir <- cfg[["outdir"]]

objects_dir <- file.path(
  rsdir,
  "objects"
)

plots_dir <- file.path(
  outdir,
  "Cancer.cells"
)

dir.create(
  objects_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  plots_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

source("R/Functions.R")

# ============================================================
# Load annotated CD45-negative object
# ============================================================

input_file <- file.path(
  objects_dir,
  "fibroblasts.CD45NP.annotated.rds"
)

if (!file.exists(input_file)) {
  stop(
    "Input object not found: ",
    input_file
  )
}

data <- readRDS(input_file)

cat(
  "\nLoaded object:",
  input_file,
  "\n"
)

cat(
  "Initial number of cells:",
  ncol(data),
  "\n"
)

# ============================================================
# Check required metadata
# ============================================================

required_metadata <- c(
  "seurat_clusters",
  "fibroblast_annotation",
  "orig.ident",
  "group"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(data@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

# Preserve the clusters from the complete CD45-negative analysis.

data$global_clusters <- as.character(
  data$seurat_clusters
)

cat("\nOriginal global clusters:\n")

print(
  table(
    data$global_clusters,
    useNA = "ifany"
  )
)

# ============================================================
# Check expected cancer-cell clusters
# ============================================================

required_cancer_clusters <- c(
  "12",
  "13"
)

missing_cancer_clusters <- setdiff(
  required_cancer_clusters,
  unique(data$global_clusters)
)

if (length(missing_cancer_clusters) > 0) {
  stop(
    "Expected global cancer clusters not found: ",
    paste(missing_cancer_clusters, collapse = ", ")
  )
}

# ============================================================
# Cancer-cell marker plots
# ============================================================

DefaultAssay(data) <- "RNA"

cancer_markers <- c(
  "Epcam",
  "Krt8",
  "Krt18",
  "Mki67",
  "Top2a",
  "Sftpc",
  "Lgals3",
  "Hmga2"
)

cancer_markers_present <- intersect(
  cancer_markers,
  rownames(data)
)

cancer_markers_missing <- setdiff(
  cancer_markers,
  cancer_markers_present
)

if (length(cancer_markers_missing) > 0) {
  warning(
    "Cancer markers absent from the object: ",
    paste(cancer_markers_missing, collapse = ", ")
  )
}

if (length(cancer_markers_present) > 0) {
  p_cancer_markers <- FeaturePlot_scCustom(
    data,
    features = cancer_markers_present
  )

  ggsave(
    filename = file.path(
      plots_dir,
      "cancer_cells_markers.pdf"
    ),
    plot = p_cancer_markers,
    width = 18,
    height = 12,
    limitsize = FALSE
  )
}

plasticity_genes <- c(
  "Hmga2",
  "Car9",
  "Muc1",
  "S100a6",
  "Clu",
  "Mki67",
  "Cpe",
  "Lgals3",
  "Vim",
  "Fn1",
  "Zeb1"
)

plasticity_genes_present <- intersect(
  plasticity_genes,
  rownames(data)
)

plasticity_genes_missing <- setdiff(
  plasticity_genes,
  plasticity_genes_present
)

if (length(plasticity_genes_missing) > 0) {
  warning(
    "Plasticity genes absent from the object: ",
    paste(plasticity_genes_missing, collapse = ", ")
  )
}

if (length(plasticity_genes_present) > 0) {
  p_plasticity <- FeaturePlot_scCustom(
    data,
    features = plasticity_genes_present
  )

  ggsave(
    filename = file.path(
      plots_dir,
      "cancer_cells_plasticity.pdf"
    ),
    plot = p_plasticity,
    width = 24,
    height = 16,
    limitsize = FALSE
  )

  p_plasticity_split <- FeaturePlot_scCustom(
    data,
    features = plasticity_genes_present,
    split.by = "orig.ident"
  )

  ggsave(
    filename = file.path(
      plots_dir,
      "cancer_cells_plasticity_by_condition.png"
    ),
    plot = p_plasticity_split,
    width = 16,
    height = 36,
    units = "in",
    dpi = 150,
    limitsize = FALSE
  )
}

# ============================================================
# Annotate the global cancer-cell clusters
#
# Global cluster 12 = proliferating cancer cells
# Global cluster 13 = EMT/Hmga2+ cancer cells
# ============================================================

data$fibroblast_annotation <- as.character(
  data$fibroblast_annotation
)

data$fibroblast_annotation[
  data$global_clusters == "12"
] <- "Prolif. Cancer Cells"

data$fibroblast_annotation[
  data$global_clusters == "13"
] <- "EMT Cancer Cells"

cat("\nAnnotation after identifying cancer cells:\n")

print(
  table(
    global_cluster = data$global_clusters,
    annotation = data$fibroblast_annotation,
    useNA = "ifany"
  )
)

cancer_cells_initial <- data$global_clusters %in% c(
  "12",
  "13"
)

cat("\nCancer cells by condition and group:\n")

print(
  table(
    condition = data$orig.ident[cancer_cells_initial],
    group = data$group[cancer_cells_initial],
    annotation = data$fibroblast_annotation[cancer_cells_initial],
    useNA = "ifany"
  )
)

# ============================================================
# Recluster Hmga2+/EMT cancer cells
# ============================================================

Idents(data) <- "fibroblast_annotation"

hmga2_cells <- subset(
  data,
  idents = "EMT Cancer Cells"
)

if (ncol(hmga2_cells) == 0) {
  stop(
    "No EMT Cancer Cells were found in global cluster 13."
  )
}

cat(
  "\nEMT cancer cells before reclustering:",
  ncol(hmga2_cells),
  "\n"
)

print(
  table(
    condition = hmga2_cells$orig.ident,
    group = hmga2_cells$group,
    useNA = "ifany"
  )
)

DefaultAssay(hmga2_cells) <- "RNA"

# Join RNA layers when using Seurat v5.

if (inherits(hmga2_cells[["RNA"]], "Assay5")) {
  hmga2_cells <- JoinLayers(
    hmga2_cells,
    assay = "RNA"
  )
}

# Original reclustering workflow.

hmga2_cells <- FindVariableFeatures(
  hmga2_cells
)

hmga2_cells <- ScaleData(
  hmga2_cells
)

hmga2_cells <- RunPCA(
  hmga2_cells
)

# ============================================================
# PCA inspection
# Original analysis selected 10 PCs
# ============================================================

number_pcs <- ncol(
  Embeddings(
    hmga2_cells,
    reduction = "pca"
  )
)

p_elbow <- ElbowPlot(
  hmga2_cells,
  ndims = min(30, number_pcs)
)

ggsave(
  filename = file.path(
    plots_dir,
    "EMT_cancer_cells_PCA_elbow.pdf"
  ),
  plot = p_elbow,
  width = 7,
  height = 5
)

if (number_pcs < 10) {
  stop(
    "Fewer than 10 principal components were generated."
  )
}

# ============================================================
# Neighbors, clustering and UMAP
#
# Original parameters:
#   PCs: 1:10
#   resolution: 0.1
# ============================================================

hmga2_cells <- FindNeighbors(
  hmga2_cells,
  dims = 1:10
)

hmga2_cells <- FindClusters(
  hmga2_cells,
  resolution = 0.1
)

hmga2_cells <- RunUMAP(
  hmga2_cells,
  dims = 1:10
)

hmga2_cells$hmga2_cluster <- as.character(
  hmga2_cells$seurat_clusters
)

cat("\nEMT internal reclustering:\n")

print(
  table(
    hmga2_cells$hmga2_cluster,
    useNA = "ifany"
  )
)

cat("\nEMT clusters by condition:\n")

print(
  table(
    condition = hmga2_cells$orig.ident,
    cluster = hmga2_cells$hmga2_cluster,
    useNA = "ifany"
  )
)

cat("\nEMT clusters by group:\n")

print(
  table(
    group = hmga2_cells$group,
    cluster = hmga2_cells$hmga2_cluster,
    useNA = "ifany"
  )
)

# ============================================================
# Plot EMT reclustering before filtering
# ============================================================

p_hmga2_clusters <- DimPlot_scCustom(
  hmga2_cells,
  reduction = "umap",
  group.by = "hmga2_cluster",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5
) +
  ggtitle("EMT cancer-cell reclustering")

ggsave(
  filename = file.path(
    plots_dir,
    "EMT_cancer_cells_reclustered_UMAP.pdf"
  ),
  plot = p_hmga2_clusters,
  width = 8,
  height = 6
)

p_hmga2_clusters_split <- DimPlot_scCustom(
  hmga2_cells,
  reduction = "umap",
  group.by = "hmga2_cluster",
  split.by = "orig.ident",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5
) +
  ggtitle("EMT cancer-cell reclustering")

ggsave(
  filename = file.path(
    plots_dir,
    "EMT_cancer_cells_reclustered_by_condition.pdf"
  ),
  plot = p_hmga2_clusters_split,
  width = 18,
  height = 6,
  limitsize = FALSE
)

# ============================================================
# Remove contaminant clusters and Naive cells
#
# Internal EMT clusters:
#   3 and 5 = epithelial/blood contaminants
# ============================================================

expected_contaminant_clusters <- c(
  "3",
  "5"
)

missing_contaminant_clusters <- setdiff(
  expected_contaminant_clusters,
  unique(hmga2_cells$hmga2_cluster)
)

if (length(missing_contaminant_clusters) > 0) {
  warning(
    "Expected contaminant clusters not found: ",
    paste(missing_contaminant_clusters, collapse = ", ")
  )
}

barcodes_cluster3 <- colnames(hmga2_cells)[
  hmga2_cells$hmga2_cluster == "3"
]

barcodes_cluster5 <- colnames(hmga2_cells)[
  hmga2_cells$hmga2_cluster == "5"
]

barcodes_naive <- colnames(hmga2_cells)[
  as.character(hmga2_cells$group) == "Naive"
]

barcodes_to_remove <- unique(c(
  barcodes_cluster3,
  barcodes_cluster5,
  barcodes_naive
))

cat("\nEMT filtering summary:\n")

cat(
  "Internal cluster 3:",
  length(barcodes_cluster3),
  "cells\n"
)

cat(
  "Internal cluster 5:",
  length(barcodes_cluster5),
  "cells\n"
)

cat(
  "Naive EMT-compartment cells:",
  length(barcodes_naive),
  "cells\n"
)

cat(
  "Unique cells removed:",
  length(barcodes_to_remove),
  "\n"
)

if (length(barcodes_to_remove) == ncol(hmga2_cells)) {
  stop(
    "Filtering would remove all EMT cancer cells."
  )
}

# ============================================================
# Create filtered EMT object
# ============================================================

cells_to_keep_hmga2 <- setdiff(
  colnames(hmga2_cells),
  barcodes_to_remove
)

hmga2_cells_filtered <- subset(
  hmga2_cells,
  cells = cells_to_keep_hmga2
)

cat(
  "Retained EMT cancer cells:",
  ncol(hmga2_cells_filtered),
  "\n"
)

print(
  table(
    condition = hmga2_cells_filtered$orig.ident,
    cluster = hmga2_cells_filtered$hmga2_cluster,
    useNA = "ifany"
  )
)

# ============================================================
# Remove the same cells from the complete object
# ============================================================

cells_to_keep_data <- setdiff(
  colnames(data),
  barcodes_to_remove
)

data <- subset(
  data,
  cells = cells_to_keep_data
)

cat(
  "\nFinal number of cells in the complete object:",
  ncol(data),
  "\n"
)

# ============================================================
# Final annotations
# ============================================================

# Keep every annotation present in the object.
# This avoids converting Mesothelial or other valid populations to NA.

data$fibroblast_annotation <- factor(
  as.character(data$fibroblast_annotation),
  levels = unique(
    as.character(data$fibroblast_annotation)
  )
)

hmga2_cells_filtered$cancer_annotation <- factor(
  rep(
    "EMT Cancer Cells",
    ncol(hmga2_cells_filtered)
  ),
  levels = c(
    "Prolif. Cancer Cells",
    "EMT Cancer Cells"
  )
)

cat("\nFinal complete annotation:\n")

print(
  table(
    data$fibroblast_annotation,
    useNA = "ifany"
  )
)

cat("\nFinal annotation by condition:\n")

print(
  table(
    condition = data$orig.ident,
    annotation = data$fibroblast_annotation,
    useNA = "ifany"
  )
)

cat("\nFinal annotation by condition and group:\n")

print(
  table(
    condition = data$orig.ident,
    group = data$group,
    annotation = data$fibroblast_annotation,
    useNA = "ifany"
  )
)

# ============================================================
# Plot retained EMT cancer cells
#
# The original reclustered UMAP coordinates are preserved.
# ============================================================

p_hmga2_filtered <- DimPlot_scCustom(
  hmga2_cells_filtered,
  reduction = "umap",
  group.by = "hmga2_cluster",
  split.by = "orig.ident",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5
) +
  ggtitle("Filtered EMT cancer cells")

ggsave(
  filename = file.path(
    plots_dir,
    "EMT_cancer_cells_filtered_UMAP.pdf"
  ),
  plot = p_hmga2_filtered,
  width = 18,
  height = 6,
  limitsize = FALSE
)

# ============================================================
# Final complete CD45-negative UMAP
# ============================================================

p_final <- DimPlot_scCustom(
  data,
  reduction = "umap",
  group.by = "fibroblast_annotation",
  split.by = "orig.ident",
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  pt.size = 0.4
) +
  ggtitle("Final CD45-negative annotation")

ggsave(
  filename = file.path(
    plots_dir,
    "CD45neg_final_annotation_UMAP.pdf"
  ),
  plot = p_final,
  width = 18,
  height = 6,
  limitsize = FALSE
)

# ============================================================
# Save objects
# ============================================================

Idents(hmga2_cells_filtered) <- "cancer_annotation"
Idents(data) <- "fibroblast_annotation"

hmga2_output_file <- file.path(
  objects_dir,
  "fibroblasts.CancerCells.EMT.filtered.rds"
)

final_output_file <- file.path(
  objects_dir,
  "fibroblasts.CD45NP.final.rds"
)

saveRDS(
  data,
  hmga2_output_file
)

saveRDS(
  data,
  final_output_file
)

cat("\nObjects saved:\n")
cat(hmga2_output_file, "\n")
cat(final_output_file, "\n")

cat("\n05_CancerCells.R completed successfully.\n")