# ============================================================
# 04_CD45neg.R
# CD45-negative annotation
#
# Starting objects:
#   fibroblasts.QC.rds
#   fibroblasts.Fibros.annotated.rds
#
# Annotation columns:
#
#   Clusters
#     Generic CD45-negative annotation:
#       Endothelial Vascular
#       Aerocytes
#       Epithelial AT1
#       Epithelial AT2
#       Epithelial Krt8+
#       Ciliated
#       Mesenchymal
#
#   fibro_annotation
#     Same complete annotation, replacing Mesenchymal with:
#       Fibroblasts Universal
#       Fibroblasts Activated 1
#       Fibroblasts Activated 2
#       Mesothelial
#
#   cancer_annotation
#     Same complete annotation, replacing tumour epithelial cells:
#       cluster 12                  -> Prolif. Cancer Cells
#       remaining Epithelial Krt8+ -> EMT Cancer Cells
#
# Naive Epithelial Krt8+ cells are removed so that no cancer
# cells remain in naive samples.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SCINA)
  library(scCustomize)
  library(ggplot2)
  library(cowplot)
})

set.seed(42)

# ============================================================
# Config
# ============================================================

config <- read.table(
  "config.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

keys <- config$key
values <- config$value
names(values) <- keys

rsdir <- values["rsdir"]
datadir <- values["datadir"]
outdir <- values["outdir"]

dir.create(
  rsdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(outdir, "objects"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(outdir, "CD45neg"),
  recursive = TRUE,
  showWarnings = FALSE
)

source("R/Functions.R")

# ============================================================
# Read objects
# ============================================================

fibroblasts.QC <- readRDS(
  file.path(
    outdir,
    "objects",
    "fibroblasts.QC.rds"
  )
)

fibroblasts.Fibros <- readRDS(
  file.path(
    outdir,
    "objects",
    "fibroblasts.Fibros.annotated.rds"
  )
)

DefaultAssay(fibroblasts.QC) <- "RNA"

# ============================================================
# Select CD45-negative cells
# ============================================================

fibros.CD45NP <- subset(
  fibroblasts.QC,
  subset = CD == "CD45+",
  invert = TRUE
)

cat(
  "\nCD45-negative cells before filtering:",
  ncol(fibros.CD45NP),
  "\n"
)

# ============================================================
# Remove unwanted global clusters
# ============================================================

excluded_clusters <- c(
  0, 1, 4, 5, 6, 7, 8, 10,
  14, 15, 19, 21, 22, 24, 25
)

cells_keep <- colnames(fibros.CD45NP)[
  !as.character(
    fibros.CD45NP$seurat_clusters
  ) %in% as.character(excluded_clusters)
]

fibros.CD45NP <- subset(
  fibros.CD45NP,
  cells = cells_keep
)

# ============================================================
# Apply original global UMAP boundary
# Keep cells with UMAP_1 <= 7
# ============================================================

umap_coordinates <- Embeddings(
  fibros.CD45NP,
  reduction = "umap"
)

cells_keep <- rownames(umap_coordinates)[
  umap_coordinates[, 1] <= 7
]

fibros.CD45NP <- subset(
  fibros.CD45NP,
  cells = cells_keep
)

cat(
  "CD45-negative cells after initial filtering:",
  ncol(fibros.CD45NP),
  "\n"
)

print(
  table(
    fibros.CD45NP$seurat_clusters
  )
)

# ============================================================
# Prepare normalized RNA matrix for SCINA
# ============================================================

DefaultAssay(fibros.CD45NP) <- "RNA"

if (inherits(fibros.CD45NP[["RNA"]], "Assay5")) {
  fibros.CD45NP <- JoinLayers(
    fibros.CD45NP,
    assay = "RNA"
  )
}

expression_matrix <- GetAssayData(
  fibros.CD45NP,
  assay = "RNA",
  layer = "data"
)

expression_matrix <- as.matrix(
  expression_matrix
)

# ============================================================
# Read SCINA signatures
# ============================================================

signature <- read.table(
  file.path(
    datadir,
    "signatures",
    "signature.endothelial.tsv"
  ),
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ============================================================
# Run SCINA
# ============================================================

scina_results <- SCINA(
  expression_matrix,
  signature,
  max_iter = 100,
  convergence_n = 10,
  convergence_rate = 0.999,
  sensitivity_cutoff = 0.9,
  rm_overlap = FALSE,
  allow_unknown = TRUE,
  log_file = file.path(
    outdir,
    "CD45neg",
    "SCINA.log"
  )
)

fibros.CD45NP$SCINA <- scina_results$cell_labels

cat("\nSCINA annotation:\n")

print(
  table(
    fibros.CD45NP$SCINA,
    useNA = "ifany"
  )
)

# ============================================================
# Remove platelet cluster
# ============================================================

fibros.CD45NP <- subset(
  fibros.CD45NP,
  subset = seurat_clusters == 23,
  invert = TRUE
)

# ============================================================
# Fibroblast labels from independent fibroblast analysis
# ============================================================

fibro_labels <- setNames(
  as.character(
    fibroblasts.Fibros$fibro_annotation
  ),
  colnames(fibroblasts.Fibros)
)

valid_fibro_cells <- names(fibro_labels)[
  !is.na(fibro_labels)
]

# ============================================================
# Remove cells discarded during fibroblast reclustering
#
# Global clusters 17 and 20 contained endothelial contaminants.
# Only cells retained in fibroblasts.Fibros.annotated.rds
# are kept.
# ============================================================

global_mesenchymal_cells <- colnames(fibros.CD45NP)[
  as.character(
    fibros.CD45NP$seurat_clusters
  ) %in% c("17", "20")
]

discarded_mesenchymal_cells <- setdiff(
  global_mesenchymal_cells,
  valid_fibro_cells
)

cat(
  "\nCells removed from global clusters 17/20:",
  length(discarded_mesenchymal_cells),
  "\n"
)

fibros.CD45NP <- subset(
  fibros.CD45NP,
  cells = setdiff(
    colnames(fibros.CD45NP),
    discarded_mesenchymal_cells
  )
)

# ============================================================
# Generic CD45-negative annotation
# ============================================================

fibros.CD45NP$Clusters <- NA_character_

# Endothelial vascular
fibros.CD45NP$Clusters[
  as.character(
    fibros.CD45NP$seurat_clusters
  ) %in% c("2", "3", "18")
] <- "Endothelial Vascular"

# Aerocytes
fibros.CD45NP$Clusters[
  as.character(
    fibros.CD45NP$seurat_clusters
  ) %in% c("9", "11")
] <- "Aerocytes"

# Mesenchymal compartment
fibros.CD45NP$Clusters[
  as.character(
    fibros.CD45NP$seurat_clusters
  ) %in% c("17", "20")
] <- "Mesenchymal"

# Krt8-positive epithelial compartment
fibros.CD45NP$Clusters[
  as.character(
    fibros.CD45NP$seurat_clusters
  ) %in% c("12", "13")
] <- "Epithelial Krt8+"

# Cluster 16 resolved with SCINA
fibros.CD45NP$Clusters[
  as.character(
    fibros.CD45NP$seurat_clusters
  ) == "16" &
    fibros.CD45NP$SCINA == "Epithelial AT1"
] <- "Epithelial AT1"

fibros.CD45NP$Clusters[
  as.character(
    fibros.CD45NP$seurat_clusters
  ) == "16" &
    fibros.CD45NP$SCINA == "Epithelial AT2"
] <- "Epithelial AT2"

fibros.CD45NP$Clusters[
  as.character(
    fibros.CD45NP$seurat_clusters
  ) == "16" &
    fibros.CD45NP$SCINA == "unknown"
] <- "Epithelial AT2"

fibros.CD45NP$Clusters[
  as.character(
    fibros.CD45NP$seurat_clusters
  ) == "16" &
    fibros.CD45NP$SCINA == "Ciliated"
] <- "Ciliated"

# Keep only cells with a generic annotation
fibros.CD45NP <- subset(
  fibros.CD45NP,
  cells = colnames(fibros.CD45NP)[
    !is.na(fibros.CD45NP$Clusters)
  ]
)

# ============================================================
# Complete annotation with fibroblast subtypes
#
# Start from Clusters and replace Mesenchymal with:
#   Fibroblasts Universal
#   Fibroblasts Activated 1
#   Fibroblasts Activated 2
#   Mesothelial
# ============================================================

fibros.CD45NP$fibro_annotation <- as.character(
  fibros.CD45NP$Clusters
)

mesenchymal_cells <- colnames(fibros.CD45NP)[
  as.character(
    fibros.CD45NP$Clusters
  ) == "Mesenchymal"
]

common_fibro_cells <- intersect(
  mesenchymal_cells,
  valid_fibro_cells
)

fibros.CD45NP$fibro_annotation[
  match(
    common_fibro_cells,
    colnames(fibros.CD45NP)
  )
] <- unname(
  fibro_labels[common_fibro_cells]
)

cat(
  "\nMesenchymal cells with fine annotation:",
  length(common_fibro_cells),
  "of",
  length(mesenchymal_cells),
  "\n"
)

# Check that all retained mesenchymal cells received
# a fibroblast or mesothelial identity
remaining_mesenchymal <- sum(
  fibros.CD45NP$fibro_annotation == "Mesenchymal",
  na.rm = TRUE
)

if (remaining_mesenchymal > 0) {
  stop(
    remaining_mesenchymal,
    " mesenchymal cells lack a final fibroblast annotation."
  )
}

# ============================================================
# Remove Krt8-positive epithelial cells from naive samples
#
# These clusters correspond to cancer cells in tumour samples.
# They are removed from naive samples so that no cancer-cell
# population remains in naive conditions.
# ============================================================

naive_krt8_cells <- colnames(fibros.CD45NP)[
  as.character(
    fibros.CD45NP$group
  ) == "Naive" &
    as.character(
      fibros.CD45NP$fibro_annotation
    ) == "Epithelial Krt8+"
]

cat(
  "\nEpithelial Krt8+ cells removed from naive samples:",
  length(naive_krt8_cells),
  "\n"
)

fibros.CD45NP <- subset(
  fibros.CD45NP,
  cells = setdiff(
    colnames(fibros.CD45NP),
    naive_krt8_cells
  )
)

# ============================================================
# Complete annotation including cancer cells
#
# Start from fibro_annotation.
#
# After removing naive Krt8-positive cells:
#   cluster 12                  -> Prolif. Cancer Cells
#   remaining Epithelial Krt8+ -> EMT Cancer Cells
# ============================================================

fibros.CD45NP$cancer_annotation <- as.character(
  fibros.CD45NP$fibro_annotation
)

global_clusters <- as.character(
  fibros.CD45NP$seurat_clusters
)

# Cluster 12: proliferating cancer cells
fibros.CD45NP$cancer_annotation[
  global_clusters == "12"
] <- "Prolif. Cancer Cells"

# All remaining Epithelial Krt8+ cells become EMT cancer cells
fibros.CD45NP$cancer_annotation[
  fibros.CD45NP$cancer_annotation ==
    "Epithelial Krt8+"
] <- "EMT Cancer Cells"

# ============================================================
# Annotation levels
# ============================================================

fibros.CD45NP$Clusters <- factor(
  fibros.CD45NP$Clusters,
  levels = c(
    "Endothelial Vascular",
    "Aerocytes",
    "Epithelial AT1",
    "Epithelial AT2",
    "Epithelial Krt8+",
    "Ciliated",
    "Mesenchymal"
  )
)

fibros.CD45NP$fibro_annotation <- factor(
  fibros.CD45NP$fibro_annotation,
  levels = c(
    "Endothelial Vascular",
    "Aerocytes",
    "Epithelial AT1",
    "Epithelial AT2",
    "Epithelial Krt8+",
    "Ciliated",
    "Fibroblasts Universal",
    "Fibroblasts Activated 1",
    "Fibroblasts Activated 2",
    "Mesothelial"
  )
)

fibros.CD45NP$cancer_annotation <- factor(
  fibros.CD45NP$cancer_annotation,
  levels = c(
    "Endothelial Vascular",
    "Aerocytes",
    "Epithelial AT1",
    "Epithelial AT2",
    "Ciliated",
    "Fibroblasts Universal",
    "Fibroblasts Activated 1",
    "Fibroblasts Activated 2",
    "Mesothelial",
    "Prolif. Cancer Cells",
    "EMT Cancer Cells"
  )
)

fibros.CD45NP$orig.ident <- factor(
  as.character(
    fibros.CD45NP$orig.ident
  ),
  levels = c(
    "WT",
    "CD169DTR",
    "CCR2KO"
  )
)

# ============================================================
# Final annotation checks
# ============================================================

remaining_krt8 <- sum(
  as.character(
    fibros.CD45NP$cancer_annotation
  ) == "Epithelial Krt8+",
  na.rm = TRUE
)

if (remaining_krt8 > 0) {
  stop(
    remaining_krt8,
    " Epithelial Krt8+ cells remain in cancer_annotation."
  )
}

cancer_labels <- c(
  "Prolif. Cancer Cells",
  "EMT Cancer Cells"
)

naive_cancer_cells <- (
  as.character(
    fibros.CD45NP$group
  ) == "Naive"
) & (
  as.character(
    fibros.CD45NP$cancer_annotation
  ) %in% cancer_labels
)

if (any(naive_cancer_cells)) {
  stop(
    "Cancer-cell annotations remain in naive samples."
  )
}

# ============================================================
# Final tables
# ============================================================

cat("\nGeneric CD45-negative annotation:\n")

print(
  table(
    fibros.CD45NP$Clusters,
    useNA = "ifany"
  )
)

cat("\nAnnotation with fibroblast subtypes:\n")

print(
  table(
    fibros.CD45NP$fibro_annotation,
    useNA = "ifany"
  )
)

cat("\nFinal annotation including cancer cells:\n")

print(
  table(
    fibros.CD45NP$cancer_annotation,
    useNA = "ifany"
  )
)

cat("\nFinal annotation by experimental group:\n")

print(
  table(
    group = fibros.CD45NP$group,
    annotation = fibros.CD45NP$cancer_annotation,
    useNA = "ifany"
  )
)

cat("\nFinal annotation by condition:\n")

print(
  table(
    condition = fibros.CD45NP$orig.ident,
    annotation = fibros.CD45NP$cancer_annotation,
    useNA = "ifany"
  )
)

cat("\nFinal annotation by condition and group:\n")

print(
  table(
    condition = fibros.CD45NP$orig.ident,
    group = fibros.CD45NP$group,
    annotation = fibros.CD45NP$cancer_annotation,
    useNA = "ifany"
  )
)

# ============================================================
# Final check plots
# ============================================================

p_generic <- DimPlot_scCustom(
  fibros.CD45NP,
  reduction = "umap",
  group.by = "Clusters",
  split.by = "orig.ident",
  label = TRUE,
  repel = TRUE,
  label.size = 3.2,
  ggplot_default_colors = TRUE,
  pt.size = 0.4
) +
  ggtitle(
    "Generic CD45-negative annotation"
  )

p_fibro <- DimPlot_scCustom(
  fibros.CD45NP,
  reduction = "umap",
  group.by = "fibro_annotation",
  split.by = "orig.ident",
  label = TRUE,
  repel = TRUE,
  label.size = 3.2,
  ggplot_default_colors = TRUE,
  pt.size = 0.4
) +
  ggtitle(
    "Annotation with fibroblast subtypes"
  )

p_cancer <- DimPlot_scCustom(
  fibros.CD45NP,
  reduction = "umap",
  group.by = "cancer_annotation",
  split.by = "orig.ident",
  label = TRUE,
  repel = TRUE,
  label.size = 3.2,
  ggplot_default_colors = TRUE,
  pt.size = 0.4
) +
  ggtitle(
    "Final annotation including cancer cells"
  )

p_final <- plot_grid(
  p_generic,
  p_fibro,
  p_cancer,
  ncol = 1
)

print(p_final)

ggsave(
  filename = file.path(
    outdir,
    "CD45neg",
    "CD45neg_final_annotation_UMAP.pdf"
  ),
  plot = p_final,
  width = 18,
  height = 18,
  limitsize = FALSE
)

# ============================================================
# Final identity
# ============================================================

Idents(fibros.CD45NP) <- "cancer_annotation"

# ============================================================
# Save final object
# ============================================================

saveRDS(
  fibros.CD45NP,
  file.path(
    outdir,
    "objects",
    "fibroblasts.CD45NP.annotated.rds"
  )
)

cat(
  "\nObject saved in:\n",
  file.path(
    outdir,
    "objects",
    "fibroblasts.CD45NP.annotated.rds"
  ),
  "\n"
)