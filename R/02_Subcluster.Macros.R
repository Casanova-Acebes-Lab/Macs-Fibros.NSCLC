# ============================================================
# Macrophage / monocyte / dendritic-cell subclustering
# and annotation
#
# Input:
#   objects/fibroblasts.QC.rds
#
# Main outputs:
#   objects/fibroblasts.Macrophage_DC.annotated.rds
#   annotations/barcodes.macros.rds
#   annotations/barcodes.macros.tsv
#
# Diagnostic plots:
#   diagnostics/Macrophage_DC/
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(SCINA)
})

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

# ============================================================
# Output directories
# ============================================================

objects_dir <- file.path(rsdir, "objects")
annotations_dir <- file.path(rsdir, "annotations")
tables_dir <- file.path(rsdir, "tables", "Macrophage_DC")
diagnostics_dir <- file.path(outdir, "diagnostics", "Macrophage_DC")

dir.create(objects_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(annotations_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Functions
# ============================================================

source("R/Functions.R")

if (!exists(
  "scale",
  envir = .GlobalEnv,
  inherits = FALSE,
  mode = "function"
)) {
  stop("The scale() function was not found in R/Functions.R")
}

if (!exists(
  "PCA",
  envir = .GlobalEnv,
  inherits = FALSE,
  mode = "function"
)) {
  stop("The PCA() function was not found in R/Functions.R")
}

# ============================================================
# Input files
# ============================================================

input_file <- file.path(
  objects_dir,
  "fibroblasts.QC.rds"
)

signature_file <- file.path(
  datadir,
  "signatures/Maria.gene.tsv"
)

if (!file.exists(input_file)) {
  stop("Input object not found: ", input_file)
}

if (!file.exists(signature_file)) {
  stop("SCINA signature file not found: ", signature_file)
}

# ============================================================
# Read QC object
# ============================================================

fibroblasts <- readRDS(input_file)

DefaultAssay(fibroblasts) <- "RNA"

required_metadata <- c(
  "seurat_clusters",
  "orig.ident"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(fibroblasts@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}

cat("\nInput object:\n")
print(fibroblasts)

cat("\nGlobal clusters:\n")
print(table(fibroblasts$seurat_clusters))

# ============================================================
# Preserve global clustering
# ============================================================

fibroblasts$global_cluster <- as.character(
  fibroblasts$seurat_clusters
)

Idents(fibroblasts) <- "seurat_clusters"

# ============================================================
# Select macrophage / monocyte / DC clusters
# ============================================================

fibros.CD45.Macros <- subset(
  fibroblasts,
  subset = seurat_clusters == 15 |
    seurat_clusters == 21
)

if (ncol(fibros.CD45.Macros) == 0) {
  stop("No cells were selected from global clusters 15 and 21.")
}

cat("\nCells selected for macrophage/DC subclustering:\n")
cat(ncol(fibros.CD45.Macros), "\n")

cat("\nCells by global cluster:\n")
print(table(fibros.CD45.Macros$global_cluster))

cat("\nCells by condition:\n")
print(table(fibros.CD45.Macros$orig.ident))

# ============================================================
# Subset preprocessing
# ============================================================
# scale() is defined in R/Functions.R:
#   - ScaleData using all genes
#   - RunPCA using the variable features already present
# ============================================================

fibros.CD45.Macros <- scale(
  fibros.CD45.Macros
)

if (length(VariableFeatures(fibros.CD45.Macros)) == 0) {
  stop("No variable features are available for subclustering.")
}

cat(
  "\nVariable features used: ",
  length(VariableFeatures(fibros.CD45.Macros)),
  "\n",
  sep = ""
)

# ============================================================
# Neighbors and clustering
# ============================================================

fibros.CD45.Macros <- FindNeighbors(
  fibros.CD45.Macros,
  features = VariableFeatures(
    object = fibros.CD45.Macros
  )
)

fibros.CD45.Macros <- FindClusters(
  fibros.CD45.Macros,
  verbose = FALSE,
  resolution = 0.6
)

cat("\nMacrophage/DC subclusters:\n")
print(table(fibros.CD45.Macros$seurat_clusters))

# ============================================================
# PCA diagnostic
# ============================================================

fibros.CD45.Macros.PCA <- PCA(
  fibros.CD45.Macros
)

print(fibros.CD45.Macros.PCA[1])
print(fibros.CD45.Macros.PCA[2])

# ============================================================
# PCA and UMAP
# ============================================================

fibros.CD45.Macros <- RunPCA(
  fibros.CD45.Macros
)

fibros.CD45.Macros <- RunUMAP(
  fibros.CD45.Macros,
  dims = 1:50
)

fibros.CD45.Macros$macro_cluster <- as.character(
  fibros.CD45.Macros$seurat_clusters
)

Idents(fibros.CD45.Macros) <- "seurat_clusters"

# ============================================================
# Diagnostic plot 1: subclusters
# ============================================================

p_macro_clusters <- DimPlot(
  object = fibros.CD45.Macros,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.5
) +
  NoLegend() +
  ggtitle("Macrophage/DC subclusters")

ggsave(
  filename = file.path(
    diagnostics_dir,
    "01_Macrophage_DC_subclusters.pdf"
  ),
  plot = p_macro_clusters,
  width = 8,
  height = 6
)

# ============================================================
# Cluster markers
# ============================================================

DefaultAssay(fibros.CD45.Macros) <- "RNA"
Idents(fibros.CD45.Macros) <- "seurat_clusters"

macro_markers <- FindAllMarkers(
  object = fibros.CD45.Macros,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

macro_top_markers <- macro_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 50,
    with_ties = FALSE
  ) %>%
  ungroup()

write.table(
  macro_markers,
  file = file.path(
    tables_dir,
    "Macrophage_DC_all_markers.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  macro_top_markers,
  file = file.path(
    tables_dir,
    "Macrophage_DC_top50_markers.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Check cluster structure
# ============================================================
# The annotation rules below are defined for clusters 0-10.
# ============================================================

expected_macro_clusters <- as.character(0:10)

observed_macro_clusters <- unique(
  as.character(fibros.CD45.Macros$seurat_clusters)
)

if (!setequal(
  observed_macro_clusters,
  expected_macro_clusters
)) {
  stop(
    paste0(
      "Unexpected macrophage/DC cluster structure.\n",
      "Expected clusters: ",
      paste(expected_macro_clusters, collapse = ", "),
      "\nObserved clusters: ",
      paste(
        sort(as.numeric(observed_macro_clusters)),
        collapse = ", "
      ),
      "\nThe annotation rules cannot be applied safely."
    )
  )
}

# ============================================================
# Diagnostic plot 2: marker panel
# ============================================================

marker_panel <- c(
  # Monocytes and monocyte-derived macrophages
  "Ccr2", "Ly6c2", "Cd14", "Itgam",
  "C1qa", "C1qb", "C1qc", "Adgre1",

  # Alveolar macrophages
  "Siglecf", "Marco", "Pparg", "Car4",

  # mregDC
  "Ccr7", "Fscn1", "Cd274",

  # DC1
  "Xcr1", "Clec9a", "Irf8",

  # DC2
  "Cd209a", "Sirpa", "Irf4",

  # pDC
  "Siglech", "Ly6d", "Tlr9"
)

marker_panel <- intersect(
  marker_panel,
  rownames(fibros.CD45.Macros)
)

p_marker_panel <- FeaturePlot(
  object = fibros.CD45.Macros,
  features = marker_panel,
  reduction = "umap",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q95",
  ncol = 4
)

ggsave(
  filename = file.path(
    diagnostics_dir,
    "02_Macrophage_DC_marker_panel.pdf"
  ),
  plot = p_marker_panel,
  width = 14,
  height = 16,
  limitsize = FALSE
)

# ============================================================
# SCINA annotation
# ============================================================

signature <- read.table(
  signature_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

signature_genes <- unique(
  trimws(
    as.character(
      unlist(signature, use.names = FALSE)
    )
  )
)

signature_genes <- signature_genes[
  !is.na(signature_genes) &
    nzchar(signature_genes)
]

signature_genes_present <- intersect(
  signature_genes,
  rownames(fibros.CD45.Macros)
)

missing_signature_genes <- setdiff(
  signature_genes,
  signature_genes_present
)

cat(
  "\nSCINA signature genes present: ",
  length(signature_genes_present),
  " of ",
  length(signature_genes),
  "\n",
  sep = ""
)

if (length(signature_genes_present) == 0) {
  stop("No SCINA signature genes were found in the object.")
}

if (length(missing_signature_genes) > 0) {
  warning(
    length(missing_signature_genes),
    " SCINA signature genes are absent from the object."
  )
}

scina_matrix <- GetAssayData(
  object = fibros.CD45.Macros,
  assay = "RNA",
  layer = "data"
)

# SCINA only requires genes included in the signature
scina_matrix <- scina_matrix[
  signature_genes_present,
  ,
  drop = FALSE
]

scina_matrix <- as.matrix(scina_matrix)

set.seed(42)

scina_results <- SCINA(
  scina_matrix,
  signature,
  max_iter = 100,
  convergence_n = 10,
  convergence_rate = 0.999,
  sensitivity_cutoff = 0.9,
  rm_overlap = TRUE,
  allow_unknown = TRUE,
  log_file = file.path(
    diagnostics_dir,
    "SCINA.log"
  )
)

scina_labels <- scina_results$cell_labels

if (!is.null(names(scina_labels))) {

  scina_labels <- scina_labels[
    colnames(fibros.CD45.Macros)
  ]

} else {

  if (length(scina_labels) != ncol(fibros.CD45.Macros)) {
    stop(
      "The number of SCINA labels does not match ",
      "the number of cells."
    )
  }

  names(scina_labels) <- colnames(fibros.CD45.Macros)
}

fibros.CD45.Macros$SCINA <- scina_labels[
  colnames(fibros.CD45.Macros)
]

cat("\nSCINA labels:\n")
print(
  table(
    fibros.CD45.Macros$SCINA,
    useNA = "ifany"
  )
)

# ============================================================
# SCINA categories used for annotation
# ============================================================

allowed_scina_labels <- c(
  "Monocytes",
  "Alveolar Macrophages",
  "Inflammatory Monocytes",
  "MoMacs",
  "mregDC",
  "DC1",
  "DC2",
  "pDC"
)

fibros.CD45.Macros$Cluster <- NA_character_

valid_scina <- fibros.CD45.Macros$SCINA %in%
  allowed_scina_labels

fibros.CD45.Macros$Cluster[valid_scina] <-
  fibros.CD45.Macros$SCINA[valid_scina]

# ============================================================
# Record excluded clusters
# ============================================================
# Cluster 5: contaminated cluster
# Cluster 7: plasmablast cluster
# ============================================================

excluded_macros <- data.frame(
  cell = colnames(fibros.CD45.Macros),
  global_cluster = as.character(
    fibros.CD45.Macros$global_cluster
  ),
  macro_cluster = as.character(
    fibros.CD45.Macros$seurat_clusters
  ),
  exclusion_reason = NA_character_,
  stringsAsFactors = FALSE
)

excluded_macros$exclusion_reason[
  excluded_macros$macro_cluster == "5"
] <- "contaminated_cluster"

excluded_macros$exclusion_reason[
  excluded_macros$macro_cluster == "7"
] <- "plasmablast_cluster"

excluded_macros <- excluded_macros[
  !is.na(excluded_macros$exclusion_reason),
  ,
  drop = FALSE
]

excluded_tsv_file <- file.path(
  annotations_dir,
  "barcodes.macros.excluded.tsv"
)

excluded_rds_file <- file.path(
  annotations_dir,
  "barcodes.macros.excluded.rds"
)

write.table(
  excluded_macros,
  file = excluded_tsv_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  excluded_macros,
  file = excluded_rds_file
)

cat("\nExcluded macrophage/DC cells:\n")
print(table(excluded_macros$exclusion_reason))

if (nrow(excluded_macros) > 0) {

  cells_to_keep <- setdiff(
    colnames(fibros.CD45.Macros),
    excluded_macros$cell
  )

  fibros.CD45.Macros <- subset(
    fibros.CD45.Macros,
    cells = cells_to_keep
  )
}

# ============================================================
# Final annotation
# ============================================================

fibros.CD45.Macros$subcluster <- NA_character_

macro_cluster <- as.character(
  fibros.CD45.Macros$seurat_clusters
)

# Monocytes
fibros.CD45.Macros$subcluster[
  macro_cluster %in% c("1", "4")
] <- "Monocytes Ly6cLo"

fibros.CD45.Macros$subcluster[
  macro_cluster == "2"
] <- "Monocytes Ly6cHi"

# Alveolar macrophages
fibros.CD45.Macros$subcluster[
  macro_cluster %in% c("0", "6")
] <- "em-AMs"

# SCINA MoMac assignment
fibros.CD45.Macros$subcluster[
  fibros.CD45.Macros$Cluster == "MoMacs"
] <- "MoMacs"

# Dendritic cells
fibros.CD45.Macros$subcluster[
  macro_cluster == "10"
] <- "mregDC"

fibros.CD45.Macros$subcluster[
  macro_cluster == "9"
] <- "pDC"

fibros.CD45.Macros$subcluster[
  macro_cluster == "3"
] <- "DC2"

fibros.CD45.Macros$subcluster[
  macro_cluster == "3" &
    fibros.CD45.Macros$Cluster == "DC1"
] <- "DC1"

# Remaining macrophages
fibros.CD45.Macros$subcluster[
  is.na(fibros.CD45.Macros$subcluster) &
    macro_cluster %in% c("0", "6", "8")
] <- "MoMacs"

# CCR2KO cluster 6 cells assigned as MoMacs are retained as em-AMs
fibros.CD45.Macros$subcluster[
  macro_cluster == "6" &
    fibros.CD45.Macros$subcluster == "MoMacs" &
    fibros.CD45.Macros$orig.ident == "CCR2KO"
] <- "em-AMs"

# Alveolar macrophages from CD169DTR are assigned as mo-AMs
fibros.CD45.Macros$subcluster[
  fibros.CD45.Macros$orig.ident == "CD169DTR" &
    fibros.CD45.Macros$subcluster == "em-AMs"
] <- "mo-AMs"

# ============================================================
# Check unassigned cells
# ============================================================

unassigned_cells <- is.na(
  fibros.CD45.Macros$subcluster
)

unassigned_macros <- data.frame(
  cell = colnames(fibros.CD45.Macros)[unassigned_cells],
  global_cluster = as.character(
    fibros.CD45.Macros$global_cluster[unassigned_cells]
  ),
  macro_cluster = as.character(
    fibros.CD45.Macros$seurat_clusters[unassigned_cells]
  ),
  SCINA = as.character(
    fibros.CD45.Macros$SCINA[unassigned_cells]
  ),
  stringsAsFactors = FALSE
)

unassigned_file <- file.path(
  annotations_dir,
  "barcodes.macros.unassigned.tsv"
)

if (nrow(unassigned_macros) > 0) {

  write.table(
    unassigned_macros,
    file = unassigned_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  warning(
    nrow(unassigned_macros),
    " macrophage/DC cells remain unassigned. ",
    "See barcodes.macros.unassigned.tsv"
  )

} else if (file.exists(unassigned_file)) {

  unlink(unassigned_file)
}

# ============================================================
# Annotation order and colors
# ============================================================

macro_levels <- c(
  "MoMacs",
  "mo-AMs",
  "em-AMs",
  "DC1",
  "DC2",
  "mregDC",
  "pDC",
  "Monocytes Ly6cLo",
  "Monocytes Ly6cHi"
)

fibros.CD45.Macros$subcluster <- factor(
  fibros.CD45.Macros$subcluster,
  levels = macro_levels
)

Nora_colors <- c(
  "MoMacs" = "#A8D5BA",
  "mo-AMs" = "#66CDAA",
  "em-AMs" = "#2D6A4F",
  "DC1" = "#7F60B6",
  "DC2" = "#D0A0D5",
  "mregDC" = "#FFB84D",
  "pDC" = "#3E8BC4",
  "Monocytes Ly6cLo" = "#FF4D4D",
  "Monocytes Ly6cHi" = "#D21F1F"
)

Idents(fibros.CD45.Macros) <- "subcluster"

cat("\nFinal macrophage/DC annotations:\n")
print(
  table(
    fibros.CD45.Macros$subcluster,
    useNA = "ifany"
  )
)

# ============================================================
# Diagnostic plot 3: annotation overview
# ============================================================

p_scina <- DimPlot(
  object = fibros.CD45.Macros,
  reduction = "umap",
  group.by = "SCINA",
  pt.size = 0.4
) +
  ggtitle("SCINA")

p_final <- DimPlot(
  object = fibros.CD45.Macros,
  reduction = "umap",
  group.by = "subcluster",
  cols = Nora_colors,
  label = TRUE,
  repel = TRUE,
  pt.size = 0.4
) +
  NoLegend() +
  ggtitle("Final annotation")

p_annotation_overview <- plot_grid(
  p_macro_clusters +
    ggtitle("Subclusters before filtering"),
  p_scina,
  p_final,
  ncol = 3
)

ggsave(
  filename = file.path(
    diagnostics_dir,
    "03_Macrophage_DC_annotation_overview.pdf"
  ),
  plot = p_annotation_overview,
  width = 18,
  height = 6
)

# ============================================================
# Save annotated object
# ============================================================

macro_object_file <- file.path(
  objects_dir,
  "fibroblasts.Macrophage_DC.annotated.rds"
)

saveRDS(
  fibros.CD45.Macros,
  file = macro_object_file
)

# ============================================================
# Save barcode annotations
# ============================================================

barcodes.macros <- data.frame(
  cell = colnames(fibros.CD45.Macros),
  global_cluster = as.character(
    fibros.CD45.Macros$global_cluster
  ),
  macro_cluster = as.character(
    fibros.CD45.Macros$seurat_clusters
  ),
  SCINA = as.character(
    fibros.CD45.Macros$SCINA
  ),
  macro_annotation = as.character(
    fibros.CD45.Macros$subcluster
  ),
  stringsAsFactors = FALSE
)

barcodes_rds_file <- file.path(
  annotations_dir,
  "barcodes.macros.rds"
)

barcodes_tsv_file <- file.path(
  annotations_dir,
  "barcodes.macros.tsv"
)

saveRDS(
  barcodes.macros,
  file = barcodes_rds_file
)

write.table(
  barcodes.macros,
  file = barcodes_tsv_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Save annotation summary
# ============================================================

macro_annotation_summary <- barcodes.macros %>%
  count(
    macro_annotation,
    name = "cells"
  ) %>%
  arrange(desc(cells))

summary_file <- file.path(
  tables_dir,
  "Macrophage_DC_annotation_summary.tsv"
)

write.table(
  macro_annotation_summary,
  file = summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Final messages
# ============================================================

cat("\nDone.\n")

cat("\nAnnotated macrophage/DC object:\n")
cat(macro_object_file, "\n")

cat("\nBarcode annotations:\n")
cat(barcodes_tsv_file, "\n")
cat(barcodes_rds_file, "\n")

cat("\nExcluded barcodes:\n")
cat(excluded_tsv_file, "\n")
cat(excluded_rds_file, "\n")

cat("\nAnnotation summary:\n")
cat(summary_file, "\n")

cat("\nDiagnostic plots:\n")
cat(
  file.path(
    diagnostics_dir,
    "01_Macrophage_DC_subclusters.pdf"
  ),
  "\n"
)
cat(
  file.path(
    diagnostics_dir,
    "02_Macrophage_DC_marker_panel.pdf"
  ),
  "\n"
)
cat(
  file.path(
    diagnostics_dir,
    "03_Macrophage_DC_annotation_overview.pdf"
  ),
  "\n"
)