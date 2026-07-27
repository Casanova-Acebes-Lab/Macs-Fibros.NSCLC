# ============================================================
# 04_CD45neg.R
# CD45-negative annotation
#
# Input objects:
#   results/objects/fibroblasts.QC.rds
#   results/objects/fibroblasts.Fibro_Mesothelial.annotated.rds
#
# Output object:
#   results/objects/fibroblasts.CD45NP.annotated.rds
#
# Annotation columns:
#
#   Clusters
#     Broad CD45-negative annotation:
#       Endothelial Vascular
#       Aerocytes
#       Epithelial AT1
#       Epithelial AT2
#       Epithelial Krt8+
#       Ciliated
#       Mesenchymal
#
#   fibroblast_annotation
#     Same annotation, replacing Mesenchymal with:
#       Fibroblasts Universal
#       Fibroblasts Activated 1
#       Fibroblasts Activated 2
#       Mesothelial
#
#   cancer_annotation
#     Same annotation, replacing tumour epithelial cells:
#       cluster 12                  -> Prolif. Cancer Cells
#       remaining Epithelial Krt8+ -> EMT Cancer Cells
#
# Epithelial Krt8+ cells are removed from naive samples.
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
# Configuration
# ============================================================

config <- read.delim(
  "config.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("key", "value") %in% colnames(config))) {
  stop(
    "config.tsv must contain the columns 'key' and 'value'."
  )
}

values <- setNames(
  config$value,
  config$key
)

required_config <- c(
  "rsdir",
  "datadir",
  "outdir"
)

missing_config <- setdiff(
  required_config,
  names(values)
)

if (length(missing_config) > 0) {
  stop(
    "Missing keys in config.tsv: ",
    paste(
      missing_config,
      collapse = ", "
    )
  )
}

rsdir <- unname(
  values[["rsdir"]]
)

datadir <- unname(
  values[["datadir"]]
)

outdir <- unname(
  values[["outdir"]]
)

objects_dir <- file.path(
  rsdir,
  "objects"
)

cd45neg_dir <- file.path(
  outdir,
  "CD45neg"
)

dir.create(
  objects_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  cd45neg_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

source("R/Functions.R")

# ============================================================
# Input and output files
# ============================================================

qc_file <- file.path(
  objects_dir,
  "fibroblasts.QC.rds"
)

fibro_meso_file <- file.path(
  objects_dir,
  "fibroblasts.Fibro_Mesothelial.annotated.rds"
)

output_file <- file.path(
  objects_dir,
  "fibroblasts.CD45NP.annotated.rds"
)

input_files <- c(
  qc_file,
  fibro_meso_file
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing input files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

# ============================================================
# Read objects
# ============================================================

fibroblasts.QC <- readRDS(
  qc_file
)

fibroblasts.Fibro_Mesothelial <- readRDS(
  fibro_meso_file
)

stopifnot(
  inherits(
    fibroblasts.QC,
    "Seurat"
  ),
  inherits(
    fibroblasts.Fibro_Mesothelial,
    "Seurat"
  )
)

DefaultAssay(
  fibroblasts.QC
) <- "RNA"

cat(
  "\nLoaded QC object:\n",
  qc_file,
  "\n",
  sep = ""
)

cat(
  "\nLoaded fibroblast/mesothelial object:\n",
  fibro_meso_file,
  "\n",
  sep = ""
)

# ============================================================
# Check required metadata
# ============================================================

required_qc_metadata <- c(
  "CD",
  "seurat_clusters",
  "group",
  "orig.ident"
)

missing_qc_metadata <- setdiff(
  required_qc_metadata,
  colnames(
    fibroblasts.QC@meta.data
  )
)

if (length(missing_qc_metadata) > 0) {
  stop(
    "Missing metadata in fibroblasts.QC: ",
    paste(
      missing_qc_metadata,
      collapse = ", "
    )
  )
}

if (
  !"fibro_meso_annotation" %in%
    colnames(
      fibroblasts.Fibro_Mesothelial@meta.data
    )
) {
  stop(
    paste0(
      "fibroblasts.Fibro_Mesothelial.annotated.rds ",
      "does not contain 'fibro_meso_annotation'."
    )
  )
}

if (
  !"umap" %in%
    Reductions(fibroblasts.QC)
) {
  stop(
    "The UMAP reduction is missing from fibroblasts.QC."
  )
}

# ============================================================
# Check fibroblast/mesothelial annotations
# ============================================================

expected_mesenchymal_labels <- c(
  "Fibroblasts Universal",
  "Fibroblasts Activated 1",
  "Fibroblasts Activated 2",
  "Mesothelial"
)

observed_mesenchymal_labels <- unique(
  as.character(
    fibroblasts.Fibro_Mesothelial$
      fibro_meso_annotation
  )
)

unexpected_mesenchymal_labels <- setdiff(
  observed_mesenchymal_labels[
    !is.na(observed_mesenchymal_labels)
  ],
  expected_mesenchymal_labels
)

if (length(unexpected_mesenchymal_labels) > 0) {
  stop(
    "Unexpected fibroblast/mesothelial labels: ",
    paste(
      unexpected_mesenchymal_labels,
      collapse = ", "
    )
  )
}

cat(
  "\nInput fibroblast/mesothelial annotation:\n"
)

print(
  table(
    fibroblasts.Fibro_Mesothelial$
      fibro_meso_annotation,
    useNA = "ifany"
  )
)

# ============================================================
# Select CD45-negative cells
# ============================================================

fibros.CD45NP <- subset(
  fibroblasts.QC,
  subset = CD == "CD45+",
  invert = TRUE
)

cat(
  "\nCD45-negative cells before filtering: ",
  ncol(fibros.CD45NP),
  "\n",
  sep = ""
)

# ============================================================
# Remove unwanted global clusters
# ============================================================

excluded_clusters <- c(
  0,
  1,
  4,
  5,
  6,
  7,
  8,
  10,
  14,
  15,
  19,
  21,
  22,
  24,
  25
)

global_clusters <- as.character(
  fibros.CD45NP$seurat_clusters
)

cells_keep <- colnames(
  fibros.CD45NP
)[
  !global_clusters %in%
    as.character(excluded_clusters)
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

cells_keep <- rownames(
  umap_coordinates
)[
  umap_coordinates[, 1] <= 7
]

fibros.CD45NP <- subset(
  fibros.CD45NP,
  cells = cells_keep
)

cat(
  "\nCD45-negative cells after initial filtering: ",
  ncol(fibros.CD45NP),
  "\n",
  sep = ""
)

cat(
  "\nGlobal clusters after initial filtering:\n"
)

print(
  table(
    fibros.CD45NP$seurat_clusters,
    useNA = "ifany"
  )
)

# ============================================================
# Prepare normalized RNA matrix for SCINA
# ============================================================

DefaultAssay(
  fibros.CD45NP
) <- "RNA"

if (
  inherits(
    fibros.CD45NP[["RNA"]],
    "Assay5"
  )
) {
  fibros.CD45NP <- JoinLayers(
    fibros.CD45NP,
    assay = "RNA"
  )
}

expression_matrix <- GetAssayData(
  object = fibros.CD45NP,
  assay = "RNA",
  layer = "data"
)

expression_matrix <- as.matrix(
  expression_matrix
)

# ============================================================
# Read SCINA signatures
# ============================================================

signature_file <- file.path(
  datadir,
  "signatures",
  "signature.endothelial.tsv"
)

if (!file.exists(signature_file)) {
  stop(
    "SCINA signature file not found: ",
    signature_file
  )
}

signature <- read.table(
  signature_file,
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
    cd45neg_dir,
    "SCINA.log"
  )
)

fibros.CD45NP$SCINA <- as.character(
  scina_results$cell_labels
)

cat(
  "\nSCINA annotation:\n"
)

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

cat(
  "\nCells after removing platelet cluster 23: ",
  ncol(fibros.CD45NP),
  "\n",
  sep = ""
)

# ============================================================
# Recover fibroblast and mesothelial labels
# ============================================================

fibro_labels <- setNames(
  as.character(
    fibroblasts.Fibro_Mesothelial$
      fibro_meso_annotation
  ),
  colnames(
    fibroblasts.Fibro_Mesothelial
  )
)

valid_fibro_cells <- names(
  fibro_labels
)[
  !is.na(fibro_labels) &
    fibro_labels %in%
      expected_mesenchymal_labels
]

cat(
  "\nCells with fine fibroblast/mesothelial annotation: ",
  length(valid_fibro_cells),
  "\n",
  sep = ""
)

print(
  table(
    fibro_labels[
      valid_fibro_cells
    ],
    useNA = "ifany"
  )
)

# ============================================================
# Remove cells discarded during fibroblast reclustering
#
# Global clusters 17 and 20 correspond to the mesenchymal
# compartment. Only cells retained in script 03 are kept.
# ============================================================

global_mesenchymal_cells <- colnames(
  fibros.CD45NP
)[
  as.character(
    fibros.CD45NP$seurat_clusters
  ) %in% c(
    "17",
    "20"
  )
]

common_mesenchymal_cells <- intersect(
  global_mesenchymal_cells,
  valid_fibro_cells
)

discarded_mesenchymal_cells <- setdiff(
  global_mesenchymal_cells,
  valid_fibro_cells
)

cat(
  "\nMesenchymal cells in global clusters 17/20: ",
  length(global_mesenchymal_cells),
  "\n",
  sep = ""
)

cat(
  "Mesenchymal cells retained from script 03: ",
  length(common_mesenchymal_cells),
  "\n",
  sep = ""
)

cat(
  "Mesenchymal cells removed: ",
  length(discarded_mesenchymal_cells),
  "\n",
  sep = ""
)

fibros.CD45NP <- subset(
  fibros.CD45NP,
  cells = setdiff(
    colnames(fibros.CD45NP),
    discarded_mesenchymal_cells
  )
)

# ============================================================
# Broad CD45-negative annotation
# ============================================================

fibros.CD45NP$Clusters <- NA_character_

global_clusters <- as.character(
  fibros.CD45NP$seurat_clusters
)

# Endothelial vascular
fibros.CD45NP$Clusters[
  global_clusters %in% c(
    "2",
    "3",
    "18"
  )
] <- "Endothelial Vascular"

# Aerocytes
fibros.CD45NP$Clusters[
  global_clusters %in% c(
    "9",
    "11"
  )
] <- "Aerocytes"

# Mesenchymal compartment
fibros.CD45NP$Clusters[
  global_clusters %in% c(
    "17",
    "20"
  )
] <- "Mesenchymal"

# Krt8-positive epithelial compartment
fibros.CD45NP$Clusters[
  global_clusters %in% c(
    "12",
    "13"
  )
] <- "Epithelial Krt8+"

# Cluster 16 resolved using SCINA
fibros.CD45NP$Clusters[
  global_clusters == "16" &
    fibros.CD45NP$SCINA ==
      "Epithelial AT1"
] <- "Epithelial AT1"

fibros.CD45NP$Clusters[
  global_clusters == "16" &
    fibros.CD45NP$SCINA ==
      "Epithelial AT2"
] <- "Epithelial AT2"

# Unknown cells in cluster 16 were originally assigned as AT2
fibros.CD45NP$Clusters[
  global_clusters == "16" &
    fibros.CD45NP$SCINA ==
      "unknown"
] <- "Epithelial AT2"

fibros.CD45NP$Clusters[
  global_clusters == "16" &
    fibros.CD45NP$SCINA ==
      "Ciliated"
] <- "Ciliated"

cat(
  "\nBroad annotation before removing unresolved cells:\n"
)

print(
  table(
    fibros.CD45NP$Clusters,
    useNA = "ifany"
  )
)

# ============================================================
# Keep only cells with a broad annotation
# ============================================================

annotated_cells <- colnames(
  fibros.CD45NP
)[
  !is.na(
    fibros.CD45NP$Clusters
  )
]

unresolved_cells <- setdiff(
  colnames(fibros.CD45NP),
  annotated_cells
)

cat(
  "\nCells removed because they lack a broad annotation: ",
  length(unresolved_cells),
  "\n",
  sep = ""
)

fibros.CD45NP <- subset(
  fibros.CD45NP,
  cells = annotated_cells
)

# ============================================================
# Fine annotation with fibroblast and mesothelial identities
#
# Start from Clusters and replace Mesenchymal with:
#   Fibroblasts Universal
#   Fibroblasts Activated 1
#   Fibroblasts Activated 2
#   Mesothelial
# ============================================================

fibros.CD45NP$fibroblast_annotation <- as.character(
  fibros.CD45NP$Clusters
)

mesenchymal_cells <- colnames(
  fibros.CD45NP
)[
  as.character(
    fibros.CD45NP$Clusters
  ) == "Mesenchymal"
]

common_fibro_cells <- intersect(
  mesenchymal_cells,
  valid_fibro_cells
)

fibros.CD45NP$fibroblast_annotation[
  match(
    common_fibro_cells,
    colnames(fibros.CD45NP)
  )
] <- unname(
  fibro_labels[
    common_fibro_cells
  ]
)

cat(
  "\nMesenchymal cells with fine annotation: ",
  length(common_fibro_cells),
  " of ",
  length(mesenchymal_cells),
  "\n",
  sep = ""
)

remaining_mesenchymal <- sum(
  fibros.CD45NP$fibroblast_annotation ==
    "Mesenchymal",
  na.rm = TRUE
)

if (remaining_mesenchymal > 0) {
  stop(
    remaining_mesenchymal,
    " mesenchymal cells lack a fine annotation."
  )
}

cat(
  "\nAnnotation with fibroblast subtypes:\n"
)

print(
  table(
    fibros.CD45NP$fibroblast_annotation,
    useNA = "ifany"
  )
)

# ============================================================
# Remove Epithelial Krt8+ cells from naive samples
#
# These cells correspond to cancer-cell populations in tumour
# samples and must not remain in naive samples.
# ============================================================

naive_krt8_cells <- colnames(
  fibros.CD45NP
)[
  as.character(
    fibros.CD45NP$group
  ) == "Naive" &
    as.character(
      fibros.CD45NP$fibroblast_annotation
    ) == "Epithelial Krt8+"
]

cat(
  "\nEpithelial Krt8+ cells removed from naive samples: ",
  length(naive_krt8_cells),
  "\n",
  sep = ""
)

fibros.CD45NP <- subset(
  fibros.CD45NP,
  cells = setdiff(
    colnames(fibros.CD45NP),
    naive_krt8_cells
  )
)

# ============================================================
# Final annotation including cancer cells
#
# Cluster 12:
#   Prolif. Cancer Cells
#
# Remaining Epithelial Krt8+:
#   EMT Cancer Cells
# ============================================================

fibros.CD45NP$cancer_annotation <- as.character(
  fibros.CD45NP$fibroblast_annotation
)

global_clusters <- as.character(
  fibros.CD45NP$seurat_clusters
)

fibros.CD45NP$cancer_annotation[
  global_clusters == "12"
] <- "Prolif. Cancer Cells"

fibros.CD45NP$cancer_annotation[
  fibros.CD45NP$cancer_annotation ==
    "Epithelial Krt8+"
] <- "EMT Cancer Cells"

# ============================================================
# Factor levels
# ============================================================

cluster_levels <- c(
  "Endothelial Vascular",
  "Aerocytes",
  "Epithelial AT1",
  "Epithelial AT2",
  "Epithelial Krt8+",
  "Ciliated",
  "Mesenchymal"
)

fibroblast_annotation_levels <- c(
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

cancer_annotation_levels <- c(
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

fibros.CD45NP$Clusters <- factor(
  as.character(
    fibros.CD45NP$Clusters
  ),
  levels = cluster_levels
)

fibros.CD45NP$fibroblast_annotation <- factor(
  as.character(
    fibros.CD45NP$fibroblast_annotation
  ),
  levels = fibroblast_annotation_levels
)

fibros.CD45NP$cancer_annotation <- factor(
  as.character(
    fibros.CD45NP$cancer_annotation
  ),
  levels = cancer_annotation_levels
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

if (anyNA(fibros.CD45NP$Clusters)) {
  stop(
    "NA values remain in Clusters."
  )
}

if (
  anyNA(
    fibros.CD45NP$fibroblast_annotation
  )
) {
  stop(
    "NA values remain in fibroblast_annotation."
  )
}

if (
  anyNA(
    fibros.CD45NP$cancer_annotation
  )
) {
  stop(
    "NA values remain in cancer_annotation."
  )
}

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

cat(
  "\nBroad CD45-negative annotation:\n"
)

print(
  table(
    fibros.CD45NP$Clusters,
    useNA = "ifany"
  )
)

cat(
  "\nAnnotation with fibroblast subtypes:\n"
)

print(
  table(
    fibros.CD45NP$fibroblast_annotation,
    useNA = "ifany"
  )
)

cat(
  "\nFinal annotation including cancer cells:\n"
)

print(
  table(
    fibros.CD45NP$cancer_annotation,
    useNA = "ifany"
  )
)

cat(
  "\nFinal annotation by experimental group:\n"
)

print(
  table(
    group = fibros.CD45NP$group,
    annotation = fibros.CD45NP$cancer_annotation,
    useNA = "ifany"
  )
)

cat(
  "\nFinal annotation by condition:\n"
)

print(
  table(
    condition = fibros.CD45NP$orig.ident,
    annotation = fibros.CD45NP$cancer_annotation,
    useNA = "ifany"
  )
)

cat(
  "\nFinal annotation by condition and group:\n"
)

print(
  table(
    condition = fibros.CD45NP$orig.ident,
    group = fibros.CD45NP$group,
    annotation = fibros.CD45NP$cancer_annotation,
    useNA = "ifany"
  )
)

# ============================================================
# Final UMAP plots
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
    "Broad CD45-negative annotation"
  )

p_fibro <- DimPlot_scCustom(
  fibros.CD45NP,
  reduction = "umap",
  group.by = "fibroblast_annotation",
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

print(
  p_final
)

ggsave(
  filename = file.path(
    cd45neg_dir,
    "CD45neg_final_annotation_UMAP.pdf"
  ),
  plot = p_final,
  width = 18,
  height = 18,
  limitsize = FALSE
)

# ============================================================
# Set final identity
# ============================================================

Idents(
  fibros.CD45NP
) <- "cancer_annotation"

# ============================================================
# Save final object
# ============================================================

saveRDS(
  fibros.CD45NP,
  output_file
)

cat(
  "\nFinal CD45-negative object saved in:\n",
  output_file,
  "\n",
  sep = ""
)