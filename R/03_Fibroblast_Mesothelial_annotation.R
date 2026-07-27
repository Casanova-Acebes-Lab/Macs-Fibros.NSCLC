# ============================================================
# 03_Fibroblast_Mesothelial_annotation.R
# Recover fibroblast and mesothelial annotations
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
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

values <- config$value
names(values) <- config$key

rsdir <- values["rsdir"]
datadir <- values["datadir"]
outdir <- values["outdir"]

# ============================================================
# Directories
# ============================================================

objects_dir <- file.path(rsdir, "objects")
annotations_dir <- file.path(datadir, "annotations")
results_dir <- file.path(outdir, "Fibroblast_Mesothelial")
tables_dir <- file.path(results_dir, "tables")

dir.create(
  objects_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  annotations_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Input objects
# ============================================================

fibroblasts <- readRDS(
  file.path(
    objects_dir,
    "fibroblasts.QC.rds"
  )
)

clusterized <- readRDS(
  file.path(
    annotations_dir,
    "fibroblasts.clusterized.rds"
  )
)

stopifnot(
  inherits(fibroblasts, "Seurat"),
  inherits(clusterized, "Seurat"),
  "Cell_Type2" %in% colnames(clusterized@meta.data),
  "umap" %in% Reductions(fibroblasts)
)

# ============================================================
# Annotation levels and colours
# ============================================================

fibro_meso_levels <- c(
  "Fibroblasts Universal",
  "Fibroblasts Activated 1",
  "Fibroblasts Activated 2",
  "Mesothelial"
)

fibro_levels <- c(
  "Fibroblasts Universal",
  "Fibroblasts Activated 1",
  "Fibroblasts Activated 2"
)

Nora <- c(
  "Mesothelial" = "#FFD9B3",
  "Fibroblasts Universal" = "#F4A825",
  "Fibroblasts Activated 1" = "#B71C1C",
  "Fibroblasts Activated 2" = "#6D1A1A"
)

expected_counts <- c(
  "Fibroblasts Universal" = 240,
  "Fibroblasts Activated 1" = 75,
  "Fibroblasts Activated 2" = 24,
  "Mesothelial" = 372
)

# ============================================================
# Recover annotations by barcode
# ============================================================

old_meta <- clusterized@meta.data

recovered <- data.frame(
  cell = rownames(old_meta),
  fibro_meso_annotation = as.character(
    old_meta$Cell_Type2
  ),
  stringsAsFactors = FALSE
)

recovered <- recovered[
  recovered$fibro_meso_annotation %in% fibro_meso_levels,
  ,
  drop = FALSE
]

stopifnot(
  nrow(recovered) == 711,
  !anyDuplicated(recovered$cell),
  all(recovered$cell %in% colnames(fibroblasts))
)

observed_counts <- table(
  recovered$fibro_meso_annotation
)

stopifnot(
  identical(
    as.integer(
      observed_counts[names(expected_counts)]
    ),
    as.integer(expected_counts)
  )
)

# ============================================================
# Add metadata and UMAP coordinates
# ============================================================

qc_meta <- fibroblasts@meta.data

recovered$global_cluster <- as.character(
  qc_meta[
    recovered$cell,
    "seurat_clusters"
  ]
)

recovered$condition <- as.character(
  qc_meta[
    recovered$cell,
    "orig.ident"
  ]
)

recovered$tag <- as.character(
  qc_meta[
    recovered$cell,
    "tag"
  ]
)

global_umap <- Embeddings(
  fibroblasts,
  reduction = "umap"
)

recovered$UMAP_1 <- global_umap[
  recovered$cell,
  1
]

recovered$UMAP_2 <- global_umap[
  recovered$cell,
  2
]

cat(
  "\nRecovered annotations before filtering:\n"
)

print(
  table(
    recovered$fibro_meso_annotation
  )
)

# ============================================================
# Remove isolated fibroblast cells
# ============================================================
# Outliers are identified using the distance to the tenth
# nearest fibroblast neighbour in the global UMAP.
# Mesothelial cells are not included in this filtering.
# ============================================================

fibro_rows <- recovered$fibro_meso_annotation %in%
  fibro_levels

fibro_umap <- as.matrix(
  recovered[
    fibro_rows,
    c(
      "UMAP_1",
      "UMAP_2"
    ),
    drop = FALSE
  ]
)

rownames(fibro_umap) <- recovered$cell[
  fibro_rows
]

fibro_distances <- as.matrix(
  stats::dist(
    fibro_umap
  )
)

diag(fibro_distances) <- Inf

k_neighbor <- 10L

knn_distance <- apply(
  fibro_distances,
  1,
  function(x) {
    sort(
      x,
      partial = k_neighbor
    )[k_neighbor]
  }
)

median_knn <- median(
  knn_distance,
  na.rm = TRUE
)

mad_knn <- stats::mad(
  knn_distance,
  na.rm = TRUE
)

mad_multiplier <- 5

if (
  !is.finite(mad_knn) ||
    mad_knn == 0
) {

  outlier_cutoff <- as.numeric(
    quantile(
      knn_distance,
      probs = 0.99,
      na.rm = TRUE
    )
  )

} else {

  outlier_cutoff <- median_knn +
    mad_multiplier * mad_knn
}

fibro_outlier_cells <- names(
  knn_distance[
    knn_distance > outlier_cutoff
  ]
)

barcodes.fibro_meso.excluded <- recovered[
  recovered$cell %in% fibro_outlier_cells,
  ,
  drop = FALSE
]

barcodes.fibro_meso.excluded$knn_distance <-
  knn_distance[
    barcodes.fibro_meso.excluded$cell
  ]

barcodes.fibro_meso.excluded$exclusion_reason <-
  "isolated_fibroblast_UMAP_cell"

barcodes.fibro_meso <- recovered[
  !recovered$cell %in% fibro_outlier_cells,
  ,
  drop = FALSE
]

stopifnot(
  !any(
    barcodes.fibro_meso.excluded$
      fibro_meso_annotation == "Mesothelial"
  ),
  !anyDuplicated(
    barcodes.fibro_meso$cell
  ),
  nrow(barcodes.fibro_meso) +
    nrow(barcodes.fibro_meso.excluded) ==
    nrow(recovered)
)

cat(
  "\nFibroblast UMAP filtering:\n"
)

cat(
  "k-neighbour:",
  k_neighbor,
  "\n"
)

cat(
  "Median kNN distance:",
  median_knn,
  "\n"
)

cat(
  "MAD:",
  mad_knn,
  "\n"
)

cat(
  "Outlier cutoff:",
  outlier_cutoff,
  "\n"
)

cat(
  "Excluded fibroblasts:",
  length(fibro_outlier_cells),
  "\n"
)

cat(
  "\nExcluded cells by annotation:\n"
)

print(
  table(
    barcodes.fibro_meso.excluded$
      fibro_meso_annotation
  )
)

cat(
  "\nFinal retained annotations:\n"
)

print(
  table(
    barcodes.fibro_meso$
      fibro_meso_annotation
  )
)

# ============================================================
# Save barcode annotations
# ============================================================

barcodes.fibro_meso$fibro_meso_annotation <- factor(
  barcodes.fibro_meso$fibro_meso_annotation,
  levels = fibro_meso_levels
)

saveRDS(
  barcodes.fibro_meso,
  file.path(
    annotations_dir,
    "barcodes.fibro_meso.rds"
  )
)

write.table(
  barcodes.fibro_meso,
  file.path(
    annotations_dir,
    "barcodes.fibro_meso.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  barcodes.fibro_meso.excluded,
  file.path(
    annotations_dir,
    "barcodes.fibro_meso.excluded.rds"
  )
)

write.table(
  barcodes.fibro_meso.excluded,
  file.path(
    annotations_dir,
    "barcodes.fibro_meso.excluded.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Create fibroblast/mesothelial object
# ============================================================

fibro_meso <- subset(
  fibroblasts,
  cells = barcodes.fibro_meso$cell
)

annotation_match <- match(
  colnames(fibro_meso),
  barcodes.fibro_meso$cell
)

fibro_meso$fibro_meso_annotation <-
  barcodes.fibro_meso$
    fibro_meso_annotation[
      annotation_match
    ]

fibro_meso$fibro_meso_annotation <- factor(
  fibro_meso$fibro_meso_annotation,
  levels = fibro_meso_levels
)

stopifnot(
  !anyNA(
    fibro_meso$fibro_meso_annotation
  )
)

Idents(fibro_meso) <- "fibro_meso_annotation"

fibro_meso_file <- file.path(
  objects_dir,
  "fibroblasts.Fibro_Mesothelial.annotated.rds"
)

saveRDS(
  fibro_meso,
  fibro_meso_file
)

# ============================================================
# Create fibroblast-only object
# ============================================================

fibro_only <- subset(
  fibro_meso,
  idents = fibro_levels
)

fibro_only$fibroblast_annotation <- factor(
  as.character(
    fibro_only$fibro_meso_annotation
  ),
  levels = fibro_levels
)

Idents(fibro_only) <- "fibroblast_annotation"

fibro_only_file <- file.path(
  objects_dir,
  "fibroblasts.Fibroblasts.annotated.rds"
)

saveRDS(
  fibro_only,
  fibro_only_file
)

# ============================================================
# Summary tables
# ============================================================

fibro_meso_summary <- barcodes.fibro_meso %>%
  count(
    fibro_meso_annotation,
    name = "cells"
  ) %>%
  mutate(
    percentage = 100 * cells / sum(cells)
  )

fibro_meso_by_condition <- barcodes.fibro_meso %>%
  count(
    condition,
    fibro_meso_annotation,
    name = "cells"
  ) %>%
  group_by(condition) %>%
  mutate(
    total_fibro_meso = sum(cells),
    percentage = 100 * cells / total_fibro_meso
  ) %>%
  ungroup()

fibro_meso_by_sample <- barcodes.fibro_meso %>%
  count(
    condition,
    tag,
    fibro_meso_annotation,
    name = "cells"
  )

fibro_only_by_condition <- barcodes.fibro_meso %>%
  filter(
    fibro_meso_annotation %in% fibro_levels
  ) %>%
  count(
    condition,
    fibro_meso_annotation,
    name = "cells"
  ) %>%
  group_by(condition) %>%
  mutate(
    total_fibroblasts = sum(cells),
    percentage = 100 * cells / total_fibroblasts
  ) %>%
  ungroup()

write.table(
  fibro_meso_summary,
  file.path(
    tables_dir,
    "Fibroblast_Mesothelial_annotation_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  fibro_meso_by_condition,
  file.path(
    tables_dir,
    "Fibroblast_Mesothelial_counts_by_condition.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  fibro_meso_by_sample,
  file.path(
    tables_dir,
    "Fibroblast_Mesothelial_counts_by_sample.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  fibro_only_by_condition,
  file.path(
    tables_dir,
    "Fibroblast_counts_by_condition.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# UMAP: fibroblast and mesothelial annotations
# ============================================================

p_fibro_meso <- DimPlot(
  object = fibro_meso,
  reduction = "umap",
  group.by = "fibro_meso_annotation",
  cols = Nora,
  label = TRUE,
  repel = TRUE,
  pt.size = 0.7
) +
  NoLegend() +
  ggtitle(
    "Fibroblast and mesothelial annotations"
  )

ggsave(
  filename = file.path(
    results_dir,
    "Fibroblast_Mesothelial_UMAP.pdf"
  ),
  plot = p_fibro_meso,
  width = 8,
  height = 6
)

# ============================================================
# UMAP split by condition
# ============================================================

p_fibro_meso_split <- DimPlot(
  object = fibro_meso,
  reduction = "umap",
  group.by = "fibro_meso_annotation",
  split.by = "orig.ident",
  cols = Nora,
  pt.size = 0.6,
  ncol = 3
)

ggsave(
  filename = file.path(
    results_dir,
    "Fibroblast_Mesothelial_UMAP_by_condition.pdf"
  ),
  plot = p_fibro_meso_split,
  width = 18,
  height = 6
)

# ============================================================
# UMAP: fibroblast states
# ============================================================

p_fibro_only <- DimPlot(
  object = fibro_only,
  reduction = "umap",
  group.by = "fibroblast_annotation",
  cols = Nora[fibro_levels],
  label = TRUE,
  repel = TRUE,
  pt.size = 0.8
) +
  NoLegend() +
  ggtitle(
    "Fibroblast states"
  )

ggsave(
  filename = file.path(
    results_dir,
    "Fibroblast_states_UMAP.pdf"
  ),
  plot = p_fibro_only,
  width = 8,
  height = 6
)

# ============================================================
# Marker panel
# ============================================================

marker_panel <- c(
  "Col1a1",
  "Col1a2",
  "Pdgfra",
  "Clec3b",
  "Pi16",
  "Dpep1",
  "Hhip",
  "Grem2",
  "Bmp5",
  "Ddit4",
  "Nampt",
  "Npr3",
  "Msln",
  "Wt1",
  "Krt8",
  "Krt18"
)

marker_panel <- intersect(
  marker_panel,
  rownames(fibro_meso)
)

p_marker_panel <- FeaturePlot(
  object = fibro_meso,
  features = marker_panel,
  reduction = "umap",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q95",
  ncol = 4
)

ggsave(
  filename = file.path(
    results_dir,
    "Fibroblast_Mesothelial_marker_panel.pdf"
  ),
  plot = p_marker_panel,
  width = 14,
  height = 12,
  limitsize = FALSE
)

# ============================================================
# Filtering diagnostics
# ============================================================

filter_summary <- data.frame(
  parameter = c(
    "k_neighbor",
    "median_knn_distance",
    "mad_knn_distance",
    "mad_multiplier",
    "outlier_cutoff",
    "recovered_fibroblasts",
    "excluded_fibroblasts",
    "retained_fibroblasts",
    "retained_mesothelial"
  ),
  value = c(
    k_neighbor,
    median_knn,
    mad_knn,
    mad_multiplier,
    outlier_cutoff,
    sum(fibro_rows),
    nrow(barcodes.fibro_meso.excluded),
    ncol(fibro_only),
    sum(
      fibro_meso$fibro_meso_annotation ==
        "Mesothelial"
    )
  )
)

write.table(
  filter_summary,
  file.path(
    tables_dir,
    "Fibroblast_UMAP_filter_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Final report
# ============================================================

cat(
  "\nDone.\n"
)

cat(
  "\nFibroblast/mesothelial object:\n"
)

cat(
  fibro_meso_file,
  "\n"
)

cat(
  "\nFibroblast-only object:\n"
)

cat(
  fibro_only_file,
  "\n"
)

cat(
  "\nFinal fibroblast/mesothelial counts:\n"
)

print(
  table(
    fibro_meso$fibro_meso_annotation
  )
)

cat(
  "\nExcluded isolated fibroblasts:\n"
)

print(
  table(
    barcodes.fibro_meso.excluded$
      fibro_meso_annotation
  )
)