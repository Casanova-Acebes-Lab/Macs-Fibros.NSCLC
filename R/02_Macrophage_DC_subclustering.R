# ============================================================
# 02_Macrophage_DC_annotation.R
# Recover macrophage, monocyte and DC annotations
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(igraph)
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

objects_dir <- file.path(
  rsdir,
  "objects"
)

annotations_dir <- file.path(
  datadir,
  "annotations"
)

results_dir <- file.path(
  outdir,
  "Macrophage_DC"
)

tables_dir <- file.path(
  results_dir,
  "tables"
)

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
# Input files
# ============================================================

qc_file <- file.path(
  objects_dir,
  "fibroblasts.QC.rds"
)

annotation_object_file <- file.path(
  annotations_dir,
  "fibroblasts.clusterized.rds"
)

stopifnot(
  file.exists(qc_file),
  file.exists(annotation_object_file)
)

# ============================================================
# Load objects
# ============================================================

fibroblasts <- readRDS(
  qc_file
)

clusterized <- readRDS(
  annotation_object_file
)

stopifnot(
  inherits(fibroblasts, "Seurat"),
  inherits(clusterized, "Seurat"),
  "Cell_Type2" %in% colnames(clusterized@meta.data),
  "umap" %in% Reductions(fibroblasts)
)

# ============================================================
# Macrophage/DC annotations
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

expected_counts <- c(
  "Monocytes Ly6cLo" = 318,
  "Monocytes Ly6cHi" = 162,
  "MoMacs" = 68,
  "em-AMs" = 142,
  "mo-AMs" = 161,
  "mregDC" = 23,
  "DC1" = 31,
  "DC2" = 98,
  "pDC" = 47
)

# ============================================================
# Recover annotations by barcode
# ============================================================

old_meta <- clusterized@meta.data

recovered_macros <- data.frame(
  cell = rownames(old_meta),
  macro_annotation = as.character(
    old_meta$Cell_Type2
  ),
  stringsAsFactors = FALSE
)

recovered_macros <- recovered_macros[
  recovered_macros$macro_annotation %in% macro_levels,
  ,
  drop = FALSE
]

stopifnot(
  nrow(recovered_macros) == 1050,
  !anyDuplicated(recovered_macros$cell),
  all(recovered_macros$cell %in% colnames(fibroblasts))
)

observed_counts <- table(
  recovered_macros$macro_annotation
)

stopifnot(
  identical(
    as.integer(
      observed_counts[names(expected_counts)]
    ),
    as.integer(expected_counts)
  )
)

cat("\nRecovered annotations before UMAP filtering:\n")

print(
  table(
    recovered_macros$macro_annotation
  )
)

# ============================================================
# Add metadata from the reconstructed QC object
# ============================================================

qc_meta <- fibroblasts@meta.data

recovered_macros$global_cluster <- as.character(
  qc_meta[
    recovered_macros$cell,
    "seurat_clusters"
  ]
)

recovered_macros$condition <- as.character(
  qc_meta[
    recovered_macros$cell,
    "orig.ident"
  ]
)

recovered_macros$tag <- as.character(
  qc_meta[
    recovered_macros$cell,
    "tag"
  ]
)

# ============================================================
# Global UMAP coordinates
# ============================================================

global_umap <- Embeddings(
  fibroblasts,
  reduction = "umap"
)

umap_macros <- global_umap[
  recovered_macros$cell,
  1:2,
  drop = FALSE
]

colnames(umap_macros) <- c(
  "UMAP_1",
  "UMAP_2"
)

recovered_macros$UMAP_1 <- umap_macros[
  recovered_macros$cell,
  "UMAP_1"
]

recovered_macros$UMAP_2 <- umap_macros[
  recovered_macros$cell,
  "UMAP_2"
]

# ============================================================
# Keep the main macrophage/DC UMAP component
# ============================================================
# A radius graph is generated from the global UMAP coordinates.
# The radius is derived from local 10-nearest-neighbour distances.
# Only the largest connected component is retained.
# ============================================================

umap_distance <- as.matrix(
  stats::dist(
    umap_macros
  )
)

diag(umap_distance) <- Inf

k_local <- 10L

kth_neighbor_distance <- apply(
  umap_distance,
  1,
  function(x) {
    sort(
      x,
      partial = k_local
    )[k_local]
  }
)

radius_quantile <- 0.95
radius_multiplier <- 1.5

connection_radius <- as.numeric(
  stats::quantile(
    kth_neighbor_distance,
    probs = radius_quantile,
    na.rm = TRUE
  )
) * radius_multiplier

umap_adjacency <- umap_distance <= connection_radius

diag(umap_adjacency) <- FALSE

umap_graph <- igraph::graph_from_adjacency_matrix(
  umap_adjacency,
  mode = "undirected",
  diag = FALSE
)

umap_components <- igraph::components(
  umap_graph
)

main_component <- which.max(
  umap_components$csize
)

keep_main_component <- (
  umap_components$membership == main_component
)

names(keep_main_component) <- rownames(
  umap_macros
)

# ============================================================
# Cells excluded from the main UMAP component
# ============================================================

barcodes.macros.umap_excluded <- recovered_macros[
  !keep_main_component[
    recovered_macros$cell
  ],
  ,
  drop = FALSE
]

barcodes.macros.umap_excluded$exclusion_reason <-
  "outside_main_macro_UMAP_component"

# ============================================================
# Final retained macrophage/DC barcodes
# ============================================================

barcodes.macros <- recovered_macros[
  keep_main_component[
    recovered_macros$cell
  ],
  ,
  drop = FALSE
]

stopifnot(
  !anyDuplicated(barcodes.macros$cell),
  !any(
    barcodes.macros$cell %in%
      barcodes.macros.umap_excluded$cell
  )
)

cat("\nUMAP filtering:\n")
cat("Connection radius:", connection_radius, "\n")
cat("Recovered before filtering:", nrow(recovered_macros), "\n")
cat("Retained in main component:", nrow(barcodes.macros), "\n")
cat(
  "Excluded outside main component:",
  nrow(barcodes.macros.umap_excluded),
  "\n"
)

cat("\nExcluded isolated cells by annotation:\n")

print(
  table(
    barcodes.macros.umap_excluded$macro_annotation
  )
)

cat("\nFinal retained annotations:\n")

print(
  table(
    barcodes.macros$macro_annotation
  )
)

if (
  nrow(barcodes.macros) <
    0.90 * nrow(recovered_macros)
) {
  warning(
    "More than 10% of recovered macrophage/DC cells ",
    "were excluded by the UMAP-component filter."
  )
}

# ============================================================
# Candidate cells not annotated as macrophage/DC
# ============================================================

macro_candidates <- colnames(fibroblasts)[
  as.character(
    fibroblasts$seurat_clusters
  ) %in% c("15", "21")
]

annotation_excluded_cells <- setdiff(
  macro_candidates,
  recovered_macros$cell
)

annotation_excluded <- data.frame(
  cell = annotation_excluded_cells,
  macro_annotation = NA_character_,
  global_cluster = as.character(
    qc_meta[
      annotation_excluded_cells,
      "seurat_clusters"
    ]
  ),
  condition = as.character(
    qc_meta[
      annotation_excluded_cells,
      "orig.ident"
    ]
  ),
  tag = as.character(
    qc_meta[
      annotation_excluded_cells,
      "tag"
    ]
  ),
  UMAP_1 = global_umap[
    annotation_excluded_cells,
    1
  ],
  UMAP_2 = global_umap[
    annotation_excluded_cells,
    2
  ],
  exclusion_reason = "excluded_during_macro_annotation",
  stringsAsFactors = FALSE
)

# ============================================================
# Combine all excluded cells
# ============================================================

excluded_columns <- c(
  "cell",
  "macro_annotation",
  "global_cluster",
  "condition",
  "tag",
  "UMAP_1",
  "UMAP_2",
  "exclusion_reason"
)

barcodes.macros.excluded <- rbind(
  annotation_excluded[
    ,
    excluded_columns,
    drop = FALSE
  ],
  barcodes.macros.umap_excluded[
    ,
    excluded_columns,
    drop = FALSE
  ]
)

stopifnot(
  !anyDuplicated(
    barcodes.macros.excluded$cell
  ),
  !any(
    barcodes.macros$cell %in%
      barcodes.macros.excluded$cell
  )
)

cat("\nAll excluded macrophage/DC candidates:\n")

print(
  table(
    barcodes.macros.excluded$exclusion_reason
  )
)

# ============================================================
# Save UMAP-filter diagnostics
# ============================================================

umap_filter_summary <- data.frame(
  parameter = c(
    "k_local",
    "radius_quantile",
    "radius_multiplier",
    "connection_radius",
    "recovered_before_filtering",
    "retained_main_component",
    "excluded_UMAP_components"
  ),
  value = c(
    k_local,
    radius_quantile,
    radius_multiplier,
    connection_radius,
    nrow(recovered_macros),
    nrow(barcodes.macros),
    nrow(barcodes.macros.umap_excluded)
  ),
  stringsAsFactors = FALSE
)

write.table(
  umap_filter_summary,
  file.path(
    tables_dir,
    "Macrophage_DC_UMAP_filter_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Save barcode annotations
# ============================================================

saveRDS(
  barcodes.macros,
  file.path(
    annotations_dir,
    "barcodes.macros.rds"
  )
)

write.table(
  barcodes.macros,
  file.path(
    annotations_dir,
    "barcodes.macros.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  barcodes.macros.excluded,
  file.path(
    annotations_dir,
    "barcodes.macros.excluded.rds"
  )
)

write.table(
  barcodes.macros.excluded,
  file.path(
    annotations_dir,
    "barcodes.macros.excluded.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Create annotated macrophage/DC object
# ============================================================

fibros.CD45.Macros <- subset(
  fibroblasts,
  cells = barcodes.macros$cell
)

annotation_match <- match(
  colnames(fibros.CD45.Macros),
  barcodes.macros$cell
)

fibros.CD45.Macros$macro_annotation <-
  barcodes.macros$macro_annotation[
    annotation_match
  ]

fibros.CD45.Macros$macro_annotation <- factor(
  fibros.CD45.Macros$macro_annotation,
  levels = macro_levels
)

stopifnot(
  !anyNA(
    fibros.CD45.Macros$macro_annotation
  )
)

Idents(fibros.CD45.Macros) <- "macro_annotation"

macro_object_file <- file.path(
  objects_dir,
  "fibroblasts.Macrophage_DC.annotated.rds"
)

saveRDS(
  fibros.CD45.Macros,
  macro_object_file
)

# ============================================================
# Summary tables
# ============================================================

macro_summary <- barcodes.macros %>%
  count(
    macro_annotation,
    name = "cells"
  ) %>%
  mutate(
    percentage = 100 * cells / sum(cells)
  ) %>%
  arrange(
    factor(
      macro_annotation,
      levels = macro_levels
    )
  )

macro_by_condition <- barcodes.macros %>%
  count(
    condition,
    macro_annotation,
    name = "cells"
  ) %>%
  group_by(condition) %>%
  mutate(
    total_macrophage_DC = sum(cells),
    percentage = 100 * cells / total_macrophage_DC
  ) %>%
  ungroup()

macro_by_sample <- barcodes.macros %>%
  count(
    condition,
    tag,
    macro_annotation,
    name = "cells"
  )

write.table(
  macro_summary,
  file.path(
    tables_dir,
    "Macrophage_DC_annotation_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  macro_by_condition,
  file.path(
    tables_dir,
    "Macrophage_DC_counts_by_condition.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  macro_by_sample,
  file.path(
    tables_dir,
    "Macrophage_DC_counts_by_sample.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Colours
# ============================================================

macro_colors <- c(
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

# ============================================================
# UMAP: final macrophage/DC annotations
# ============================================================

p_macros <- DimPlot(
  object = fibros.CD45.Macros,
  reduction = "umap",
  group.by = "macro_annotation",
  cols = macro_colors,
  label = TRUE,
  repel = TRUE,
  pt.size = 0.7
) +
  NoLegend() +
  ggtitle("Macrophage and DC annotations")

ggsave(
  filename = file.path(
    results_dir,
    "Macrophage_DC_UMAP.pdf"
  ),
  plot = p_macros,
  width = 8,
  height = 6
)

# ============================================================
# UMAP split by condition
# ============================================================

p_macros_split <- DimPlot(
  object = fibros.CD45.Macros,
  reduction = "umap",
  group.by = "macro_annotation",
  split.by = "orig.ident",
  cols = macro_colors,
  pt.size = 0.6,
  ncol = 3
)

ggsave(
  filename = file.path(
    results_dir,
    "Macrophage_DC_UMAP_by_condition.pdf"
  ),
  plot = p_macros_split,
  width = 18,
  height = 6
)

# ============================================================
# Final report
# ============================================================

cat("\nDone.\n")

cat("\nAnnotated object:\n")
cat(macro_object_file, "\n")

cat("\nRetained macrophage/DC cells:\n")
cat(ncol(fibros.CD45.Macros), "\n")

cat("\nFinal annotation counts:\n")
print(
  table(
    fibros.CD45.Macros$macro_annotation
  )
)

cat("\nExcluded cells:\n")
print(
  table(
    barcodes.macros.excluded$exclusion_reason
  )
)