# ============================================================
# 06_Tcells.R
# Selección, reclustering y anotación de T cells
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(scCustomize)
  library(RPresto)
})

set.seed(42)

# ============================================================
# Configuración
# ============================================================

config <- read.table(
  "config.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("key", "value") %in% colnames(config))) {
  stop("config.tsv debe contener las columnas 'key' y 'value'.")
}

values <- setNames(config$value, config$key)

required_config <- c("rsdir", "datadir", "outdir")
missing_config <- setdiff(required_config, names(values))

if (length(missing_config) > 0) {
  stop(
    "Faltan entradas en config.tsv: ",
    paste(missing_config, collapse = ", ")
  )
}

rsdir <- unname(values[["rsdir"]])
datadir <- unname(values[["datadir"]])
outdir <- unname(values[["outdir"]])

source("R/Functions.R")
rsdir <- unname(values[["rsdir"]])
datadir <- unname(values[["datadir"]])
outdir <- unname(values[["outdir"]])

source("R/Functions.R")

# ============================================================
# Directorios
# ============================================================

objects_dir <- file.path(
  rsdir,
  "objects"
)

plots_dir <- file.path(
  outdir,
  "CD45+",
  "Tcells"
)

tables_dir <- file.path(
  rsdir,
  "tables",
  "CD45+",
  "Tcells"
)

for (directory in c(
  objects_dir,
  plots_dir,
  tables_dir
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

cat("\nDirectorios del análisis de T cells:\n")
cat("Objetos: ", objects_dir, "\n")
cat("Plots:   ", plots_dir, "\n")
cat("Tablas:  ", tables_dir, "\n")

# ============================================================
# Funciones auxiliares
# ============================================================


save_plot_pdf_png <- function(
    plot,
    filename,
    width,
    height
) {
  ggsave(
    filename = file.path(plots_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    useDingbats = FALSE,
    limitsize = FALSE
  )

  ggsave(
    filename = file.path(plots_dir, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
}

find_markers_and_top <- function(
    object,
    identity_column,
    top_n,
    prefix,
    min_pct = 0.20,
    logfc_threshold = 0.25
) {
  DefaultAssay(object) <- "RNA"
  object <- JoinLayers(object, assay = "RNA")
  Idents(object) <- identity_column

  markers <- FindAllMarkers(
    object = object,
    assay = "RNA",
    only.pos = TRUE,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    test.use = "wilcox"
  )

  markers <- markers %>%
    arrange(cluster, desc(avg_log2FC))

  top_markers <- markers %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    slice_max(
      order_by = avg_log2FC,
      n = top_n,
      with_ties = FALSE
    ) %>%
    ungroup()

  write.table(
    markers,
    file = file.path(
      tables_dir,
      paste0(prefix, "_FindAllMarkers_all.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    top_markers,
    file = file.path(
      tables_dir,
      paste0(prefix, "_FindAllMarkers_top", top_n, ".tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  list(
    object = object,
    markers = markers,
    top_markers = top_markers
  )
}

validate_annotation_map <- function(observed, annotation_map, object_name) {
  observed <- sort(unique(as.character(observed)))
  missing_annotations <- setdiff(observed, names(annotation_map))

  if (length(missing_annotations) > 0) {
    stop(
      "Faltan anotaciones para ",
      object_name,
      ": ",
      paste(missing_annotations, collapse = ", ")
    )
  }
}
# ============================================================
# Leer objeto de QC
# ============================================================

input_file <- file.path(
  objects_dir,
  "fibroblasts.QC.rds"
)

if (!file.exists(input_file)) {
  stop(
    "No se encuentra el objeto de QC: ",
    input_file
  )
}

data <- readRDS(
  input_file
)

cat("\nObjeto de entrada cargado:\n")
cat(input_file, "\n")

# ============================================================
# Guardar clustering global original
# ============================================================

data$global_clusters <- as.character(data$seurat_clusters)

# ============================================================
# Subsets originales
# ============================================================

data <- subset(
  data,
  subset = (
    CD == "CD45+" |
      seurat_clusters == 21 |
      seurat_clusters == 15
  )
)

data <- subset(
  data,
  subset = (
    seurat_clusters == 0 |
      seurat_clusters == 1 |
      seurat_clusters == 7 |
      seurat_clusters == 8 |
      seurat_clusters == 6 |
      seurat_clusters == 19 |
      seurat_clusters == 14
  )
)

Idents(data) <- "seurat_clusters"

cat("\nCélulas después de los subsets originales:\n")
print(ncol(data))

cat("\nClusters globales seleccionados:\n")
print(table(data$global_clusters))

cat("\nCélulas por genotipo:\n")
print(table(data$orig.ident))

# ============================================================
# UMAP del subset con coordenadas globales
# ============================================================

p_subset <- DimPlot(
  data,
  reduction = "umap",
  group.by = "global_clusters",
  label = TRUE,
  repel = TRUE,
  label.size = 4
) +
  NoLegend() +
  ggtitle("Initial T-cell subset: global clusters")

save_plot_pdf_png(
  plot = p_subset,
  filename = "Tcells_initial_subset",
  width = 8,
  height = 7
)

# ============================================================
# Separar genotipos
# ============================================================

WT <- subset(data, subset = orig.ident == "WT")
CD169DTR <- subset(data, subset = orig.ident == "CD169DTR")
CCR2KO <- subset(data, subset = orig.ident == "CCR2KO")

objects_by_genotype <- list(
  WT = WT,
  CD169DTR = CD169DTR,
  CCR2KO = CCR2KO
)

empty_objects <- names(objects_by_genotype)[
  vapply(objects_by_genotype, ncol, numeric(1)) == 0
]

if (length(empty_objects) > 0) {
  stop(
    "No hay células para: ",
    paste(empty_objects, collapse = ", ")
  )
}

cat("\nCélulas por objeto:\n")
print(vapply(objects_by_genotype, ncol, numeric(1)))

# ============================================================
# Primer reclustering por genotipo
# ============================================================

recluster_tcells <- function(
    object,
    resolution,
    object_name,
    seed = 42
) {
  DefaultAssay(object) <- "RNA"

  if (length(VariableFeatures(object)) == 0) {
    stop(
      object_name,
      " no conserva genes variables. No se puede reproducir ",
      "el PCA original sin recalcularlos."
    )
  }

  object <- ScaleData(
    object = object,
    features = VariableFeatures(object),
    verbose = TRUE
  )

  object <- RunPCA(
    object = object,
    features = VariableFeatures(object),
    npcs = 50,
    seed.use = seed,
    verbose = TRUE
  )

  available_pcs <- ncol(
    Embeddings(object, reduction = "pca")
  )

  if (available_pcs < 50) {
    stop(
      object_name,
      " solo ha generado ",
      available_pcs,
      " PCs; el workflow histórico necesita 50."
    )
  }

  object <- FindNeighbors(
    object = object,
    reduction = "pca",
    dims = 1:50,
    verbose = TRUE
  )

  object <- FindClusters(
    object = object,
    resolution = resolution,
    random.seed = seed,
    verbose = TRUE
  )

  object <- RunUMAP(
    object = object,
    reduction = "pca",
    dims = 1:50,
    seed.use = seed,
    verbose = TRUE
  )

  object$initial_tcell_cluster <- as.character(
    object$seurat_clusters
  )

  Idents(object) <- "initial_tcell_cluster"
  object
}

WT <- recluster_tcells(
  object = WT,
  resolution = 1.3,
  object_name = "WT"
)

CD169DTR <- recluster_tcells(
  object = CD169DTR,
  resolution = 1.0,
  object_name = "CD169DTR"
)

CCR2KO <- recluster_tcells(
  object = CCR2KO,
  resolution = 1.5,
  object_name = "CCR2KO"
)

cat("\nClusters iniciales de WT:\n")
print(table(WT$initial_tcell_cluster))

cat("\nClusters iniciales de CD169DTR:\n")
print(table(CD169DTR$initial_tcell_cluster))

cat("\nClusters iniciales de CCR2KO:\n")
print(table(CCR2KO$initial_tcell_cluster))

# ============================================================
# UMAP del primer reclustering
# ============================================================

make_cluster_plot <- function(object, title) {
  DimPlot(
    object,
    reduction = "umap",
    group.by = "initial_tcell_cluster",
    label = TRUE,
    repel = TRUE,
    label.size = 4
  ) +
    NoLegend() +
    ggtitle(title)
}

p_WT <- make_cluster_plot(
  WT,
  "WT: initial T-cell reclustering"
)

p_CD169DTR <- make_cluster_plot(
  CD169DTR,
  "CD169DTR: initial T-cell reclustering"
)

p_CCR2KO <- make_cluster_plot(
  CCR2KO,
  "CCR2KO: initial T-cell reclustering"
)

p_all <- cowplot::plot_grid(
  p_WT,
  p_CD169DTR,
  p_CCR2KO,
  nrow = 1,
  align = "h"
)

save_plot_pdf_png(
  plot = p_all,
  filename = "Tcells_initial_reclustering_by_genotype",
  width = 24,
  height = 7
)

# ============================================================
# Marcadores del primer reclustering
# ============================================================

marker_result_WT <- find_markers_and_top(
  object = WT,
  identity_column = "seurat_clusters",
  top_n = 20,
  prefix = "WT_Tcells"
)
WT <- marker_result_WT$object
markers.WT <- marker_result_WT$markers
table.WT <- marker_result_WT$top_markers

marker_result_CD169DTR <- find_markers_and_top(
  object = CD169DTR,
  identity_column = "seurat_clusters",
  top_n = 20,
  prefix = "CD169DTR_Tcells"
)
CD169DTR <- marker_result_CD169DTR$object
markers.CD169DTR <- marker_result_CD169DTR$markers
table.CD169DTR <- marker_result_CD169DTR$top_markers

marker_result_CCR2KO <- find_markers_and_top(
  object = CCR2KO,
  identity_column = "seurat_clusters",
  top_n = 20,
  prefix = "CCR2KO_Tcells"
)
CCR2KO <- marker_result_CCR2KO$object
markers.CCR2KO <- marker_result_CCR2KO$markers
table.CCR2KO <- marker_result_CCR2KO$top_markers

cat("\nTop markers WT:\n")
print(table.WT)

cat("\nTop markers CD169DTR:\n")
print(table.CD169DTR)

cat("\nTop markers CCR2KO:\n")
print(table.CCR2KO)

# ============================================================
# Anotación inicial por genotipo
# ============================================================

WT$Cluster <- dplyr::recode(
  as.character(WT$seurat_clusters),
  "0"  = "Immature T Cells",
  "1"  = "Immature T Cells",
  "2"  = "Immature T Cells",
  "3"  = "Stem like T Cells",
  "4"  = "Tpex Transitional",
  "5"  = "CD4 Naive",
  "6"  = "Stem like T Cells",
  "7"  = "CD4 Early Activated",
  "8"  = "CD8 Memory",
  "9"  = "CD4 Naive",
  "10" = "NK",
  "11" = "Cd4 IL7R Memory Like",
  "12" = "Treg",
  "13" = "CD8 Memory",
  "14" = "Tgd",
  "15" = "Stem like T Cells",
  "16" = "NKT",
  "17" = "Cd8 Exhausted",
  "18" = "Tgd",
  "19" = "Remove",
  "20" = "Remove"
)

CD169DTR$Cluster <- dplyr::recode(
  as.character(CD169DTR$seurat_clusters),
  "0"  = "Immature T Cells",
  "1"  = "Immature T Cells",
  "2"  = "Immature T Cells",
  "3"  = "Stem like T Cells",
  "4"  = "CD4 Naive",
  "5"  = "Tpex Transitional",
  "6"  = "Stem like T Cells",
  "7"  = "CD4 Early Activated",
  "8"  = "CD8 Memory",
  "9"  = "NK",
  "10" = "Refine",
  "11" = "Tgd",
  "12" = "Tpex Transitional"
)

CCR2KO$Cluster <- dplyr::recode(
  as.character(CCR2KO$seurat_clusters),
  "0"  = "Immature T Cells",
  "1"  = "Immature T Cells",
  "2"  = "CD4 Naive",
  "3"  = "Refine",
  "4"  = "Tpex Transitional",
  "5"  = "NK",
  "6"  = "CD8 Memory",
  "7"  = "Treg",
  "8"  = "Cd4 IL7R Memory Like",
  "9"  = "CD8 Memory",
  "10" = "Stem like T Cells",
  "11" = "Stem like T Cells",
  "12" = "Cd8 Exhausted",
  "13" = "Stem like T Cells",
  "14" = "NKT",
  "15" = "CD4 Naive",
  "16" = "Tgd",
  "17" = "Remove"
)

cat("\nAnotación inicial de WT:\n")
print(table(WT$seurat_clusters, WT$Cluster))

cat("\nAnotación inicial de CD169DTR:\n")
print(table(CD169DTR$seurat_clusters, CD169DTR$Cluster))

cat("\nAnotación inicial de CCR2KO:\n")
print(table(CCR2KO$seurat_clusters, CCR2KO$Cluster))

# ============================================================
# Refinamiento del cluster 3 de CCR2KO
# Subcluster 0 = CD8 Memory
# Subcluster 1 = CD4 Naive
# ============================================================

Idents(CCR2KO) <- "seurat_clusters"

CCR2KO_cluster3 <- subset(
  CCR2KO,
  idents = "3"
)

DefaultAssay(CCR2KO_cluster3) <- "RNA"

cat(
  "\nCélulas del cluster 3 de CCR2KO: ",
  ncol(CCR2KO_cluster3),
  "\n",
  sep = ""
)

set.seed(42)

CCR2KO_cluster3 <- FindVariableFeatures(
  CCR2KO_cluster3,
  selection.method = "vst",
  nfeatures = 2000
)

CCR2KO_cluster3 <- ScaleData(
  CCR2KO_cluster3,
  features = VariableFeatures(CCR2KO_cluster3)
)

CCR2KO_cluster3 <- RunPCA(
  CCR2KO_cluster3,
  features = VariableFeatures(CCR2KO_cluster3),
  npcs = 30,
  seed.use = 42
)

available_pcs_cluster3 <- ncol(
  Embeddings(CCR2KO_cluster3, reduction = "pca")
)

dims_cluster3 <- seq_len(
  min(15, available_pcs_cluster3)
)

CCR2KO_cluster3 <- FindNeighbors(
  CCR2KO_cluster3,
  reduction = "pca",
  dims = dims_cluster3
)

CCR2KO_cluster3 <- FindClusters(
  CCR2KO_cluster3,
  resolution = 0.4,
  random.seed = 42
)

CCR2KO_cluster3 <- RunUMAP(
  CCR2KO_cluster3,
  reduction = "pca",
  dims = dims_cluster3,
  seed.use = 42
)

CCR2KO_cluster3$cluster3_subcluster <- as.character(
  CCR2KO_cluster3$seurat_clusters
)

Idents(CCR2KO_cluster3) <- "cluster3_subcluster"

cat("\nSubclusters del cluster 3 de CCR2KO:\n")
print(table(CCR2KO_cluster3$cluster3_subcluster))

p_CCR2KO_cluster3 <- DimPlot(
  CCR2KO_cluster3,
  reduction = "umap",
  group.by = "cluster3_subcluster",
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  pt.size = 0.7
) +
  NoLegend() +
  ggtitle("CCR2KO cluster 3 refinement") +
  theme_classic()

save_plot_pdf_png(
  plot = p_CCR2KO_cluster3,
  filename = "CCR2KO_cluster3_refinement_UMAP",
  width = 7,
  height = 6
)

cluster3_features <- c(
  "Cd4",
  "Cd8a",
  "Cd8b1",
  "Ccr7",
  "Sell",
  "Il7r",
  "Tcf7",
  "Lef1",
  "Klf2",
  "Bcl2",
  "Ccl5",
  "Nkg7",
  "Gzmk",
  "Ctsw",
  "Prf1",
  "Gzmb",
  "Socs3",
  "Cxcl2"
)

cluster3_features <- intersect(
  cluster3_features,
  rownames(CCR2KO_cluster3)
)

p_CCR2KO_cluster3_dotplot <- DotPlot(
  CCR2KO_cluster3,
  features = cluster3_features,
  group.by = "cluster3_subcluster",
  dot.scale = 7
) +
  RotatedAxis() +
  labs(
    x = NULL,
    y = "Cluster 3 subcluster",
    color = "Average expression",
    size = "Percent expressed"
  ) +
  theme_classic(base_size = 11)

save_plot_pdf_png(
  plot = p_CCR2KO_cluster3_dotplot,
  filename = "CCR2KO_cluster3_refinement_DotPlot",
  width = 12,
  height = 5
)

marker_result_CCR2KO_cluster3 <- find_markers_and_top(
  object = CCR2KO_cluster3,
  identity_column = "cluster3_subcluster",
  top_n = 30,
  prefix = "CCR2KO_cluster3"
)

CCR2KO_cluster3 <- marker_result_CCR2KO_cluster3$object
markers.CCR2KO_cluster3 <- marker_result_CCR2KO_cluster3$markers
top_markers.CCR2KO_cluster3 <- marker_result_CCR2KO_cluster3$top_markers

cluster3_annotation_map <- c(
  "0" = "CD8 Memory",
  "1" = "CD4 Naive"
)

validate_annotation_map(
  observed = CCR2KO_cluster3$cluster3_subcluster,
  annotation_map = cluster3_annotation_map,
  object_name = "los subclusters del cluster 3 de CCR2KO"
)

CCR2KO_cluster3$Cluster_refined <- unname(
  cluster3_annotation_map[
    as.character(CCR2KO_cluster3$cluster3_subcluster)
  ]
)

cells_cluster3 <- intersect(
  colnames(CCR2KO_cluster3),
  colnames(CCR2KO)
)

CCR2KO$cluster3_subcluster <- NA_character_

CCR2KO@meta.data[
  cells_cluster3,
  "cluster3_subcluster"
] <- as.character(
  CCR2KO_cluster3@meta.data[
    cells_cluster3,
    "cluster3_subcluster"
  ]
)

CCR2KO@meta.data[
  cells_cluster3,
  "Cluster"
] <- as.character(
  CCR2KO_cluster3@meta.data[
    cells_cluster3,
    "Cluster_refined"
  ]
)

if (any(CCR2KO$Cluster == "Refine")) {
  stop("Quedan células del cluster 3 de CCR2KO sin anotar.")
}

cat("\nAnotación final del antiguo cluster 3 de CCR2KO:\n")
print(
  table(
    CCR2KO@meta.data[cells_cluster3, "cluster3_subcluster"],
    CCR2KO@meta.data[cells_cluster3, "Cluster"]
  )
)

CCR2KO_cluster3_annotation <- data.frame(
  barcode = cells_cluster3,
  cluster3_subcluster = CCR2KO@meta.data[
    cells_cluster3,
    "cluster3_subcluster"
  ],
  Cluster = CCR2KO@meta.data[
    cells_cluster3,
    "Cluster"
  ],
  stringsAsFactors = FALSE
)

write.table(
  CCR2KO_cluster3_annotation,
  file = file.path(
    tables_dir,
    "CCR2KO_cluster3_refined_annotation.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Refinamiento del cluster 10 de CD169DTR
# ============================================================

Idents(CD169DTR) <- "seurat_clusters"

DTR_cluster10 <- subset(
  CD169DTR,
  idents = "10"
)

DefaultAssay(DTR_cluster10) <- "RNA"

set.seed(42)

DTR_cluster10 <- FindVariableFeatures(
  DTR_cluster10,
  selection.method = "vst",
  nfeatures = 2000
)

DTR_cluster10 <- ScaleData(
  DTR_cluster10,
  features = VariableFeatures(DTR_cluster10)
)

DTR_cluster10 <- RunPCA(
  DTR_cluster10,
  features = VariableFeatures(DTR_cluster10),
  npcs = 30,
  seed.use = 42
)

available_pcs_DTR10 <- ncol(
  Embeddings(DTR_cluster10, reduction = "pca")
)

dims_DTR10 <- seq_len(min(20, available_pcs_DTR10))

DTR_cluster10 <- FindNeighbors(
  DTR_cluster10,
  reduction = "pca",
  dims = dims_DTR10
)

DTR_cluster10 <- FindClusters(
  DTR_cluster10,
  resolution = 1.3,
  random.seed = 42
)

DTR_cluster10 <- RunUMAP(
  DTR_cluster10,
  reduction = "pca",
  dims = dims_DTR10,
  seed.use = 42
)

DTR_cluster10$dtr_cluster10_subcluster <- as.character(
  DTR_cluster10$seurat_clusters
)

Idents(DTR_cluster10) <- "dtr_cluster10_subcluster"

p_DTR_cluster10 <- DimPlot(
  DTR_cluster10,
  reduction = "umap",
  group.by = "dtr_cluster10_subcluster",
  label = TRUE,
  repel = TRUE
) +
  NoLegend() +
  ggtitle("CD169DTR cluster 10 refinement")

save_plot_pdf_png(
  plot = p_DTR_cluster10,
  filename = "DTR_cluster10_refinement",
  width = 7,
  height = 6
)

DTR_cluster10_features <- c(
  "Cd4",
  "Il7r",
  "Ccr7",
  "Ccr6",
  "Il17re",
  "Foxp3",
  "Ctla4",
  "Il2ra",
  "Lrrc32",
  "Trdc",
  "Tcrg-C1",
  "Zbtb16"
)

DTR_cluster10_features <- intersect(
  DTR_cluster10_features,
  rownames(DTR_cluster10)
)

p_DTR_cluster10_dotplot <- DotPlot(
  DTR_cluster10,
  features = DTR_cluster10_features,
  group.by = "dtr_cluster10_subcluster"
) +
  RotatedAxis() +
  theme_classic(base_size = 11)

save_plot_pdf_png(
  plot = p_DTR_cluster10_dotplot,
  filename = "DTR_cluster10_annotation_DotPlot",
  width = 9,
  height = 5
)

marker_result_DTR_cluster10 <- find_markers_and_top(
  object = DTR_cluster10,
  identity_column = "dtr_cluster10_subcluster",
  top_n = 50,
  prefix = "DTR_cluster10"
)

DTR_cluster10 <- marker_result_DTR_cluster10$object
markers.DTR_cluster10 <- marker_result_DTR_cluster10$markers
table.DTR_cluster10 <- marker_result_DTR_cluster10$top_markers

DTR_cluster10_annotation_map <- c(
  "0" = "Cd4 IL7R Memory Like",
  "1" = "Treg",
  "2" = "Remove",
  "3" = "Tgd",
  "4" = "Cd4 IL7R Memory Like"
)

validate_annotation_map(
  observed = DTR_cluster10$dtr_cluster10_subcluster,
  annotation_map = DTR_cluster10_annotation_map,
  object_name = "los subclusters del cluster 10 de CD169DTR"
)

DTR_cluster10$Cluster_refined <- unname(
  DTR_cluster10_annotation_map[
    as.character(DTR_cluster10$dtr_cluster10_subcluster)
  ]
)

cells_DTR10 <- intersect(
  colnames(DTR_cluster10),
  colnames(CD169DTR)
)

CD169DTR$refined_cluster <- NA_character_

CD169DTR@meta.data[
  cells_DTR10,
  "refined_cluster"
] <- as.character(
  DTR_cluster10@meta.data[
    cells_DTR10,
    "dtr_cluster10_subcluster"
  ]
)

CD169DTR@meta.data[
  cells_DTR10,
  "Cluster"
] <- as.character(
  DTR_cluster10@meta.data[
    cells_DTR10,
    "Cluster_refined"
  ]
)

if (any(CD169DTR$Cluster == "Refine")) {
  stop("Quedan células del cluster 10 de CD169DTR sin anotar.")
}

cat("\nAnotación final del antiguo cluster 10 de CD169DTR:\n")
print(
  table(
    CD169DTR@meta.data[cells_DTR10, "refined_cluster"],
    CD169DTR@meta.data[cells_DTR10, "Cluster"]
  )
)

# ============================================================
# Eliminar poblaciones Remove
# ============================================================

WT <- subset(WT, subset = Cluster != "Remove")
CD169DTR <- subset(CD169DTR, subset = Cluster != "Remove")
CCR2KO <- subset(CCR2KO, subset = Cluster != "Remove")

Idents(WT) <- "Cluster"
Idents(CD169DTR) <- "Cluster"
Idents(CCR2KO) <- "Cluster"

# ============================================================
# Plot final con UMAP individuales
# ============================================================

p_WT_final <- DimPlot(
  WT,
  reduction = "umap",
  group.by = "Cluster",
  label = TRUE,
  repel = TRUE
) +
  NoLegend() +
  ggtitle("WT")

p_CD169DTR_final <- DimPlot(
  CD169DTR,
  reduction = "umap",
  group.by = "Cluster",
  label = TRUE,
  repel = TRUE
) +
  NoLegend() +
  ggtitle("CD169DTR")

p_CCR2KO_final <- DimPlot(
  CCR2KO,
  reduction = "umap",
  group.by = "Cluster",
  label = TRUE,
  repel = TRUE
) +
  NoLegend() +
  ggtitle("CCR2KO")

p_Tcells_annotation <- cowplot::plot_grid(
  p_WT_final,
  p_CD169DTR_final,
  p_CCR2KO_final,
  nrow = 1,
  align = "h"
)

save_plot_pdf_png(
  plot = p_Tcells_annotation,
  filename = "Tcells_final_annotation_split",
  width = 21,
  height = 7
)

# ============================================================
# Unir los tres genotipos
# ============================================================

Tcells <- merge(
  WT,
  y = list(CD169DTR, CCR2KO),
  merge.data = TRUE
)

Idents(Tcells) <- "Cluster"

cat("\nAnotación final por genotipo:\n")
print(table(Tcells$orig.ident, Tcells$Cluster))

if (anyNA(Tcells$Cluster)) {
  stop("El objeto final contiene células sin anotación Cluster.")
}

if (any(Tcells$Cluster %in% c("Refine", "Remove"))) {
  stop("El objeto final conserva etiquetas Refine o Remove.")
}


# ============================================================
# Guardar objeto y anotación
# ============================================================

output_object <- file.path(
  objects_dir,
  "fibroblasts.CD45.Tcells.annotated.rds"
)

saveRDS(
  Tcells,
  output_object
)

Tcells_annotation <- data.frame(
  barcode = colnames(Tcells),
  orig.ident = as.character(Tcells$orig.ident),
  seurat_clusters = as.character(Tcells$seurat_clusters),
  Cluster = as.character(Tcells$Cluster),
  stringsAsFactors = FALSE
)

output_annotation <- file.path(
  tables_dir,
  "fibroblasts.CD45.Tcells.annotation.tsv"
)

write.table(
  Tcells_annotation,
  file = output_annotation,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    tables_dir,
    "06_Tcells_sessionInfo.txt"
  )
)

cat(
  "\nAnotación final de T cells guardada correctamente.\n",
  "Objeto: ", output_object, "\n",
  "Anotación: ", output_annotation, "\n",
  "Plots: ", plots_dir, "\n",
  "Tablas: ", tables_dir, "\n",
  sep = ""
)