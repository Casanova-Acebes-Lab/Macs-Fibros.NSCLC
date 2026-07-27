# ============================================================
# 06_Tcells.R
# Selección inicial y primer reclustering de T cells
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(scCustomize)
library(RPresto)

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

values <- config$value
names(values) <- config$key

rsdir <- values["rsdir"]
datadir <- values["datadir"]
outdir <- values["outdir"]

source("R/Functions.R")

objects_dir <- file.path(rsdir, "objects")
plots_dir <- file.path(outdir, "CD45+", "Tcells")

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

# ============================================================
# Leer objeto original
# ============================================================

input_file <- file.path(
  objects_dir,
  "fibroblasts.Not.Integrated.rds"
)

if (!file.exists(input_file)) {
  stop("No se encuentra el objeto: ", input_file)
}

data <- readRDS(input_file)

required_metadata <- c(
  "CD",
  "seurat_clusters",
  "orig.ident"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(data@meta.data)
)

if (length(missing_metadata) > 0) {
  stop(
    "Faltan columnas de metadatos: ",
    paste(missing_metadata, collapse = ", ")
  )
}

# ============================================================
# Guardar el clustering global original
#
# FindClusters() sobrescribirá seurat_clusters posteriormente.
# ============================================================

data$global_clusters <- as.character(
  data$seurat_clusters
)

# ============================================================
# Primer subset original
#
# Se mantiene literalmente para reproducir la selección usada
# en el workflow histórico.
# ============================================================

data <- subset(
  data,
  subset = (
    CD == "CD45+" |
      seurat_clusters == 21 |
      seurat_clusters == 15
  )
)

# ============================================================
# Segundo subset original
#
# Este segundo filtro elimina finalmente los clusters 15 y 21.
# El conjunto resultante contiene:
# 0, 1, 6, 7, 8, 14 y 19.
# ============================================================

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
# UMAP del subset usando las coordenadas globales originales
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

ggsave(
  filename = file.path(
    plots_dir,
    "Tcells_initial_subset.pdf"
  ),
  plot = p_subset,
  width = 8,
  height = 7
)

ggsave(
  filename = file.path(
    plots_dir,
    "Tcells_initial_subset.png"
  ),
  plot = p_subset,
  width = 8,
  height = 7,
  dpi = 300
)

# ============================================================
# Separar los tres genotipos
# ============================================================

WT <- subset(
  data,
  subset = orig.ident == "WT"
)

CD169DTR <- subset(
  data,
  subset = orig.ident == "CD169DTR"
)

CCR2KO <- subset(
  data,
  subset = orig.ident == "CCR2KO"
)

objects_by_genotype <- list(
  WT = WT,
  CD169DTR = CD169DTR,
  CCR2KO = CCR2KO
)

empty_objects <- names(objects_by_genotype)[
  vapply(
    objects_by_genotype,
    ncol,
    numeric(1)
  ) == 0
]

if (length(empty_objects) > 0) {
  stop(
    "No hay células para: ",
    paste(empty_objects, collapse = ", ")
  )
}

cat("\nCélulas por objeto:\n")

print(
  vapply(
    objects_by_genotype,
    ncol,
    numeric(1)
  )
)

# ============================================================
# Función para reproducir el primer reclustering
#
# No se ejecuta FindVariableFeatures(), porque el workflow
# original reutilizaba las variables ya presentes en el objeto.
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
      " no conserva genes variables. ",
      "No se puede reproducir el PCA original sin decidir ",
      "si deben recalcularse."
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
    Embeddings(
      object,
      reduction = "pca"
    )
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

  # Guardar explícitamente el resultado del primer reclustering.
  object$initial_tcell_cluster <- as.character(
    object$seurat_clusters
  )

  Idents(object) <- "initial_tcell_cluster"

  return(object)
}

# ============================================================
# Primer reclustering por genotipo
#
# Resoluciones originales:
#   WT        = 1.3
#   CD169DTR  = 1.0
#   CCR2KO    = 1.5
# ============================================================

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

ggsave(
  filename = file.path(
    plots_dir,
    "Tcells_initial_reclustering_by_genotype.pdf"
  ),
  plot = p_all,
  width = 24,
  height = 7,
  limitsize = FALSE
)





# ============================================================
# Marcadores del primer reclustering
# ============================================================

DefaultAssay(WT) <- "RNA"
DefaultAssay(CD169DTR) <- "RNA"
DefaultAssay(CCR2KO) <- "RNA"

WT <- JoinLayers(WT, assay = "RNA")
CD169DTR <- JoinLayers(CD169DTR, assay = "RNA")
CCR2KO <- JoinLayers(CCR2KO, assay = "RNA")

# WT
Idents(WT) <- "seurat_clusters"

markers.WT <- FindAllMarkers(
  WT,
  only.pos = TRUE,
  min.pct = 0.2,
  logfc.threshold = 0.25
)

table.WT <- markers.WT %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = abs(avg_log2FC),
    n = 20,
    with_ties = TRUE
  ) %>%
  dplyr::ungroup() %>%
  as.data.frame()

table.WT

# CD169DTR
Idents(CD169DTR) <- "seurat_clusters"

markers.CD169DTR <- FindAllMarkers(
  CD169DTR,
  only.pos = TRUE,
  min.pct = 0.2,
  logfc.threshold = 0.25
)

table.CD169DTR <- markers.CD169DTR %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = abs(avg_log2FC),
    n = 20,
    with_ties = TRUE
  ) %>%
  dplyr::ungroup() %>%
  as.data.frame()

table.CD169DTR

# CCR2KO
Idents(CCR2KO) <- "seurat_clusters"

markers.CCR2KO <- FindAllMarkers(
  CCR2KO,
  only.pos = TRUE,
  min.pct = 0.2,
  logfc.threshold = 0.25
)

table.CCR2KO <- markers.CCR2KO %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = abs(avg_log2FC),
    n = 20,
    with_ties = TRUE
  ) %>%
  dplyr::ungroup() %>%
  as.data.frame()

table.CCR2KO





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

Idents(WT) <- "Cluster"

table(WT$seurat_clusters, WT$Cluster)



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

table(CD169DTR$seurat_clusters, CD169DTR$Cluster)



CCR2KO$Cluster <- dplyr::recode(
  as.character(CCR2KO$seurat_clusters),
  "0"  = "Immature T Cells",
  "1"  = "Immature T Cells",
  "2"  = "CD4 Naive",
  "3"  = "CD4 Naive",
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

table(CCR2KO$seurat_clusters, CCR2KO$Cluster)




# ============================================================
# Refinement del cluster 10 de CD169DTR
# ============================================================

Idents(CD169DTR) <- "seurat_clusters"

DTR_cluster10 <- subset(
  CD169DTR,
  idents = "10"
)

DefaultAssay(DTR_cluster10) <- "RNA"

DTR_cluster10 <- FindVariableFeatures(DTR_cluster10)
DTR_cluster10 <- ScaleData(DTR_cluster10)
DTR_cluster10 <- RunPCA(DTR_cluster10, npcs = 30)

DTR_cluster10 <- FindNeighbors(
  DTR_cluster10,
  dims = 1:20
)

DTR_cluster10 <- FindClusters(
  DTR_cluster10,
  resolution = 1.3
)

DTR_cluster10 <- RunUMAP(
  DTR_cluster10,
  dims = 1:20
)

DimPlot(
  DTR_cluster10,
  label = TRUE,
  repel = TRUE
) +
  NoLegend()

# Marcadores
Idents(DTR_cluster10) <- "seurat_clusters"

markers.DTR_cluster10 <- FindAllMarkers(
  DTR_cluster10,
  only.pos = TRUE,
  min.pct = 0.2,
  logfc.threshold = 0.25
)

table.DTR_cluster10 <- markers.DTR_cluster10 %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = abs(avg_log2FC),
    n = 50,
    with_ties = TRUE
  ) %>%
  dplyr::ungroup() %>%
  as.data.frame()

View(table.DTR_cluster10)




# ============================================================
# Refinement del cluster 10 de CD169DTR
# ============================================================

Idents(CD169DTR) <- "seurat_clusters"

DTR_cluster10 <- subset(
  CD169DTR,
  idents = "10"
)

DefaultAssay(DTR_cluster10) <- "RNA"

DTR_cluster10 <- FindVariableFeatures(DTR_cluster10)
DTR_cluster10 <- ScaleData(DTR_cluster10)
DTR_cluster10 <- RunPCA(DTR_cluster10, npcs = 30)

DTR_cluster10 <- FindNeighbors(
  DTR_cluster10,
  dims = 1:20
)

DTR_cluster10 <- FindClusters(
  DTR_cluster10,
  resolution = 1.3
)

DTR_cluster10 <- RunUMAP(
  DTR_cluster10,
  dims = 1:20
)

# UMAP
p_DTR_cluster10 <- DimPlot(
  DTR_cluster10,
  label = TRUE,
  repel = TRUE
) +
  NoLegend()

p_DTR_cluster10

ggsave(
  file.path(plots_dir, "DTR_cluster10_refinement.pdf"),
  plot = p_DTR_cluster10,
  width = 7,
  height = 6
)

# Marcadores
Idents(DTR_cluster10) <- "seurat_clusters"

markers.DTR_cluster10 <- FindAllMarkers(
  DTR_cluster10,
  only.pos = TRUE,
  min.pct = 0.2,
  logfc.threshold = 0.25
)

table.DTR_cluster10 <- markers.DTR_cluster10 %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = abs(avg_log2FC),
    n = 50,
    with_ties = TRUE
  ) %>%
  dplyr::ungroup() %>%
  as.data.frame()


table.DTR_cluster10



DotPlot(
  DTR_cluster10,
  features = c(
    "Cd4", "Il7r", "Ccr7", "Ccr6", "Il17re",
    "Foxp3", "Ctla4", "Il2ra", "Lrrc32",
    "Trdc", "Tcrg-C1", "Zbtb16"
  )
) +
  RotatedAxis()

ggsave(
  file.path(plots_dir, "DTR_cluster10_annotation_DotPlot.pdf"),
  width = 9,
  height = 5
)



DTR_cluster10$Cluster_refined <- dplyr::recode(
  as.character(DTR_cluster10$seurat_clusters),
  "0" = "Cd4 IL7R Memory Like",
  "1" = "Treg",
  "2" = "Remove",
  "3" = "Tgd",
  "4" = "Cd4 IL7R Memory Like"
)





# ============================================================
# Transferir refinement al objeto CD169DTR
# ============================================================

cells_common <- intersect(
  colnames(DTR_cluster10),
  colnames(CD169DTR)
)

CD169DTR$refined_cluster <- NA_character_

CD169DTR@meta.data[cells_common, "refined_cluster"] <-
  as.character(
    DTR_cluster10@meta.data[cells_common, "seurat_clusters"]
  )

CD169DTR@meta.data[cells_common, "Cluster"] <-
  as.character(
    DTR_cluster10@meta.data[cells_common, "Cluster_refined"]
  )

table(CD169DTR$Cluster, useNA = "ifany")


# ============================================================
# Eliminar poblaciones Remove
# ============================================================

WT <- subset(WT, subset = Cluster != "Remove")
CD169DTR <- subset(CD169DTR, subset = Cluster != "Remove")
CCR2KO <- subset(CCR2KO, subset = Cluster != "Remove")


# ============================================================
# Plot final usando los UMAP individuales
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

p_Tcells_annotation

ggsave(
  file.path(plots_dir, "Tcells_final_annotation_split.pdf"),
  plot = p_Tcells_annotation,
  width = 21,
  height = 7,
  limitsize = FALSE
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

table(Tcells$orig.ident, Tcells$Cluster)


# ============================================================
# Guardar objeto y anotación
# ============================================================

saveRDS(
  Tcells,
  file.path(
    objects_dir,
    "fibroblasts.CD45.Tcells.annotated.rds"
  )
)

Tcells_annotation <- data.frame(
  barcode = colnames(Tcells),
  orig.ident = Tcells$orig.ident,
  seurat_clusters = Tcells$seurat_clusters,
  Cluster = Tcells$Cluster,
  stringsAsFactors = FALSE
)

write.table(
  Tcells_annotation,
  file = file.path(
    objects_dir,
    "fibroblasts.CD45.Tcells.annotation.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nAnotación final de T cells guardada correctamente.\n")