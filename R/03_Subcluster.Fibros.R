# ============================================================
# 03_Subcluster.Fibros.R
# Historical fibroblast/mesothelial subclustering
#
# Starting object:
#   fibroblasts.QC.rds
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(cowplot)
  library(ggplot2)
  library(scCustomize)
  library(dplyr)
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

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(rsdir, "objects"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(outdir, "Fibros"),
  recursive = TRUE,
  showWarnings = FALSE
)

source("R/Functions.R")

# ============================================================
# Read QC object
# ============================================================

fibroblasts <- readRDS(
  file.path(
    rsdir,
    "objects",
    "fibroblasts.QC.rds"
  )
)

DefaultAssay(fibroblasts) <- "RNA"
Idents(fibroblasts) <- "seurat_clusters"

# ============================================================
# Extract CD45-negative compartment
#
# As in the historical script:
# - remove CD45-positive cells
# - remove global macrophage clusters 15 and 21
# ============================================================

fibros.CD45N <- subset(
  fibroblasts,
  subset = (
    CD == "CD45+" |
      seurat_clusters == 15 |
      seurat_clusters == 21
  ),
  invert = TRUE
)



fibros.CD45N <- subset(fibros.CD45N, subset=Cell_type.Image %in% c("B cells", "B cells, pro", "Basophils"
, "DC", "Eosinophils", "ILC", "Macrophages",  "Mast cells", "Monocytes", "Neutrophils",
"NKT",  'T cells', "Tgd", "NK"), invert = TRUE)


# ============================================================
# Extract global fibroblast and mesothelial clusters
# ============================================================

fibros.CD45N <- subset(
  fibros.CD45N,
  subset = (
    seurat_clusters == 17 |
      seurat_clusters == 20
  )
)

fibros.CD45N$global_cluster <- as.character(
  fibros.CD45N$seurat_clusters
)

cat(
  "\nFibroblast/mesothelial cells selected:",
  ncol(fibros.CD45N),
  "\n"
)

print(
  table(
    fibros.CD45N$global_cluster
  )
)



# ============================================================
# First historical clustering
#
# scale():
#   ScaleData using all genes
#   RunPCA using existing VariableFeatures
#
# Clustering:
#   PCA dimensions 1:10
#   resolution 0.2
# ============================================================

fibros.CD45N <- scale(
  fibros.CD45N
)

fibros.CD45N.PCA <- PCA(
  fibros.CD45N
)

cat("\nFirst PCA criteria:\n")
print(fibros.CD45N.PCA[[1]])
print(fibros.CD45N.PCA[[2]])

fibros.CD45N <- FindNeighbors(
  fibros.CD45N,
  dims = 1:10
)

fibros.CD45N <- FindClusters(
  fibros.CD45N,
  verbose = FALSE,
  resolution = 0.2
)

fibros.CD45N <- RunUMAP(
  fibros.CD45N,
  dims = 1:10
)

fibros.CD45N$fibro_first_cluster <- as.character(
  fibros.CD45N$seurat_clusters
)

cat("\nFirst fibroblast clustering:\n")
print(table(fibros.CD45N$fibro_first_cluster))


png(paste0(outdir,"/Umap.fibroblasts.CD45.Fibros.png"), width=1200, height=800)
p <- DimPlot(
  fibros.CD45N, reduction="umap", pt.size = 2) 
p
dev.off()



# ============================================================
# FeaturePlot: fibroblast and mesothelial markers
# ============================================================

DefaultAssay(fibros.CD45N) <- "RNA"

markers_fibro_meso <- c(
  # Fibroblastos
  "Pdgfra",
  "Col1a1",
  "Col1a2",
  "Col3a1",
  "Dcn",
  "Lum",
  "Pi16",
  "Dpt",
  "Hhip",
  "Grem2",
  "Ddit4",
  "Postn",

  # Mesoteliales
  "Msln",
  "Upk3b",
  "Wt1",
  "Lrrn4",
  "Krt8",
  "Krt18",
  "Krt19",
  "Pdpn"
)

markers_fibro_meso <- intersect(
  markers_fibro_meso,
  rownames(fibros.CD45N)
)

p_feature_fibro_meso <- FeaturePlot_scCustom(
  fibros.CD45N,
  features = markers_fibro_meso,
  reduction = "umap",
  colors_use = c("lightgrey", "blue"),
  pt.size = 0.5
)

print(p_feature_fibro_meso)

ggsave(
  filename = file.path(
    outdir,
    "Fibros",
    "FeaturePlot_fibroblast_mesothelial_markers.pdf"
  ),
  plot = p_feature_fibro_meso,
  width = 16,
  height = 13,
  limitsize = FALSE
)


# ============================================================
# FindAllMarkers: fibroblast/mesothelial populations
# ============================================================

DefaultAssay(fibros.CD45N) <- "RNA"
Idents(fibros.CD45N) <- "seurat_clusters"

markers_fibros <- FindAllMarkers(
  fibros.CD45N,
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

markers_fibros <- markers_fibros[order(
  markers_fibros$cluster,
  -markers_fibros$avg_log2FC
), ]

head(markers_fibros)

write.table(
  markers_fibros,
  file.path(
    outdir,
    "Fibros",
    "FindAllMarkers_fibroblast_mesothelial.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Top 50 markers por población
top50_markers_fibros <- markers_fibros %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(
    order_by = avg_log2FC,
    n = 50,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

write.table(
  top50_markers_fibros,
  file.path(
    outdir,
    "Fibros",
    "FindAllMarkers_fibroblast_mesothelial_top50.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(
  table(markers_fibros$cluster)
)

# ============================================================
# Historical filtering
#
# Keep first-round clusters:
#   0, 2 and 4
# ============================================================

fibros.CD45N <- subset(
  fibros.CD45N,
  subset = (
    seurat_clusters == 0 |
      seurat_clusters == 2 |
      seurat_clusters == 4
  )
)

cat(
  "\nCells retained after first clustering:",
  ncol(fibros.CD45N),
  "\n"
)

print(
  table(
    fibros.CD45N$fibro_first_cluster
  )
)

# ============================================================
# Second historical clustering
#
# Historical order:
#   scale
#   FindNeighbors using VariableFeatures
#   FindClusters resolution 0.3
#   PCA criteria
#   UMAP dimensions 1:10
# ============================================================

fibros.CD45N <- scale(
  fibros.CD45N
)

fibros.CD45N <- FindNeighbors(
  fibros.CD45N,
  features = VariableFeatures(
    object = fibros.CD45N
  )
)

fibros.CD45N <- FindClusters(
  fibros.CD45N,
  verbose = FALSE,
  resolution = 0.3
)

fibros.CD45N.PCA <- PCA(
  fibros.CD45N
)

cat("\nSecond PCA criteria:\n")
print(fibros.CD45N.PCA[[1]])
print(fibros.CD45N.PCA[[2]])

fibros.CD45N <- RunUMAP(
  fibros.CD45N,
  dims = 1:10
)

fibros.CD45N$fibro_cluster <- as.character(
  fibros.CD45N$seurat_clusters
)

cat("\nFinal fibroblast clusters:\n")
print(table(fibros.CD45N$fibro_cluster))



png(paste0(outdir,"/Umap.fibroblasts.CD45.Fibros2.png"), width=1200, height=800)
p <- DimPlot(
  fibros.CD45N, reduction="umap", pt.size = 2) 
p
dev.off()


# ============================================================
# Historical manual annotation
# ============================================================

fibros.CD45N$fibros <- NA_character_

fibros.CD45N$fibros[
  fibros.CD45N$seurat_clusters == 0
] <- "Fibroblasts Universal"

fibros.CD45N$fibros[
  fibros.CD45N$seurat_clusters == 1
] <- "Mesothelial"

fibros.CD45N$fibros[
  fibros.CD45N$seurat_clusters == 2
] <- "Mesothelial"

fibros.CD45N$fibros[
  fibros.CD45N$seurat_clusters == 3
] <- "Fibroblasts Activated 1"

fibros.CD45N$fibros[
  fibros.CD45N$seurat_clusters == 4
] <- "Fibroblasts Activated 2"

fibros.CD45N$fibro_annotation <- factor(
  fibros.CD45N$fibros,
  levels = c(
    "Fibroblasts Universal",
    "Fibroblasts Activated 1",
    "Fibroblasts Activated 2",
    "Mesothelial"
  )
)

fibros.CD45N$orig.ident <- factor(
  as.character(fibros.CD45N$orig.ident),
  levels = c(
    "WT",
    "CD169DTR",
    "CCR2KO"
  )
)

Idents(fibros.CD45N) <- "fibro_annotation"

cat("\nFinal fibroblast annotation:\n")
print(table(fibros.CD45N$fibro_annotation))

cat("\nAnnotation by condition:\n")
print(
  table(
    fibros.CD45N$orig.ident,
    fibros.CD45N$fibro_annotation
  )
)

# ============================================================
# Final check plot
# ============================================================

Nora <- c(
  "Mesothelial" = "#FFD9B3",
  "Fibroblasts Universal" = "#F4A825",
  "Fibroblasts Activated 1" = "#B71C1C",
  "Fibroblasts Activated 2" = "#6D1A1A"
)

p1 <- DimPlot(
  fibros.CD45N,
  reduction = "umap",
  group.by = "fibro_annotation",
  cols = Nora,
  label = TRUE,
  repel = TRUE,
  pt.size = 1.4
) +
  NoLegend() +
  ggtitle("Fibroblast annotation")

p2 <- DimPlot(
  fibros.CD45N,
  reduction = "umap",
  group.by = "fibro_annotation",
  split.by = "orig.ident",
  cols = Nora,
  pt.size = 1
) +
  ggtitle("Fibroblast annotation by condition")

p_final <- plot_grid(
  p1,
  p2,
  ncol = 1,
  rel_heights = c(1, 1)
)

png(
  filename = file.path(
    outdir,
    "Fibros",
    "Fibroblasts_final_annotation.png"
  ),
  width = 2400,
  height = 1600,
  res = 180
)

print(p_final)

dev.off()

# ============================================================
# Save final object
# ============================================================

saveRDS(
  fibros.CD45N,
  file.path(
    rsdir,
    "objects",
    "fibroblasts.Fibros.annotated.rds"
  )
)

cat(
  "\nObject saved in:\n",
  file.path(
    rsdir,
    "objects",
    "fibroblasts.Fibros.annotated.rds"
  ),
  "\n"
)