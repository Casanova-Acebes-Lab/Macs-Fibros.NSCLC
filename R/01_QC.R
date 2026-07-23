# ============================================================
# QC / Preprocessing
# Fibroblasts object
# Raw counts + original doublet barcodes 
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(cowplot)
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

keys <- config$key
values <- config$value
names(values) <- keys

rsdir <- values["rsdir"]
datadir <- values["datadir"]
outdir <- values["outdir"]

dir.create(rsdir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "objects"), recursive = TRUE, showWarnings = FALSE)


source("R/Functions.R")

# ============================================================
# Historical HVGs
# ============================================================

hvg_file <- paste0(datadir,"/HVG/fibroblasts_historical_HVGs.rds")

if (!file.exists(hvg_file)) {
  stop("No encuentro el archivo de HVGs históricos: ", hvg_file)
}

historical_hvg <- readRDS(hvg_file)

cat("Historical HVGs loaded from:", hvg_file, "\n")
cat("Number of historical HVGs:", length(historical_hvg), "\n")

stopifnot(length(historical_hvg) == 2000)

# ============================================================
# Sample table
# ============================================================

samples <- data.frame(
  sample_id = c(
    "m1", "m2", "m3", "m4",
    "m10", "m12", "m13", "m14",
    "m15", "m16", "m17", "m18"
  ),
  condition = c(
    rep("CCR2KO", 4),
    rep("WT", 4),
    rep("CD169DTR", 4)
  ),
  raw_file = c(
    "m1.Bis.rds", "m2.Bis.rds", "m3.Bis.rds", "m4.Bis.rds",
    "m10.rds", "m12.rds", "m13.rds", "m14.rds",
    "m15.rds", "m16.rds", "m17.rds", "m18.rds"
  ),
  tag = c(
    "m1-Tag9", "m2-Tag10", "m3-Tag11", "m4-Tag12",
    "m10-Tag2", "m12-Tag4", "m13-Tag1", "m14-Tag3",
    "m15-Tag7", "m16-Tag5", "m17-Tag6", "m18-Tag8"
  ),
  stringsAsFactors = FALSE
)



# ============================================================
# Read raw samples and remove doublets
# ============================================================

object_list <- lapply(seq_len(nrow(samples)), function(i) {
  
  read_sample_remove_doublets(
    sample_id = samples$sample_id[i],
    condition = samples$condition[i],
    raw_file = samples$raw_file[i],
    tag = samples$tag[i]
  )
})

names(object_list) <- samples$sample_id

# ============================================================
# Merge all samples
# ============================================================

fibroblasts <- Reduce(
  f = function(x, y) merge(x = x, y = y),
  x = object_list
)

# Canonical old-style cell names: exactly one leading underscore
fibroblasts <- collapse_to_single_leading_underscore(fibroblasts)

# Join Seurat v5 layers
fibroblasts <- JoinLayers(
  object = fibroblasts,
  assay = "RNA"
)

# Rebuild clean object from joined counts
counts_clean <- LayerData(
  object = fibroblasts,
  assay = "RNA",
  layer = "counts"
)

meta_clean <- fibroblasts@meta.data[colnames(counts_clean), , drop = FALSE]

fibroblasts <- CreateSeuratObject(
  counts = counts_clean,
  assay = "RNA"
)

fibroblasts <- AddMetaData(
  object = fibroblasts,
  metadata = meta_clean
)

fibroblasts <- collapse_to_single_leading_underscore(fibroblasts)

cat("\nMerged object before QC:\n")
print(fibroblasts)

cat("\nCells by condition before QC:\n")
print(table(fibroblasts$orig.ident))

cat("\nCells by tag before QC:\n")
print(table(fibroblasts$tag))

# ============================================================
# QC
# ============================================================

DefaultAssay(fibroblasts) <- "RNA"

fibroblasts$percent.mt <- PercentageFeatureSet(
  object = fibroblasts,
  pattern = "^mt-"
)

fibroblasts <- subset(
  x = fibroblasts,
  subset = nFeature_RNA > 300 &
    nFeature_RNA < 6000 &
    percent.mt < 20
)

cat("\nObject after QC:\n")
print(fibroblasts)

cat("\nCells by condition after QC:\n")
print(table(fibroblasts$orig.ident))

cat("\nCells by tag after QC:\n")
print(table(fibroblasts$tag))

if (ncol(fibroblasts) != 51195) {
  warning(
    "Expected 51195 cells after QC, but found ",
    ncol(fibroblasts),
    ". Check doublet barcodes, QC thresholds or input files."
  )
}

if (!all(historical_hvg %in% rownames(fibroblasts))) {
  missing_hvg <- setdiff(historical_hvg, rownames(fibroblasts))
  stop(
    "Some historical HVGs are missing from the object:\n",
    paste(missing_hvg, collapse = ", ")
  )
}

# ============================================================
# Normalize
# ============================================================

fibroblasts <- NormalizeData(
  object = fibroblasts,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

# ============================================================
# Historical HVGs
# ============================================================
# To ensure exact reproducibility of the  clustering, we use the
# original HVGs calculated in the first instance and stored from the reference
# fibroblasts object.


VariableFeatures(fibroblasts) <- historical_hvg

cat("\nVariable features fixed to historical HVGs:\n")
cat("Number of HVGs:", length(VariableFeatures(fibroblasts)), "\n")

# ============================================================
# Scale
# ============================================================


fibroblasts <- ScaleData(
  object = fibroblasts,
  features = rownames(fibroblasts),
  verbose = TRUE
)

# ============================================================
# PCA
# Parameters:
# npcs = 50
# seed.use = 42
# ============================================================

fibroblasts <- RunPCA(
  object = fibroblasts,
  features = historical_hvg,
  npcs = 50,
  seed.use = 42,
  verbose = TRUE
)

# ============================================================
# Neighbors, clustering and UMAP
# Parameters:
# dims = 1:24
# resolution = 0.6
# clustering random.seed = 0
# UMAP seed.use = 42
# ============================================================

pcs_use <- 1:24

fibroblasts <- FindNeighbors(
  object = fibroblasts,
  reduction = "pca",
  dims = pcs_use,
  k.param = 20,
  prune.SNN = 0.06666667,
  nn.method = "annoy",
  n.trees = 50,
  annoy.metric = "euclidean",
  l2.norm = FALSE,
  verbose = TRUE
)

fibroblasts <- FindClusters(
  object = fibroblasts,
  graph.name = "RNA_snn",
  resolution = 0.6,
  algorithm = 1,
  n.start = 10,
  n.iter = 10,
  random.seed = 0,
  group.singletons = TRUE,
  verbose = TRUE
)

fibroblasts <- RunUMAP(
  object = fibroblasts,
  reduction = "pca",
  dims = pcs_use,
  umap.method = "uwot",
  metric = "cosine",
  n.neighbors = 30,
  min.dist = 0.3,
  spread = 1,
  seed.use = 42,
  verbose = TRUE
)

# ============================================================
# Save QC object
# ============================================================

saveRDS(
  fibroblasts,
  file = file.path(outdir, "objects/fibroblasts.QC.rds")
)


