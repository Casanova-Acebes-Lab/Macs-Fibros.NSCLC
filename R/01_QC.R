# ============================================================
# 01_QC.R
# QC, preprocessing and global clustering
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(cowplot)
  library(ggplot2)
  library(SingleR)
  library(celldex)
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

source("R/Functions.R")

dir.create(
  file.path(rsdir, "objects"),
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Historical highly variable genes
# ============================================================

hvg_file <- file.path(
  datadir,
  "HVG",
  "fibroblasts_historical_HVGs.rds"
)

if (!file.exists(hvg_file)) {
  stop("Historical HVG file not found: ", hvg_file)
}

hvg <- readRDS(hvg_file)

stopifnot(
  is.character(hvg),
  length(hvg) == 2000,
  length(unique(hvg)) == 2000
)

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
    "m1.Bis.rds",
    "m2.Bis.rds",
    "m3.Bis.rds",
    "m4.Bis.rds",
    "m10.rds",
    "m12.rds",
    "m13.rds",
    "m14.rds",
    "m15.rds",
    "m16.rds",
    "m17.rds",
    "m18.rds"
  ),
  tag = c(
    "m1-Tag9",
    "m2-Tag10",
    "m3-Tag11",
    "m4-Tag12",
    "m10-Tag2",
    "m12-Tag4",
    "m13-Tag1",
    "m14-Tag3",
    "m15-Tag7",
    "m16-Tag5",
    "m17-Tag6",
    "m18-Tag8"
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# Functions
# ============================================================

rename_like_old <- function(obj, sample_id) {

  if (sample_id == "m1") {
    return(obj)
  }

  suffix <- sub("^m", "", sample_id)

  new_names <- gsub(
    "1",
    suffix,
    colnames(obj)
  )

  RenameCells(
    object = obj,
    new.names = new_names
  )
}

collapse_to_single_leading_underscore <- function(obj) {

  new_names <- sub(
    "^_+",
    "",
    colnames(obj)
  )

  new_names <- paste0(
    "_",
    new_names
  )

  if (anyDuplicated(new_names)) {
    stop("Duplicated cell names after renaming.")
  }

  if (!identical(colnames(obj), new_names)) {
    obj <- RenameCells(
      object = obj,
      new.names = new_names
    )
  }

  obj
}

read_sample_remove_doublets <- function(
    sample_id,
    condition,
    raw_file,
    tag
) {

  message("Reading ", sample_id, " / ", condition)

  obj <- readRDS(
    file.path(
      datadir,
      "raw",
      condition,
      raw_file
    )
  )

  doublets <- readRDS(
    file.path(
      datadir,
      "barcodes.doublets",
      condition,
      paste0("barcodes.", sample_id, ".rds")
    )
  )

  obj <- subset(
    x = obj,
    cells = doublets,
    invert = TRUE
  )

  obj <- rename_like_old(
    obj = obj,
    sample_id = sample_id
  )

  obj$orig.ident <- condition
  obj$tag <- tag

  obj
}

# ============================================================
# Read samples and remove doublets
# ============================================================

object_list <- lapply(
  seq_len(nrow(samples)),
  function(i) {

    read_sample_remove_doublets(
      sample_id = samples$sample_id[i],
      condition = samples$condition[i],
      raw_file = samples$raw_file[i],
      tag = samples$tag[i]
    )
  }
)

names(object_list) <- samples$sample_id

# ============================================================
# Merge samples
# ============================================================

fibroblasts <- Reduce(
  f = function(x, y) {
    merge(
      x = x,
      y = y
    )
  },
  x = object_list
)

fibroblasts <- collapse_to_single_leading_underscore(
  fibroblasts
)

fibroblasts <- JoinLayers(
  object = fibroblasts,
  assay = "RNA"
)

# ============================================================
# Rebuild clean object
# ============================================================

counts_clean <- LayerData(
  object = fibroblasts,
  assay = "RNA",
  layer = "counts"
)

meta_clean <- fibroblasts@meta.data[
  colnames(counts_clean),
  ,
  drop = FALSE
]

fibroblasts <- CreateSeuratObject(
  counts = counts_clean,
  assay = "RNA"
)

fibroblasts <- AddMetaData(
  object = fibroblasts,
  metadata = meta_clean
)

fibroblasts <- collapse_to_single_leading_underscore(
  fibroblasts
)

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
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 6000 &
    percent.mt < 20
)

cat("\nObject after QC:\n")
print(fibroblasts)

cat("\nCells by condition:\n")
print(table(fibroblasts$orig.ident))

stopifnot(
  nrow(fibroblasts) == 19061,
  ncol(fibroblasts) == 51195,
  sum(fibroblasts$nCount_RNA) == 141223989,
  sum(fibroblasts$nFeature_RNA) == 72720248
)

expected_conditions <- c(
  CCR2KO = 24601,
  CD169DTR = 10582,
  WT = 16012
)

observed_conditions <- table(
  fibroblasts$orig.ident
)

stopifnot(
  identical(
    as.integer(observed_conditions[names(expected_conditions)]),
    as.integer(expected_conditions)
  )
)

# ============================================================
# Normalize
# ============================================================

fibroblasts <- NormalizeData(
  object = fibroblasts,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

VariableFeatures(fibroblasts) <- hvg

stopifnot(
  identical(
    VariableFeatures(fibroblasts),
    hvg
  )
)

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
# ============================================================
# PCA is recalculated from the historical set of 2,000 HVGs.
# Seurat automatically excludes features that cannot be used
# for PCA (for example, zero-variance features).
# ============================================================

fibroblasts <- RunPCA(
  object = fibroblasts,
  features = hvg,
  npcs = 50,
  seed.use = 42,
  verbose = TRUE
)

pca_used <- rownames(
  Loadings(
    fibroblasts,
    reduction = "pca"
  )
)

cat(
  "
Historical HVGs supplied:",
  length(hvg),
  "
"
)

cat(
  "Features actually used by PCA:",
  length(pca_used),
  "
"
)

stopifnot(
  length(pca_used) == 1801,
  all(pca_used %in% hvg)
)

# ============================================================
# Neighbors and clustering
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

cat("\nClusters:\n")
print(table(fibroblasts$seurat_clusters))

expected_clusters <- c(
  "0" = 8960,
  "1" = 3999,
  "2" = 3981,
  "3" = 3709,
  "4" = 3360,
  "5" = 3338,
  "6" = 3307,
  "7" = 3052,
  "8" = 2778,
  "9" = 2521,
  "10" = 2443,
  "11" = 1791,
  "12" = 1227,
  "13" = 1201,
  "14" = 920,
  "15" = 842,
  "16" = 790,
  "17" = 684,
  "18" = 468,
  "19" = 398,
  "20" = 382,
  "21" = 351,
  "22" = 322,
  "23" = 170,
  "24" = 118,
  "25" = 83
)

observed_clusters <- table(
  fibroblasts$seurat_clusters
)

stopifnot(
  identical(
    names(observed_clusters),
    names(expected_clusters)
  ),
  identical(
    as.integer(observed_clusters),
    as.integer(expected_clusters)
  )
)

# ============================================================
# UMAP recalculated in the current environment
# ============================================================

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
  reduction.name = "umap",
  reduction.key = "UMAP_",
  verbose = TRUE
)


# ============================================================
# SingleR annotation
# Used for broad lineage QC and removal of immune contaminants
# ============================================================

## SingleR Cell identification


fibroblasts <- JoinLayers(fibroblasts)


sce <- as.SingleCellExperiment(DietSeurat(fibroblasts))

ref <- celldex::ImmGenData()

ref.main <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.main)

ref.fine <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.fine)


fibroblasts@meta.data$Cell_type.Image <- ref.main$pruned.labels
fibroblasts@meta.data$Cell_type.Image.fine <- ref.fine$pruned.labels




fibroblasts@meta.data$Cell_type.Image <- gsub("Endothelial cells.*","Endothelial",
                            gsub("Epithelial cells.*","Epithelial", 
                            gsub("NK cells.*","NK",
                            gsub("Stromal cells.*","Stromal",
                            gsub("NA.*","Ind",
                            fibroblasts@meta.data$Cell_type.Image)))))






# ============================================================
# Save final QC object
# ============================================================

qc_object_file <- file.path(rsdir, "objects",
  "fibroblasts.QC.rds"
)

# Basic integrity checks before saving
stopifnot(
  ncol(fibroblasts) == nrow(fibroblasts@meta.data),
  identical(
    colnames(fibroblasts),
    rownames(fibroblasts@meta.data)
  ),
  "pca" %in% Reductions(fibroblasts),
  "umap" %in% Reductions(fibroblasts),
  "seurat_clusters" %in% colnames(fibroblasts@meta.data),
  "Cell_type.Image" %in% colnames(fibroblasts@meta.data),
  "Cell_type.Image.fine" %in% colnames(fibroblasts@meta.data)
)

saveRDS(
  object = fibroblasts,
  file = qc_object_file
)

if (!file.exists(qc_object_file)) {
  stop(
    "The final QC object could not be saved at: ",
    qc_object_file
  )
}

cat("\n========================================\n")
cat("Final QC object saved successfully\n")
cat("Path:", qc_object_file, "\n")
cat("Cells:", ncol(fibroblasts), "\n")
cat("Features:", nrow(fibroblasts), "\n")
cat(
  "File size:",
  round(file.info(qc_object_file)$size / 1024^3, 2),
  "GB\n"
)
cat("========================================\n")