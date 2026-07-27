# ============================================================
# Functions.R
# Macs-Fibros.NSCLC
# ============================================================



# ============================================================
# SAMPLE TABLE
# ============================================================

get_sample_table <- function() {

  data.frame(
    sample_id = c(
      "m1", "m2", "m3", "m4",
      "m10", "m12", "m13", "m14",
      "m15", "m16", "m17", "m18"
    ),
    file_id = c(
      "m1.Bis", "m2.Bis", "m3.Bis", "m4.Bis",
      "m10", "m12", "m13", "m14",
      "m15", "m16", "m17", "m18"
    ),
    condition = c(
      rep("CCR2KO", 4),
      rep("WT", 4),
      rep("CD169DTR", 4)
    ),
    pcs = c(
      19, 16, 12, 18,
      18, 16, 20, 15,
      14, 24, 24, 12
    ),
    stringsAsFactors = FALSE
  )
}


# ============================================================
# READ SAMPLE
# No RenameCells here.
# DoubletFinder must see the original per-sample object.
# ============================================================

read_raw_sample <- function(
    datadir,
    condition,
    file_id
) {

  rds_file <- file.path(
    datadir,
    "raw",
    condition,
    paste0(file_id, ".rds")
  )

  if (!file.exists(rds_file)) {

    stop(
      "Input file not found: ",
      rds_file
    )
  }

  message("Reading: ", rds_file)

  obj <- readRDS(rds_file)

  if (!inherits(obj, "Seurat")) {

    stop(
      "Object is not Seurat: ",
      rds_file
    )
  }

  obj
}


# ============================================================
# LEGACY DOUBLET PREPROCESSING
# Exact order:
# MT -> QC -> Normalize -> HVG -> Scale -> PCA
# ============================================================

preprocess_doubletfinder <- function(obj) {

  DefaultAssay(obj) <- "RNA"

  obj$percent.mt <- PercentageFeatureSet(
    obj,
    pattern = "^mt-"
  )

  obj <- subset(
    obj,
    subset =
      nFeature_RNA > 300 &
      nFeature_RNA < 6000 &
      percent.mt < 20
  )

  obj <- NormalizeData(obj)

  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000
  )

  obj <- ScaleData(obj)

  obj <- RunPCA(obj)

  obj
}


# ============================================================
# LEGACY SAMPLE CLUSTERING
# ============================================================

prepare_doubletfinder_graph <- function(
    obj,
    pc
) {

  dims_use <- seq_len(pc)

  obj <- FindNeighbors(
    obj,
    dims = dims_use
  )

  obj <- FindClusters(
    obj,
    resolution = 0.5
  )

  obj <- RunUMAP(
    obj,
    dims = dims_use
  )

  obj
}


# ============================================================
# LEGACY DOUBLET RATE
# Reproduces cell() from old Functions.R
# ============================================================

legacy_doublet_rate <- function(n_cells) {

  if (
    n_cells > 1000 &&
    n_cells < 2000
  ) {

    return(0.015)
  }

  if (
    n_cells > 2000 &&
    n_cells < 3000
  ) {

    return(0.023)
  }

  if (
    n_cells > 3000 &&
    n_cells < 4000
  ) {

    return(0.03)
  }

  if (
    n_cells > 4000 &&
    n_cells < 5000
  ) {

    return(0.038)
  }

  if (
    n_cells > 3000 &&
    n_cells < 6000
  ) {

    return(0.046)
  }

  if (
    n_cells > 6000 &&
    n_cells < 7000
  ) {

    return(0.053)
  }

  if (
    n_cells > 7000 &&
    n_cells < 8000
  ) {

    return(0.061)
  }

  if (
    n_cells > 7000 &&
    n_cells < 9000
  ) {

    return(0.068)
  }

  if (n_cells > 10000) {

    return(0.08)
  }

  stop(
    "No legacy doublet rate defined for ",
    n_cells,
    " cells"
  )
}


# ============================================================
# pK
# ============================================================

select_doubletfinder_pk <- function(
    obj,
    pc
) {

  dims_use <- seq_len(pc)

  sweep.res.list <- DoubletFinder::paramSweep(
    obj,
    PCs = dims_use,
    sct = FALSE
  )

  sweep.stats <- DoubletFinder::summarizeSweep(
    sweep.res.list,
    GT = FALSE
  )

  bcmvn <- DoubletFinder::find.pK(
    sweep.stats
  )

  best_row <- which.max(
    bcmvn$BCmetric
  )[1]

  pk_max <- as.numeric(
    as.character(
      bcmvn$pK[best_row]
    )
  )

  list(
    pK = pk_max,
    bcmvn = bcmvn
  )
}


# ============================================================
# LEGACY DOUBLET FINDER
# ============================================================

run_doubletfinder <- function(
    obj,
    pc
) {

  obj <- preprocess_doubletfinder(
    obj
  )

  message(
    "Cells after QC: ",
    ncol(obj)
  )

  obj <- prepare_doubletfinder_graph(
    obj = obj,
    pc = pc
  )

  pk_result <- select_doubletfinder_pk(
    obj = obj,
    pc = pc
  )

  per <- legacy_doublet_rate(
    ncol(obj)
  )

  homotypic.prop <-
    DoubletFinder::modelHomotypic(
      obj$seurat_clusters
    )

  nExp_poi <- round(
    per * ncol(obj)
  )

  nExp_poi.adj <- round(
    nExp_poi *
      (1 - homotypic.prop)
  )

  message("PCs: 1:", pc)
  message("pK: ", pk_result$pK)
  message("Legacy doublet rate: ", per)
  message("nExp_poi: ", nExp_poi)
  message(
    "nExp_poi.adj: ",
    nExp_poi.adj,
    " [calculated but NOT used]"
  )

  # Historical workflow used nExp_poi.
  #
  # Legacy code used reuse.pANN = FALSE.
  # DoubletFinder 2.0.6 in this environment requires NULL
  # for the initial run.

  obj <- DoubletFinder::doubletFinder(
    obj,
    PCs = seq_len(pc),
    pN = 0.25,
    pK = pk_result$pK,
    nExp = nExp_poi,
    reuse.pANN = NULL,
    sct = FALSE
  )

  classification_col <- grep(
    "^DF.classifications",
    colnames(obj@meta.data),
    value = TRUE
  )

  classification_col <- tail(
    classification_col,
    1
  )

  if (
    length(classification_col) != 1
  ) {

    stop(
      "DoubletFinder classification column not found"
    )
  }

  doublet_cells <- rownames(
    obj@meta.data
  )[
    obj@meta.data[[classification_col]] ==
      "Doublet"
  ]

  list(
    object = obj,
    doublets = doublet_cells,
    bcmvn = pk_result$bcmvn,
    pK = pk_result$pK,
    rate = per,
    nExp = nExp_poi
  )
}


# ============================================================
# MATCH ORIGINAL DOUBLET BARCODES AFTER RenameCells
# ============================================================

match_sample_doublets <- function(
    obj,
    sample_id,
    doublet_barcodes
) {

  sample_cells <- colnames(obj)[
    obj$sample_id == sample_id
  ]

  original_names <- sub(
    paste0(
      "^",
      sample_id,
      "_"
    ),
    "",
    sample_cells
  )

  doublet_barcodes <- as.character(
    doublet_barcodes
  )

  matched <- sample_cells[
    original_names %in%
      doublet_barcodes
  ]

  if (
    length(matched) !=
    length(unique(doublet_barcodes))
  ) {

    stop(
      sample_id,
      ": expected ",
      length(unique(doublet_barcodes)),
      " doublets but matched ",
      length(matched)
    )
  }

  matched
}




PCA <- function(data){
## Selecting PCAs

# First metric 
#The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.

# Determine percent of variation associated with each PC
pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Second metric
# The point where the percent change in variation between the consequtive PCs is less than 0.1%.

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.

co2
list <- list(co1,co2)
return(list)

}
