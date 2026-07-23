# ============================================================
# 00_doublets.R
# DoubletFinder per sample
# Macs-Fibros.NSCLC
# Rewritten from the original DoubletFinder script
# ============================================================


# ============================================================
# Packages
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(dplyr)
  library(scCustomize)
  library(ggplot2)
  library(cowplot)
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
# Directories
# ============================================================

dir.create(rsdir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(outdir, "pk"), recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(rsdir, "barcodes", "CCR2KO"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(rsdir, "barcodes", "WT"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(rsdir, "barcodes", "CD169DTR"), recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(rsdir, "doubletfinder"), recursive = TRUE, showWarnings = FALSE)


# ============================================================
# Source functions
# ============================================================

source("R/Functions.R")


# ============================================================
# Sample metadata
# PCs are the manually selected values from the original script
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
  sample_file = c(
    "m1.Bis", "m2.Bis", "m3.Bis", "m4.Bis",
    "m10", "m12", "m13", "m14",
    "m15", "m16", "m17", "m18"
  ),
  pc = c(
    19, 16, 12, 18,
    18, 16, 20, 15,
    14, 24, 24, 12
  ),
  plot_label = c(
    "M1", "M2", "M3", "M4",
    "M10", "M12", "M13", "M14",
    "M15", "M16", "M17", "M18"
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# Check input files
# ============================================================

samples$input_file <- file.path(
  datadir,
  "raw",
  samples$condition,
  paste0(samples$sample_file, ".rds")
)

missing_files <- samples$input_file[!file.exists(samples$input_file)]

if (length(missing_files) > 0) {
  stop(
    "Missing input files:\n",
    paste(missing_files, collapse = "\n")
  )
}


# ============================================================
# Run DoubletFinder sample by sample
# ============================================================

plots <- list()
summary_list <- list()
pca_metrics_list <- list()

for (i in seq_len(nrow(samples))) {
  
  sample_id <- samples$sample_id[i]
  condition <- samples$condition[i]
  sample_file <- samples$sample_file[i]
  pc <- samples$pc[i]
  plot_label <- samples$plot_label[i]
  
  message("========================================")
  message("Processing sample: ", sample_id, " | ", condition)
  message("Using PCs: 1:", pc)
  message("========================================")
  
  # ----------------------------
  # Preprocess
  # ----------------------------
  
  obj <- preprocess(
    path = condition,
    sample = sample_file
  )
  
  # ----------------------------
  # PCA metrics
  # Manual PCs are used, as in the original script
  # ----------------------------
  
  pca_metrics <- PCA(obj)
  
  pca_metrics_list[[sample_id]] <- data.frame(
    sample_id = sample_id,
    condition = condition,
    co1 = pca_metrics[[1]],
    co2 = pca_metrics[[2]],
    manual_pc_used = pc,
    stringsAsFactors = FALSE
  )
  
  # ----------------------------
  # UMAP / clustering
  # ----------------------------
  
  obj <- umap(
    data = obj,
    pc = pc
  )
  
  # ----------------------------
  # DoubletFinder
  # ----------------------------
  
  obj <- pk(
    data = obj,
    pc = pc
  )
  
  # ----------------------------
  # Plot
  # ----------------------------
  
  p <- plot_doublets(
    data = obj,
    sample = plot_label
  )
  
  plots[[sample_id]] <- p
  
  png(
    filename = file.path(outdir, "pk", paste0(sample_id, ".png")),
    width = 800,
    height = 600
  )
  print(p)
  dev.off()
  
  # ----------------------------
  # Extract and save doublet barcodes
  # ----------------------------
  
  doublet_barcodes <- barcodes(obj)
  
  saveRDS(
    doublet_barcodes,
    file = file.path(
      rsdir,
      "barcodes",
      condition,
      paste0("barcodes.", sample_id, ".rds")
    )
  )
  
  # ----------------------------
  # Summary in memory only
  # ----------------------------
  
  df_class_col <- get_doubletfinder_class_column(obj)
  obj$multi <- obj@meta.data[[df_class_col]]
  
  tab <- table(obj$multi)
  
  n_doublets <- ifelse("Doublet" %in% names(tab), tab[["Doublet"]], 0)
  n_singlets <- ifelse("Singlet" %in% names(tab), tab[["Singlet"]], 0)
  
  summary_list[[sample_id]] <- data.frame(
    sample_id = sample_id,
    condition = condition,
    cells_after_qc = ncol(obj),
    manual_pc_used = pc,
    doublets = as.integer(n_doublets),
    singlets = as.integer(n_singlets),
    percent_doublets = round(
      100 * as.integer(n_doublets) /
        (as.integer(n_doublets) + as.integer(n_singlets)),
      4
    ),
    stringsAsFactors = FALSE
  )
  
  # ----------------------------
  # Save annotated object for inspection
  # ----------------------------
  
  saveRDS(
    obj,
    file = file.path(
      rsdir,
      "doubletfinder",
      paste0(sample_id, "_doubletfinder_annotated.rds")
    )
  )
  
  rm(obj, p, doublet_barcodes)
  gc()
}


# ============================================================
# Combined DoubletFinder plots
# ============================================================

png(
  filename = file.path(outdir, "pk", "all.png"),
  width = 1800,
  height = 1200
)

print(
  plot_grid(
    plots[["m1"]],
    plots[["m2"]],
    plots[["m3"]],
    plots[["m4"]],
    plots[["m10"]],
    plots[["m12"]],
    plots[["m13"]],
    plots[["m14"]],
    plots[["m15"]],
    plots[["m16"]],
    plots[["m17"]],
    plots[["m18"]],
    ncol = 4
  )
)

dev.off()


# ============================================================
# Summaries in memory only
# ============================================================

doubletfinder_summary <- dplyr::bind_rows(summary_list)
pca_metrics_summary <- dplyr::bind_rows(pca_metrics_list)

print(doubletfinder_summary)
print(pca_metrics_summary)


# ============================================================
# Final messages
# ============================================================

message("00_doublets.R completed successfully.")
message("Doublet barcodes saved in: ", file.path(rsdir, "barcodes"))
message("Annotated DoubletFinder objects saved in: ", file.path(rsdir, "doubletfinder"))
message("DoubletFinder plots saved in: ", file.path(outdir, "pk"))