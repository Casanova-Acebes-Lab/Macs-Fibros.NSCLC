# ============================================================
# 17_QC_batch_robustness.R
# Batch / sample-mixing QC and integration-robustness check
#
# Purpose (reviewer-proofing for the "why not integrate?" question):
#   1. Show there is no per-sample technical batch to correct
#      (replicates of the SAME genotype co-embed; per-cluster
#       sample entropy is high; no cluster is driven by one mouse).
#   2. Show the fibroblast/mesenchymal STATE STRUCTURE is preserved
#      after Harmony integration (ARI vs the original annotation),
#      i.e. the states are biology, not batch artifacts.
#   3. Export per-sample composition (proportions per mouse) so that
#      genotype effects on composition can be tested at the level of
#      biological replicates (not per cell).
#
# Input:
#   results/objects/fibroblasts.Not.Integrated.annotated.rds
#   (override with:  Rscript R/17_QC_batch_robustness.R /path/to/object.rds)
#
# Outputs:
#   results/plots/QC_batch/*.png|pdf
#   results/tables/QC_batch/*.tsv
#
# Env: same conda env as the rest of the repo (Seurat 5.4, harmony 1.2.4).
# No extra installs required.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(harmony)
})

set.seed(42)

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

config <- read.delim("config.tsv", stringsAsFactors = FALSE)
config <- setNames(config$value, config$key)
rsdir  <- config[["rsdir"]]

input_file <- if (length(args) >= 1L) {
  args[[1]]
} else {
  file.path(rsdir, "objects", "fibroblasts.Not.Integrated.annotated.rds")
}

if (!file.exists(input_file)) {
  stop("Input object not found: ", input_file,
       "\nPass the path as an argument if it lives elsewhere:",
       "\n  Rscript R/17_QC_batch_robustness.R /path/to/object.rds")
}

plot_dir  <- file.path(rsdir, "plots",  "QC_batch")
table_dir <- file.path(rsdir, "tables", "QC_batch")
dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load
# ------------------------------------------------------------

obj <- readRDS(input_file)
stopifnot(inherits(obj, "Seurat"))
DefaultAssay(obj) <- "RNA"

md <- obj@meta.data
cat("\nMetadata columns:\n"); print(colnames(md))

# ------------------------------------------------------------
# Detect the relevant metadata columns
# ------------------------------------------------------------

pick_col <- function(candidates, md) {
  hit <- candidates[candidates %in% colnames(md)]
  if (length(hit) >= 1L) hit[[1]] else NA_character_
}

genotype_col  <- pick_col(c("orig.ident", "genotype", "condition"), md)
condition_col <- pick_col(c("group", "Group", "state", "naive_tumor"), md)
sample_col    <- pick_col(c("sample_id", "tag", "Sample", "hash.ID", "orig.sample"), md)
state_col     <- pick_col(c("Clustering", "annotation", "cell_type"), md)
compartment_col <- pick_col(c("Clustering.Mesenchymal", "compartment"), md)

# Fallback: derive per-sample id from the "mX_" cell-name prefix used in 01_QC.R
if (is.na(sample_col)) {
  message("No explicit sample column found; deriving from cell-name prefix (mX_).")
  derived <- sub("_.*$", "", colnames(obj))          # e.g. "m15" from "m15_AAAC..."
  if (all(grepl("^m[0-9]+$", derived))) {
    obj$sample_derived <- derived
    sample_col <- "sample_derived"
  } else {
    stop("Could not identify a per-sample column. ",
         "Set `sample_col` manually near the top of the script.")
  }
}

cat("\nUsing columns:\n")
cat("  genotype   =", genotype_col,  "\n")
cat("  condition  =", condition_col, "\n")
cat("  sample     =", sample_col,    "\n")
cat("  state      =", state_col,     "\n")
cat("  compartment=", compartment_col, "\n\n")

stopifnot(!is.na(genotype_col), !is.na(sample_col), !is.na(state_col))
stopifnot("umap" %in% Reductions(obj), "pca" %in% Reductions(obj))

obj$.geno   <- as.character(md[[genotype_col]])
obj$.samp   <- as.character(obj@meta.data[[sample_col]])
obj$.cond   <- if (!is.na(condition_col)) as.character(md[[condition_col]]) else NA_character_
obj$.state  <- as.character(md[[state_col]])

cat("Cells per sample:\n"); print(table(obj$.samp))
cat("\nSample x genotype:\n"); print(table(obj$.samp, obj$.geno))

# ============================================================
# PART A. Visual mixing: UMAP by sample / genotype / condition
# ============================================================

save_umap <- function(group.by, file, cols = NULL, split.by = NULL) {
  p <- DimPlot(obj, reduction = "umap", group.by = group.by,
               split.by = split.by, shuffle = TRUE, seed = 42,
               raster = FALSE, pt.size = 0.2) +
    ggtitle(paste0("UMAP by ", group.by))
  ggsave(file.path(plot_dir, file), p, width = if (is.null(split.by)) 9 else 15,
         height = 7, dpi = 200, bg = "white")
}

save_umap(genotype_col, "UMAP_by_genotype.png")
save_umap(sample_col,   "UMAP_by_sample.png")
save_umap(sample_col,   "UMAP_by_sample_split_genotype.png", split.by = ".geno")
if (!is.na(condition_col)) save_umap(condition_col, "UMAP_by_condition.png")

# ============================================================
# PART B. Quantitative mixing: per-cluster sample entropy
#   Normalized Shannon entropy in [0,1]; 1 = perfectly mixed.
#   Computed (i) over all samples and (ii) over replicates
#   WITHIN each genotype (the clean, biology-free batch signal).
# ============================================================

norm_entropy <- function(counts) {
  p <- counts[counts > 0] / sum(counts)
  if (length(p) <= 1L) return(0)
  -sum(p * log(p)) / log(length(counts))   # normalize by max possible
}

# (i) global: mixing of all samples within each state
comp_global <- md %>%
  transmute(state = obj$.state, samp = obj$.samp) %>%
  count(state, samp, name = "n") %>%
  group_by(state) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

entropy_global <- comp_global %>%
  group_by(state) %>%
  summarise(
    n_cells      = sum(n),
    n_samples    = dplyr::n_distinct(samp),
    top_sample_pct = 100 * max(prop),
    entropy_norm = norm_entropy(n),
    .groups = "drop"
  ) %>%
  mutate(single_sample_dominated = top_sample_pct > 50) %>%
  arrange(entropy_norm)

write.table(entropy_global, file.path(table_dir, "per_state_sample_entropy_global.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(comp_global, file.path(table_dir, "per_state_sample_composition_global.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# (ii) within-genotype replicate mixing (removes the biology confound)
entropy_within <- md %>%
  transmute(state = obj$.state, samp = obj$.samp, geno = obj$.geno) %>%
  count(geno, state, samp, name = "n") %>%
  group_by(geno, state) %>%
  summarise(
    n_replicates = dplyr::n_distinct(samp),
    top_sample_pct = 100 * max(n) / sum(n),
    entropy_norm = norm_entropy(n),
    .groups = "drop"
  ) %>%
  arrange(geno, entropy_norm)

write.table(entropy_within, file.path(table_dir, "per_state_replicate_entropy_within_genotype.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nGlobal per-state sample entropy (0=one sample, 1=fully mixed):\n")
print(entropy_global)
cat("\nWithin-genotype replicate entropy (lowest = worst mixed):\n")
print(head(entropy_within, 15))

# ============================================================
# PART C. Integration robustness (Harmony)
#   Aggressively remove ALL per-sample variation (this also removes
#   genotype, since sample is nested in genotype), re-cluster, and
#   test whether the SAME states are recovered (ARI + best Jaccard).
#   Interpretation: high ARI => states are structure, not batch.
#   (We do NOT claim proportions are unchanged here.)
# ============================================================

adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  sum_comb <- function(x) sum(choose(x, 2))
  a_i <- rowSums(tab); b_j <- colSums(tab)
  index      <- sum_comb(tab)
  expected   <- sum_comb(a_i) * sum_comb(b_j) / choose(n, 2)
  maximum    <- 0.5 * (sum_comb(a_i) + sum_comb(b_j))
  (index - expected) / (maximum - expected)
}

pcs_use <- 1:24  # matches 01_QC.R

obj <- RunHarmony(
  obj,
  group.by.vars = sample_col,
  reduction.use = "pca",
  dims.use = pcs_use,
  reduction.save = "harmony",
  verbose = FALSE
)

obj <- FindNeighbors(obj, reduction = "harmony", dims = pcs_use, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.6, random.seed = 0, verbose = FALSE)
obj$harmony_clusters <- as.character(Idents(obj))

obj <- RunUMAP(obj, reduction = "harmony", dims = pcs_use,
               reduction.name = "umap_harmony", seed.use = 42, verbose = FALSE)

ari <- adjusted_rand_index(obj$.state, obj$harmony_clusters)
cat("\nAdjusted Rand Index (original annotation vs Harmony re-clustering): ",
    round(ari, 3), "\n")

# best-overlap (Jaccard) mapping of each original state to a Harmony cluster
cont <- as.data.frame.matrix(table(obj$.state, obj$harmony_clusters))
jacc <- data.frame(state = rownames(cont), stringsAsFactors = FALSE)
jacc$best_harmony_cluster <- colnames(cont)[apply(cont, 1, which.max)]
jacc$best_jaccard <- vapply(seq_len(nrow(cont)), function(i) {
  hc <- jacc$best_harmony_cluster[i]
  inter <- cont[i, hc]
  union <- sum(cont[i, ]) + sum(cont[, hc]) - inter
  inter / union
}, numeric(1))

write.table(data.frame(metric = "adjusted_rand_index", value = ari),
            file.path(table_dir, "harmony_robustness_ARI.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(jacc, file.path(table_dir, "harmony_state_best_jaccard.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cbind(state = rownames(cont), cont),
            file.path(table_dir, "harmony_vs_original_contingency.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p_h <- DimPlot(obj, reduction = "umap_harmony", group.by = ".state",
               shuffle = TRUE, seed = 42, raster = FALSE, pt.size = 0.2) +
  ggtitle(paste0("Harmony embedding, original labels (ARI = ", round(ari, 3), ")"))
ggsave(file.path(plot_dir, "UMAP_harmony_original_labels.png"),
       p_h, width = 10, height = 7, dpi = 200, bg = "white")

p_hs <- DimPlot(obj, reduction = "umap_harmony", group.by = sample_col,
                shuffle = TRUE, seed = 42, raster = FALSE, pt.size = 0.2) +
  ggtitle("Harmony embedding by sample")
ggsave(file.path(plot_dir, "UMAP_harmony_by_sample.png"),
       p_hs, width = 9, height = 7, dpi = 200, bg = "white")

# ============================================================
# PART D. Per-sample composition for replicate-level testing
#   Proportions of each state within each sample -> lets you test
#   genotype effects at the level of mice (addresses pseudoreplication).
# ============================================================

comp_sample <- md %>%
  transmute(samp = obj$.samp, geno = obj$.geno,
            cond = obj$.cond, state = obj$.state) %>%
  count(samp, geno, cond, state, name = "n") %>%
  group_by(samp) %>%
  mutate(total = sum(n), prop = n / total) %>%
  ungroup()

write.table(comp_sample, file.path(table_dir, "per_sample_state_proportions.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Optional quick test: per state, compare proportions across genotype
# (Kruskal-Wallis on per-sample proportions; within tumour if available).
logit <- function(p) { p <- pmin(pmax(p, 1e-4), 1 - 1e-4); log(p / (1 - p)) }

test_df <- comp_sample
if (!all(is.na(test_df$cond))) {
  test_df <- test_df %>% filter(cond %in% c("Tumor", "Tumour", "tumor"))
}

comp_tests <- test_df %>%
  group_by(state) %>%
  filter(dplyr::n_distinct(geno) >= 2, dplyr::n() >= 4) %>%
  summarise(
    kruskal_p = tryCatch(
      kruskal.test(logit(prop) ~ factor(geno))$p.value,
      error = function(e) NA_real_),
    n_samples = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(kruskal_fdr = p.adjust(kruskal_p, method = "fdr")) %>%
  arrange(kruskal_p)

write.table(comp_tests, file.path(table_dir, "per_sample_composition_kruskal.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nComposition test across genotype (per-sample proportions, tumour):\n")
print(comp_tests)

# ============================================================
# Session info
# ============================================================

writeLines(capture.output(sessionInfo()),
           file.path(table_dir, "sessionInfo_QC_batch.txt"))

cat("\n========================================\n")
cat("Batch-QC / robustness outputs written to:\n")
cat("  ", plot_dir, "\n")
cat("  ", table_dir, "\n")
cat("Key numbers to report: per-state entropy (Part B),",
    "Harmony ARI (Part C).\n")
cat("========================================\n")