# ============================================================
# 12_Figure4_DFK.R
#
# Figure 4D, 4F and 4K
#
# Ejecutar desde la raíz del repositorio.
# Figure 4K se reconstruye desde counts y metadata/QC; no utiliza
# resultados DESeq2 o GSEA precomputados.
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(pheatmap)
})

set.seed(123)

# ============================================================
# Configuración del repositorio
# ============================================================

if (!file.exists("config.tsv")) {
  stop(
    "No se encuentra config.tsv. Ejecuta el script desde la raíz ",
    "del repositorio Macs-Fibros.NSCLC."
  )
}

config <- read.delim(
  "config.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("key", "value") %in% colnames(config))) {
  stop("config.tsv debe contener las columnas 'key' y 'value'.")
}

config <- setNames(config$value, config$key)
required_config_keys <- c("rsdir", "datadir", "outdir")
missing_config_keys <- setdiff(required_config_keys, names(config))

if (length(missing_config_keys) > 0L) {
  stop(
    "Faltan claves en config.tsv: ",
    paste(missing_config_keys, collapse = ", ")
  )
}

rsdir <- unname(config[["rsdir"]])
datadir <- unname(config[["datadir"]])
outdir <- unname(config[["outdir"]])

table_dir <- file.path(rsdir, "tables", "Figure4")
bulk_table_dir <- file.path(table_dir, "BulkRNAseq")
figure_dir <- file.path(outdir, "Figure4")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bulk_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(data, filename) {
  write.table(
    data,
    file = filename,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# ============================================================
# Figure 4D - Axl FeaturePlot
# ============================================================

# ============================================================
# Single-cell input and validation for Figure 4D and Figure 4F
# ============================================================

if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Figures 4D and 4F require the installed package Seurat.")
}
if (!requireNamespace("scCustomize", quietly = TRUE)) {
  stop("Figures 4D and 4F require the installed package scCustomize.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Figures 4D and 4F require the installed package ggplot2.")
}

seurat_file <- file.path(
  rsdir,
  "objects",
  "fibroblasts.Not.Integrated.annotated.rds"
)

if (!file.exists(seurat_file)) {
  stop("No se encuentra el objeto Seurat final: ", seurat_file)
}

data <- readRDS(seurat_file)

if (!"RNA" %in% SeuratObject::Assays(data)) {
  stop("El objeto Seurat final no contiene el assay RNA.")
}
if (!"data" %in% SeuratObject::Layers(data[["RNA"]])) {
  stop("El objeto Seurat final no contiene la layer normalizada RNA/data.")
}
if (!"Clustering" %in% colnames(data[[]])) {
  stop("El objeto Seurat final no contiene la anotación 'Clustering'.")
}

required_single_cell_features <- c("Axl", "Pros1")
missing_single_cell_features <- setdiff(
  required_single_cell_features,
  rownames(data[["RNA"]])
)
if (length(missing_single_cell_features) > 0L) {
  stop(
    "Faltan genes requeridos en el assay RNA: ",
    paste(missing_single_cell_features, collapse = ", ")
  )
}

global_umap_reduction <- "umap"
if (!global_umap_reduction %in% SeuratObject::Reductions(data)) {
  stop(
    "No se encuentra la reducción UMAP global histórica: ",
    global_umap_reduction
  )
}

if (!exists(
  "FeaturePlot_scCustom",
  where = asNamespace("scCustomize"),
  inherits = FALSE
)) {
  stop(
    "scCustomize::FeaturePlot_scCustom no está disponible en la versión ",
    "instalada. Figure 4D no se generará con otra función."
  )
}

SeuratObject::DefaultAssay(data) <- "RNA"

umap_coordinates <- SeuratObject::Embeddings(
  data,
  reduction = global_umap_reduction
)
if (ncol(umap_coordinates) < 2L) {
  stop("La reducción UMAP global debe contener al menos dos dimensiones.")
}
if (!setequal(rownames(umap_coordinates), colnames(data))) {
  stop("Las células del UMAP global no coinciden con las del objeto final.")
}
umap_coordinates <- umap_coordinates[colnames(data), , drop = FALSE]

rna_data <- SeuratObject::LayerData(
  data[["RNA"]],
  layer = "data"
)
axl_expression <- as.numeric(rna_data["Axl", colnames(data)])
names(axl_expression) <- colnames(data)

Figure4D_source_data <- data.frame(
  cell = colnames(data),
  UMAP_1 = umap_coordinates[colnames(data), 1],
  UMAP_2 = umap_coordinates[colnames(data), 2],
  Axl = axl_expression[colnames(data)],
  row.names = NULL,
  check.names = FALSE
)

if (anyNA(Figure4D_source_data)) {
  stop("Figure 4D: los datos fuente contienen valores NA.")
}

write_tsv(
  Figure4D_source_data,
  file.path(table_dir, "Figure4D_source_data.tsv")
)

Figure4D <- scCustomize::FeaturePlot_scCustom(
  seurat_object = data,
  features = "Axl",
  order = TRUE,
  reduction = global_umap_reduction,
  layer = "data",
  label = FALSE
) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    panel.background = ggplot2::element_rect(
      fill = "white",
      color = NA
    ),
    plot.background = ggplot2::element_rect(
      fill = "white",
      color = NA
    )
  )

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "Figure4D_Axl_FeaturePlot.pdf"
  ),
  plot = Figure4D,
  width = 6,
  height = 5.5,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

# ============================================================
# Figure 4F - Pros1 violin plot
# ============================================================

if (!exists(
  "VlnPlot_scCustom",
  where = asNamespace("scCustomize"),
  inherits = FALSE
)) {
  stop(
    "scCustomize::VlnPlot_scCustom no está disponible en la versión ",
    "instalada."
  )
}

# Exact values present in data$Clustering, in the required biological order.
figure4F_clustering_order <- c(
  "MoMacs",
  "em-AMs",
  "mo-AMs",
  "Endothelial Vascular",
  "Aerocytes",
  "Prolif. Cancer Cells",
  "EMT Cancer Cells",
  "Fibroblasts Universal",
  "Fibroblasts Activated 1",
  "Fibroblasts Activated 2"
)

present_clustering <- sort(unique(as.character(data$Clustering)))
message(
  "Valores presentes en data$Clustering:\n",
  paste(present_clustering, collapse = "\n")
)
missing_figure4F_populations <- setdiff(
  figure4F_clustering_order,
  present_clustering
)
if (length(missing_figure4F_populations) > 0L) {
  stop(
    "Figure 4F: no se identificaron inequívocamente estas poblaciones: ",
    paste(missing_figure4F_populations, collapse = ", ")
  )
}

# Presentation labels are local to Figure 4F and preserve the exact
# em-AMs and mo-AMs annotations present in data$Clustering.
figure4F_display_labels <- c(
  "MoMacs" = "MoMacs",
  "em-AMs" = "em-AMs",
  "mo-AMs" = "mo-AMs",
  "Endothelial Vascular" = "Endothelial Vascular",
  "Aerocytes" = "Aerocytes",
  "Prolif. Cancer Cells" = "Prolif. Cancer Cells",
  "EMT Cancer Cells" = "EMT Cancer Cells",
  "Fibroblasts Universal" = "Fibroblasts Universal",
  "Fibroblasts Activated 1" = "Fibroblasts Activated 1",
  "Fibroblasts Activated 2" = "Fibroblasts Activated 2"
)
figure4F_celltype_order <- unname(
  figure4F_display_labels[figure4F_clustering_order]
)

# Fixed Figure 4F palette supplied for these presentation labels.
Figure4F_colors <- c(
  "MoMacs" = "#2A9D8F",
  "em-AMs" = "#1B4332",
  "mo-AMs" = "#52B788",
  "Endothelial Vascular" = "#4A90E2",
  "Aerocytes" = "#1B4F72",
  "Prolif. Cancer Cells" = "#a8bfce",
  "EMT Cancer Cells" = "#bb6666",
  "Fibroblasts Universal" = "#FBCDC7",
  "Fibroblasts Activated 1" = "#F94144",
  "Fibroblasts Activated 2" = "#6D1A1A"
)

clustering_character <- as.character(data$Clustering)
figure4F_cells <- colnames(data)[
  clustering_character %in% figure4F_clustering_order
]
if (length(figure4F_cells) == 0L) {
  stop("Figure 4F: ninguna célula pertenece a las poblaciones seleccionadas.")
}

figure4F_labels_all <- unname(
  figure4F_display_labels[clustering_character]
)
data$Figure4F_celltype <- factor(
  figure4F_labels_all,
  levels = figure4F_celltype_order
)

missing_figure4F_colors <- setdiff(
  levels(data$Figure4F_celltype),
  names(Figure4F_colors)
)
unexpected_figure4F_colors <- setdiff(
  names(Figure4F_colors),
  levels(data$Figure4F_celltype)
)
if (
  length(missing_figure4F_colors) > 0L ||
    length(unexpected_figure4F_colors) > 0L
) {
  stop(
    "Figure 4F: los niveles y la paleta no coinciden. Colores ausentes: ",
    paste(missing_figure4F_colors, collapse = ", "),
    "; categorías inesperadas en la paleta: ",
    paste(unexpected_figure4F_colors, collapse = ", ")
  )
}
if (!identical(
  levels(data$Figure4F_celltype),
  names(Figure4F_colors)
)) {
  stop(
    "Figure 4F: el orden de Figure4F_colors no coincide con el orden ",
    "del factor Figure4F_celltype."
  )
}

pros1_expression <- as.numeric(rna_data["Pros1", figure4F_cells])
names(pros1_expression) <- figure4F_cells

Figure4F_source_data <- data.frame(
  cell = figure4F_cells,
  Clustering = clustering_character[match(figure4F_cells, colnames(data))],
  Figure4F_celltype = as.character(
    data$Figure4F_celltype[match(figure4F_cells, colnames(data))]
  ),
  Pros1 = pros1_expression[figure4F_cells],
  row.names = NULL,
  check.names = FALSE
)

if (anyNA(Figure4F_source_data)) {
  stop("Figure 4F: los datos fuente contienen valores NA.")
}

write_tsv(
  Figure4F_source_data,
  file.path(table_dir, "Figure4F_source_data.tsv")
)

data_figure4F <- subset(data, cells = figure4F_cells)

Figure4F <- scCustomize::VlnPlot_scCustom(
  seurat_object = data_figure4F,
  features = "Pros1",
  colors_use = Figure4F_colors,
  pt.size = 0,
  group.by = "Figure4F_celltype",
  plot_boxplot = FALSE,
  add.noise = FALSE,
  layer = "data"
) +
  Seurat::NoLegend() +
  ggplot2::labs(
    x = NULL,
    y = "Expression level",
    title = "Pros1"
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    plot.background = ggplot2::element_rect(
      fill = "white",
      color = NA
    )
  )

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "Figure4F_Pros1_Violin.pdf"
  ),
  plot = Figure4F,
  width = 11,
  height = 4.5,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

rm(rna_data)

# ============================================================
# Figure 4K - Bulk RNA-seq pathway heatmap
# ============================================================

# ============================================================
# Inputs
# ============================================================

counts_file <- file.path("data", "Figure4K", "matrix.tsv")

metadata_file <- file.path(
  "data",
  "Figure4K",
  "metadata_with_qc.tsv"
)

missing_input_files <- c(counts_file, metadata_file)[
  !file.exists(c(counts_file, metadata_file))
]

if (length(missing_input_files) > 0L) {
  stop(
    "No se encuentran los archivos de input: ",
    paste(missing_input_files, collapse = ", ")
  )
}

# ============================================================
# Lectura histórica de counts
# ============================================================

counts <- read.table(
  counts_file,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

if ("sample" %in% colnames(counts)) {
  colnames(counts)[colnames(counts) == "sample"] <- "gene"
}

if (!"gene" %in% colnames(counts)) {
  stop(
    "matrix.tsv no contiene la columna histórica 'sample' ",
    "(renombrada a 'gene')."
  )
}

if (anyNA(counts$gene) || any(counts$gene == "")) {
  stop("La columna de símbolos génicos contiene NA o nombres vacíos.")
}

duplicated_genes <- duplicated(counts$gene)
message(
  "Genes antes de eliminar símbolos duplicados: ",
  nrow(counts),
  "; genes conservados: ",
  sum(!duplicated_genes)
)
counts <- counts[!duplicated_genes, , drop = FALSE]
rownames(counts) <- counts$gene
counts$gene <- NULL

numeric_columns <- vapply(counts, is.numeric, logical(1))
if (!all(numeric_columns)) {
  stop(
    "Las siguientes columnas de counts no son numéricas: ",
    paste(names(numeric_columns)[!numeric_columns], collapse = ", ")
  )
}

counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"

if (anyNA(counts)) {
  stop("Se introdujeron NA al convertir los counts a matriz numérica.")
}
if (any(!is.finite(counts)) || any(counts < 0)) {
  stop("Los counts deben ser valores finitos y no negativos.")
}
if (any(counts != round(counts))) {
  stop("DESeq2 requiere counts enteros, pero matrix.tsv contiene decimales.")
}
storage.mode(counts) <- "integer"

# ============================================================
# Metadata y correspondencia de muestras
# ============================================================

metadata <- read.table(
  metadata_file,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

required_metadata_columns <- c(
  "condition",
  "condition.run",
  "star_uniquely_mapped_percent",
  "umi_mean_umi_per_pos"
)
missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(metadata)
)

if (length(missing_metadata_columns) > 0L) {
  stop(
    "Faltan columnas obligatorias en metadata: ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

missing_count_samples <- setdiff(colnames(counts), rownames(metadata))
if (length(missing_count_samples) > 0L) {
  stop(
    "Hay muestras de counts ausentes en metadata: ",
    paste(missing_count_samples, collapse = ", ")
  )
}

metadata <- metadata[colnames(counts), , drop = FALSE]
if (!all(colnames(counts) == rownames(metadata))) {
  stop("No se pudo alinear metadata con las columnas de counts.")
}

# ============================================================
# Technical correction between sequencing runs
# ============================================================
#
# This RNA-seq experiment was generated in two sequencing runs.
#
# The original run contained Control, CM, Bencetinib and CM + Bencetinib
# samples. Technical QC subsequently showed a systematic difference between
# the original Control samples and the repeated Control samples.
#
# Control and CM samples were therefore repeated in a second run.
#
# For the CM vs Control comparison, only the repeated Control and CM samples
# from the same run are used. Because both biological groups were generated
# together in the rerun, no technical covariate is included in that DESeq2
# design.
#
# Bencetinib and CM + Bencetinib samples, however, originate from the original
# run, whereas the valid Control samples used for comparison originate from
# the rerun. These comparisons therefore span the technical run difference.
#
# To account for this measured technical variation, a continuous technical
# covariate (`qcPC1_ctrl`) is reconstructed from QC metrics using the six
# biological Control samples available across the two runs.
#
# The technical axis is learned only from Controls so that it captures the
# reproducible technical difference between runs without using treatment
# identity to define the correction.
#
# The two QC variables used historically are:
# - star_uniquely_mapped_percent
# - umi_mean_umi_per_pos
#
# PC1 of these scaled QC metrics is then projected onto all samples and used
# as a covariate in the DESeq2 design for:
#
# Bencetinib vs Control
# CM + Bencetinib vs Control
#
# IMPORTANT:
# The raw count matrix itself is NOT batch-corrected or transformed to remove
# the run effect. Raw counts are supplied directly to DESeq2. The measured
# technical variation is modeled statistically using the historical design:
#
# ~ qcPC1_ctrl + condition_simple
#
# `qcPC1_ctrl` should not be interpreted as a perfect batch correction or as
# eliminating every possible source of confounding between run and treatment.
# This distinction should be preserved for reproducibility and interpretation.

# ============================================================
# Reconstrucción histórica de qcPC1_ctrl
# ============================================================
#
# The three original Controls are retained temporarily at this stage because,
# together with the three rerun Controls, they are required to learn the
# historical technical axis. The original Controls are removed only afterwards
# from the downstream differential-expression analysis.

qc_variables <- c(
  "star_uniquely_mapped_percent",
  "umi_mean_umi_per_pos"
)

technical_controls <- c(
  "1-Fibros-control-n3",
  "2-Fibros-control-n3",
  "3-Fibros-control-n3",
  "01Fib-cntrln3",
  "02Fib-cntrln3",
  "03Fib-cntrln3"
)

missing_technical_controls <- setdiff(
  technical_controls,
  rownames(metadata)
)
if (length(missing_technical_controls) > 0L) {
  stop(
    "Faltan controles para reconstruir qcPC1_ctrl: ",
    paste(missing_technical_controls, collapse = ", ")
  )
}

qc_ctrl <- metadata[technical_controls, qc_variables, drop = FALSE]
qc_ctrl[] <- lapply(qc_ctrl, as.numeric)
if (anyNA(qc_ctrl)) {
  stop("Hay NA en las métricas técnicas de los seis controles.")
}

qc_center <- vapply(qc_ctrl, mean, numeric(1))
qc_scale <- vapply(qc_ctrl, stats::sd, numeric(1))
if (any(!is.finite(qc_scale)) || any(qc_scale <= 0)) {
  stop("Las SD de las dos métricas técnicas deben ser mayores que cero.")
}

qc_ctrl_scaled <- sweep(
  sweep(as.matrix(qc_ctrl), 2, qc_center, FUN = "-"),
  2,
  qc_scale,
  FUN = "/"
)
qc_ctrl_pca <- prcomp(qc_ctrl_scaled, scale. = FALSE)

X_all <- metadata[, qc_variables, drop = FALSE]
X_all[] <- lapply(X_all, as.numeric)
if (anyNA(X_all)) {
  stop("Hay NA en las métricas técnicas de las muestras completas.")
}
X_all_scaled <- sweep(
  sweep(as.matrix(X_all), 2, qc_center, FUN = "-"),
  2,
  qc_scale,
  FUN = "/"
)
metadata$qcPC1_ctrl <- as.numeric(
  X_all_scaled %*% qc_ctrl_pca$rotation[, 1]
)

qc_audit <- data.frame(
  sample = rownames(metadata),
  condition = metadata$condition,
  condition.run = metadata$condition.run,
  star_uniquely_mapped_percent = metadata$star_uniquely_mapped_percent,
  umi_mean_umi_per_pos = metadata$umi_mean_umi_per_pos,
  qcPC1_ctrl = metadata$qcPC1_ctrl,
  check.names = FALSE
)
write_tsv(
  qc_audit,
  file.path(bulk_table_dir, "Figure4K_qcPC1_ctrl.tsv")
)

# ============================================================
# Limpieza histórica de muestras
# ============================================================

samples_to_remove <- c(
  "1-Fibros-control-n3",
  "2-Fibros-control-n3",
  "3-Fibros-control-n3",
  "9-Fibros-AM-CM-n3",
  "10-Fibros-AM-CM-n3",
  "11-Fibros-AM-CM-n3",
  "12Fib-AMCMn3"
)

missing_samples_to_remove <- setdiff(
  samples_to_remove,
  rownames(metadata)
)
if (length(missing_samples_to_remove) > 0L) {
  stop(
    "Faltan muestras históricas que debían eliminarse: ",
    paste(missing_samples_to_remove, collapse = ", ")
  )
}

message(
  "Muestras históricas eliminadas tras construir qcPC1_ctrl: ",
  paste(samples_to_remove, collapse = ", ")
)
samples_used <- setdiff(rownames(metadata), samples_to_remove)
metadata_clean <- metadata[samples_used, , drop = FALSE]
counts_clean <- counts[, samples_used, drop = FALSE]

samples_used_table <- data.frame(
  sample = rownames(metadata_clean),
  condition = metadata_clean$condition,
  condition.run = metadata_clean$condition.run,
  qcPC1_ctrl = metadata_clean$qcPC1_ctrl,
  check.names = FALSE
)
write_tsv(
  samples_used_table,
  file.path(bulk_table_dir, "Figure4K_samples_used.tsv")
)

# ============================================================
# Filtros y funciones de DESeq2
# ============================================================

strict_keep_by_group <- function(
  count_matrix,
  group,
  min_count,
  min_reps
) {
  if (ncol(count_matrix) != length(group)) {
    stop("count_matrix y group no tienen el mismo número de muestras.")
  }
  group <- droplevels(factor(group))
  keep_by_group <- vapply(
    levels(group),
    function(current_group) {
      rowSums(
        count_matrix[, group == current_group, drop = FALSE] >= min_count
      ) >= min_reps
    },
    logical(nrow(count_matrix))
  )
  if (is.null(dim(keep_by_group))) {
    return(keep_by_group)
  }
  rowSums(keep_by_group) >= 1L
}

historical_gene_name_keep <- function(genes) {
  !grepl(
    "^(Rpl|Rps|Mrpl|Mrps)|-ps|^Gm[0-9]|Rik$",
    genes
  )
}

filter_counts <- function(
  count_matrix,
  group,
  min_count,
  min_reps,
  contrast_name
) {
  message(contrast_name, ": genes antes del filtrado: ", nrow(count_matrix))
  keep_expression <- strict_keep_by_group(
    count_matrix,
    group,
    min_count,
    min_reps
  )
  keep_name <- historical_gene_name_keep(rownames(count_matrix))
  filtered <- count_matrix[keep_expression & keep_name, , drop = FALSE]
  message(contrast_name, ": genes después del filtrado: ", nrow(filtered))
  if (nrow(filtered) == 0L) {
    stop(contrast_name, ": ningún gen superó los filtros.")
  }
  filtered
}

format_deseq_result <- function(result) {
  data.frame(
    gene = rownames(result),
    as.data.frame(result),
    row.names = NULL,
    check.names = FALSE
  )
}

make_stat_ranking <- function(result, contrast_name) {
  ranking <- result$stat
  names(ranking) <- rownames(result)
  ranking <- ranking[!is.na(ranking) & is.finite(ranking)]
  ranking <- ranking[!duplicated(names(ranking))]
  ranking <- sort(ranking, decreasing = TRUE)
  if (length(ranking) == 0L) {
    stop(contrast_name, ": no hay valores STAT válidos para GSEA.")
  }
  ranking
}

# ============================================================
# DESeq2: CM vs Control
# ============================================================

cm_keep <- metadata_clean$condition.run %in% c("Control 2", "CM 2")
metadata_cm <- metadata_clean[cm_keep, , drop = FALSE]
counts_cm <- counts_clean[, rownames(metadata_cm), drop = FALSE]

unexpected_cm_conditions <- setdiff(
  unique(as.character(metadata_cm$condition)),
  c("Control", "CM")
)
if (length(unexpected_cm_conditions) > 0L) {
  stop(
    "El contraste CM contiene condiciones biológicas inesperadas: ",
    paste(unexpected_cm_conditions, collapse = ", ")
  )
}
if (!setequal(unique(as.character(metadata_cm$condition)), c("Control", "CM"))) {
  stop("El contraste CM debe contener exclusivamente Control y CM.")
}

metadata_cm$condition_simple <- factor(
  metadata_cm$condition,
  levels = c("Control", "CM")
)
message(
  "CM vs Control: número de muestras por condición: ",
  paste(names(table(metadata_cm$condition_simple)), table(metadata_cm$condition_simple), collapse = "; ")
)
if (any(table(metadata_cm$condition_simple) < 3L)) {
  stop("CM vs Control requiere al menos tres réplicas por condición.")
}

counts_cm <- filter_counts(
  counts_cm,
  metadata_cm$condition_simple,
  min_count = 15,
  min_reps = 3,
  contrast_name = "CM vs Control"
)

# Control and CM are restricted to the rerun, so both groups share the same
# technical run and the historical model uses condition only.
dds_cm <- DESeqDataSetFromMatrix(
  countData = counts_cm,
  colData = metadata_cm,
  design = ~ condition_simple
)
dds_cm <- DESeq(dds_cm)
result_cm <- results(
  dds_cm,
  contrast = c("condition_simple", "CM", "Control")
)
write_tsv(
  format_deseq_result(result_cm),
  file.path(bulk_table_dir, "Figure4K_CM_vs_Control_DESeq2.tsv")
)

# ============================================================
# DESeq2: Bencetinib and CM + Bencetinib vs Control
# ============================================================

bct_conditions <- c("Control", "Bencetinib", "CM + Bencetinib")
bct_keep <- metadata_clean$condition %in% bct_conditions
metadata_bct <- metadata_clean[bct_keep, , drop = FALSE]
counts_bct <- counts_clean[, rownames(metadata_bct), drop = FALSE]

if (!setequal(unique(as.character(metadata_bct$condition)), bct_conditions)) {
  stop(
    "El análisis BCT debe contener Control, Bencetinib y ",
    "CM + Bencetinib."
  )
}

condition_simple_bct <- as.character(metadata_bct$condition)
condition_simple_bct[condition_simple_bct == "CM + Bencetinib"] <- "CMplusBenc"
metadata_bct$condition_simple <- factor(
  condition_simple_bct,
  levels = c("Control", "Bencetinib", "CMplusBenc")
)
message(
  "BCT/CM+BCT vs Control: número de muestras por condición: ",
  paste(names(table(metadata_bct$condition_simple)), table(metadata_bct$condition_simple), collapse = "; ")
)
if (any(table(metadata_bct$condition_simple) < 2L)) {
  stop("BCT/CM+BCT vs Control requiere al menos dos réplicas por condición.")
}

counts_bct <- filter_counts(
  counts_bct,
  metadata_bct$condition_simple,
  min_count = 15,
  min_reps = 2,
  contrast_name = "BCT/CM+BCT vs Control"
)

# These contrasts compare treatments from the original run with valid Controls
# from the rerun. `qcPC1_ctrl` is therefore included as the historical
# continuous technical covariate used to model the measured technical
# variation between runs.
bct_model_matrix <- model.matrix(
  ~ qcPC1_ctrl + condition_simple,
  data = metadata_bct
)
if (qr(bct_model_matrix)$rank != ncol(bct_model_matrix)) {
  stop(
    "El diseño ~ qcPC1_ctrl + condition_simple no tiene rango completo. ",
    "No se ha modificado automáticamente la fórmula."
  )
}

dds_bct <- DESeqDataSetFromMatrix(
  countData = counts_bct,
  colData = metadata_bct,
  design = ~ qcPC1_ctrl + condition_simple
)
dds_bct <- DESeq(dds_bct)
result_bct <- results(
  dds_bct,
  contrast = c("condition_simple", "Bencetinib", "Control")
)
result_cmplbct <- results(
  dds_bct,
  contrast = c("condition_simple", "CMplusBenc", "Control")
)
write_tsv(
  format_deseq_result(result_bct),
  file.path(bulk_table_dir, "Figure4K_BCT_vs_Control_DESeq2.tsv")
)
write_tsv(
  format_deseq_result(result_cmplbct),
  file.path(bulk_table_dir, "Figure4K_CMplusBCT_vs_Control_DESeq2.tsv")
)

# ============================================================
# Rankings DESeq2 STAT
# ============================================================

rankings <- list(
  "CM vs Ctrl" = make_stat_ranking(result_cm, "CM vs Control"),
  "BCT vs Ctrl" = make_stat_ranking(result_bct, "BCT vs Control"),
  "CM+BCT vs Ctrl" = make_stat_ranking(
    result_cmplbct,
    "CM+BCT vs Control"
  )
)

# ============================================================
# Hallmark GSEA
# ============================================================

hallmark_sets_table <- msigdbr(
  species = "Mus musculus",
  collection = "H"
)
required_hallmark_columns <- c("gs_name", "gene_symbol")
missing_hallmark_columns <- setdiff(
  required_hallmark_columns,
  colnames(hallmark_sets_table)
)
if (length(missing_hallmark_columns) > 0L) {
  stop(
    "Faltan columnas necesarias en msigdbr Hallmark: ",
    paste(missing_hallmark_columns, collapse = ", ")
  )
}
hallmark_sets <- split(
  hallmark_sets_table$gene_symbol,
  hallmark_sets_table$gs_name
)

run_hallmark_gsea <- function(ranking, contrast_name) {
  set.seed(123)
  output <- fgsea(
    pathways = hallmark_sets,
    stats = ranking,
    minSize = 10,
    maxSize = 500,
    eps = 0
  )
  output <- as.data.frame(output)
  required_columns <- c("pathway", "NES", "padj")
  missing_columns <- setdiff(required_columns, colnames(output))
  if (length(missing_columns) > 0L) {
    stop(
      contrast_name,
      ": faltan columnas necesarias en Hallmark GSEA: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  if ("leadingEdge" %in% colnames(output)) {
    output$leadingEdge <- vapply(
      output$leadingEdge,
      paste,
      collapse = ",",
      FUN.VALUE = character(1)
    )
  }
  output
}

hallmark_results <- Map(run_hallmark_gsea, rankings, names(rankings))
names(hallmark_results) <- names(rankings)
hallmark_filenames <- c(
  "CM vs Ctrl" = "Figure4K_CM_vs_Control_Hallmark_GSEA.tsv",
  "BCT vs Ctrl" = "Figure4K_BCT_vs_Control_Hallmark_GSEA.tsv",
  "CM+BCT vs Ctrl" = "Figure4K_CMplusBCT_vs_Control_Hallmark_GSEA.tsv"
)
for (contrast_name in names(hallmark_results)) {
  write_tsv(
    hallmark_results[[contrast_name]],
    file.path(bulk_table_dir, hallmark_filenames[[contrast_name]])
  )
}

# ============================================================
# GO Biological Process GSEA
# ============================================================

make_go_ranking <- function(ranking, contrast_name) {
  mapping <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = names(ranking),
    keytype = "SYMBOL",
    columns = "ENTREZID"
  )
  required_columns <- c("SYMBOL", "ENTREZID")
  missing_columns <- setdiff(required_columns, colnames(mapping))
  if (length(missing_columns) > 0L) {
    stop(
      contrast_name,
      ": faltan columnas en la conversión SYMBOL -> ENTREZID: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  mapping <- mapping[!is.na(mapping$ENTREZID), , drop = FALSE]
  mapping$STAT <- unname(ranking[mapping$SYMBOL])
  mapping <- mapping[is.finite(mapping$STAT), , drop = FALSE]
  mapping <- mapping[
    order(abs(mapping$STAT), decreasing = TRUE),
    ,
    drop = FALSE
  ]
  mapping <- mapping[!duplicated(mapping$ENTREZID), , drop = FALSE]
  gene_list <- mapping$STAT
  names(gene_list) <- mapping$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  if (length(gene_list) == 0L) {
    stop(contrast_name, ": no quedan ENTREZID válidos para GO GSEA.")
  }
  gene_list
}

run_go_gsea <- function(ranking, contrast_name) {
  geneList_go <- make_go_ranking(ranking, contrast_name)
  output <- gseGO(
    geneList = geneList_go,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
  output <- as.data.frame(output)
  required_columns <- c("Description", "NES", "p.adjust")
  missing_columns <- setdiff(required_columns, colnames(output))
  if (length(missing_columns) > 0L) {
    stop(
      contrast_name,
      ": faltan columnas necesarias en GO BP GSEA: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  output
}

go_results <- Map(run_go_gsea, rankings, names(rankings))
names(go_results) <- names(rankings)
go_filenames <- c(
  "CM vs Ctrl" = "Figure4K_CM_vs_Control_GO_BP_GSEA.tsv",
  "BCT vs Ctrl" = "Figure4K_BCT_vs_Control_GO_BP_GSEA.tsv",
  "CM+BCT vs Ctrl" = "Figure4K_CMplusBCT_vs_Control_GO_BP_GSEA.tsv"
)
for (contrast_name in names(go_results)) {
  write_tsv(
    go_results[[contrast_name]],
    file.path(bulk_table_dir, go_filenames[[contrast_name]])
  )
}

# ============================================================
# Pathways fijos y source data de Figure 4K
# ============================================================

go_selected <- c(
  "cell junction assembly",
  "actin filament-based process",
  "cell projection morphogenesis",
  "actin filament organization",
  "oxidative phosphorylation",
  "actomyosin structure organization",
  "toll-like receptor signaling pathway",
  "response to growth factor",
  "regulation of cell-matrix adhesion",
  "cell-matrix adhesion",
  "blood vessel development"
)

hallmark_selected <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_APICAL_JUNCTION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_HEDGEHOG_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"
)

extract_selected <- function(
  result,
  selected,
  pathway_column,
  padj_column,
  source,
  contrast_name
) {
  required_columns <- c(pathway_column, "NES", padj_column)
  missing_columns <- setdiff(required_columns, colnames(result))
  if (length(missing_columns) > 0L) {
    stop(
      contrast_name,
      ": faltan columnas para seleccionar pathways de ",
      source,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
  match_index <- match(selected, result[[pathway_column]])
  missing_pathways <- selected[is.na(match_index)]
  if (length(missing_pathways) > 0L) {
    warning(
      contrast_name,
      " - ",
      source,
      ": pathways seleccionados ausentes (se conservarán como NA): ",
      paste(missing_pathways, collapse = ", ")
    )
  }
  pathway_label <- selected
  if (source == "Hallmark") {
    pathway_label <- gsub("_", " ", sub("^HALLMARK_", "", selected))
  }
  data.frame(
    source = source,
    pathway_raw = selected,
    pathway = pathway_label,
    contrast = contrast_name,
    NES = result$NES[match_index],
    padj = result[[padj_column]][match_index],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

source_data_parts <- list()
for (contrast_name in names(rankings)) {
  source_data_parts[[paste("GO", contrast_name)]] <- extract_selected(
    go_results[[contrast_name]],
    go_selected,
    pathway_column = "Description",
    padj_column = "p.adjust",
    source = "GO",
    contrast_name = contrast_name
  )
  source_data_parts[[paste("Hallmark", contrast_name)]] <- extract_selected(
    hallmark_results[[contrast_name]],
    hallmark_selected,
    pathway_column = "pathway",
    padj_column = "padj",
    source = "Hallmark",
    contrast_name = contrast_name
  )
}
source_data <- do.call(rbind, source_data_parts)
rownames(source_data) <- NULL
write_tsv(
  source_data,
  file.path(table_dir, "Figure4K_source_data.tsv")
)

# ============================================================
# Heatmap Figure 4K
# ============================================================

contrast_order <- c("CM vs Ctrl", "BCT vs Ctrl", "CM+BCT vs Ctrl")
pathway_key_table <- unique(
  source_data[, c("source", "pathway_raw", "pathway"), drop = FALSE]
)
pathway_key_table$key <- paste(
  pathway_key_table$source,
  pathway_key_table$pathway_raw,
  sep = "::"
)

nes_mat <- matrix(
  NA_real_,
  nrow = nrow(pathway_key_table),
  ncol = length(contrast_order),
  dimnames = list(pathway_key_table$key, contrast_order)
)
for (current_row in seq_len(nrow(source_data))) {
  current_key <- paste(
    source_data$source[current_row],
    source_data$pathway_raw[current_row],
    sep = "::"
  )
  nes_mat[current_key, source_data$contrast[current_row]] <-
    source_data$NES[current_row]
}

if (all(is.na(nes_mat))) {
  stop("Todos los NES seleccionados son NA; no se puede crear Figure 4K.")
}

row_distance <- as.dist(
  1 - cor(
    t(nes_mat),
    method = "pearson"
  )
)
if (any(!is.finite(as.vector(row_distance)))) {
  stop(
    "La distancia de clustering contiene valores no finitos. ",
    "Revise los pathways ausentes/NA sin sustituirlos automáticamente."
  )
}

lim <- max(abs(nes_mat), na.rm = TRUE)
if (!is.finite(lim) || lim <= 0) {
  stop("No se puede construir una escala simétrica con los NES disponibles.")
}

annotation_row <- data.frame(
  Source = pathway_key_table$source,
  row.names = pathway_key_table$key,
  check.names = FALSE
)
annotation_colors <- list(
  Source = c(
    GO = "#B4EEB4",
    Hallmark = "#FFA07A"
  )
)

# The analytical NES matrix is retained with pathways in rows and contrasts in
# columns. It is transposed only for the horizontal presentation of Figure 4K.
nes_mat_horizontal <- t(nes_mat)

Figure4K <- pheatmap(
  nes_mat_horizontal,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  clustering_distance_cols = row_distance,
  clustering_method = "average",
  color = colorRampPalette(
    c("#0571B0", "#FFFFFF", "#CA0020")
  )(101),
  breaks = seq(-lim, lim, length.out = 102),
  na_col = "grey90",
  annotation_col = annotation_row,
  annotation_colors = annotation_colors,
  labels_col = pathway_key_table$pathway,
  display_numbers = FALSE,
  border_color = NA,
  angle_col = 45,
  fontsize_col = 9,
  silent = TRUE
)



figure_file <- file.path(
  figure_dir,
  "Figure4K_selected_GO_Hallmark_3contrasts.pdf"
)

grDevices::pdf(
  figure_file,
  width = 18,
  height = 5
)

grid::grid.newpage()

# Leave additional space on the left for the first rotated pathway label.
grid::pushViewport(
  grid::viewport(
    x = 0.55,
    y = 0.54,
    width = 0.88,
    height = 0.88
  )
)

grid::grid.draw(Figure4K$gtable)

grid::popViewport()

grDevices::dev.off()

message("Figure 4K completada.")
