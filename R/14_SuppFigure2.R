# ============================================================
# Supplementary Figure 2
# CD45-negative compartment and external fibroblast validation
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(grid)
})

set.seed(123)

# ============================================================
# Configuration and output directories
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
config_values <- setNames(config$value, config$key)
rsdir <- unname(config_values[["rsdir"]])
outdir <- unname(config_values[["outdir"]])

plot_dir <- file.path(outdir, "SuppFigure2")
table_dir <- file.path(rsdir, "tables", "SuppFigure2")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(data, filename) {
  write.table(
    data,
    file = filename,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

# ============================================================
# Input for Supplementary Figure 2A
# ============================================================

main_obj <- readRDS(
  file.path(
    "results",
    "objects",
    "fibroblasts.Not.Integrated.annotated.rds"
  )
)

required_meta <- c("Clustering.Mesenchymal", "orig.ident")
missing_meta <- setdiff(required_meta, colnames(main_obj[[]]))
if (length(missing_meta) > 0L) {
  stop("Faltan columnas en el objeto anotado final: ", paste(missing_meta, collapse = ", "))
}
if (!("umap" %in% Reductions(main_obj))) {
  stop("El objeto anotado final no contiene la reducción UMAP global 'umap'.")
}

# ============================================================
# Global palette used by the historical workflow
# ============================================================

S2A_colors <- c(
  "Endothelial Vascular" = "#4A90E2",
  "Aerocytes" = "#1B4F72",
  "Epithelial AT1" = "#2D6A4F",
  "Epithelial AT2" = "#95D5B2",
  "Ciliated" = "#9B5DE5",
  "Prolif. Cancer Cells" = "#a8bfce",
  "EMT Cancer Cells" = "#bb6666",
  "Mesenchymal" = "#FFD6A5"
)

# ============================================================
# Panel A - CD45-negative UMAP split by genotype
# ============================================================

S2A_annotation_order <- c(
  "Endothelial Vascular",
  "Aerocytes",
  "Epithelial AT1",
  "Epithelial AT2",
  "Ciliated",
  "Mesenchymal",
  "Prolif. Cancer Cells",
  "EMT Cancer Cells"
)

main_meta <- main_obj[[]]
S2A_cells <- rownames(main_meta)[
  !is.na(main_meta$Clustering.Mesenchymal) &
    as.character(main_meta$Clustering.Mesenchymal) %in% S2A_annotation_order
]

cd45neg <- subset(main_obj, cells = S2A_cells)
cd45neg_meta <- cd45neg[[]] %>%
  mutate(
    cell = rownames(.),
    genotype = factor(as.character(orig.ident), levels = c("WT", "CD169DTR", "CCR2KO")),
    annotation = as.character(Clustering.Mesenchymal)
  )

annotation_order <- S2A_annotation_order

umap_df <- as.data.frame(Embeddings(cd45neg, reduction = "umap"))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$cell <- rownames(umap_df)
umap_df <- left_join(
  umap_df,
  cd45neg_meta %>% select(cell, genotype, annotation),
  by = "cell"
)
umap_df$annotation <- factor(umap_df$annotation, levels = annotation_order)

figure2a <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = annotation)) +
  geom_point(size = 0.18, alpha = 0.75) +
  facet_wrap(~genotype, nrow = 1) +
  scale_color_manual(values = S2A_colors, drop = FALSE) +
  coord_equal() +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 3)
    )
  ) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Annotation")

ggsave(
  file.path(plot_dir, "SuppFigure2A_CD45neg_UMAP_by_Genotype.pdf"),
  figure2a, width = 17, height = 5.5, units = "in", device = grDevices::cairo_pdf
)
write_tsv(umap_df, file.path(table_dir, "SuppFigure2A_source_data.tsv"))

# ============================================================
# Supplementary Figure 2B - CD45-negative marker bubbleplot
# ============================================================

cd45neg_clusters <- c(
  "Endothelial Vascular",
  "Aerocytes",
  "Epithelial AT1",
  "Epithelial AT2",
  "Ciliated",
  "Prolif. Cancer Cells",
  "EMT Cancer Cells",
  "Mesenchymal"
)

main_meta_s2b <- main_obj[[]]
S2B_cells <- rownames(main_meta_s2b)[
  !is.na(main_meta_s2b$Clustering.Mesenchymal) &
    as.character(main_meta_s2b$Clustering.Mesenchymal) %in% cd45neg_clusters
]

cd45neg_s2b <- subset(main_obj, cells = S2B_cells)
cd45neg_s2b$Cell_Type <- factor(
  as.character(cd45neg_s2b$Clustering.Mesenchymal),
  levels = cd45neg_clusters
)

genes_interes <- c(
  "Ptprb", "Vwf",
  "Car4", "Ednrb",
  "Ager", "Pdpn",
  "Sftpb", "Sftpc",
  "Foxj1", "Dynlrb2",
  "Clu", "Krt18", "Cdkn2a",
  "Vim", "Lgals1", "Hmga2", "Cpe",
  "Upk3b", "Col1a1", "Dcn"
)

genes_present <- genes_interes[
  genes_interes %in% rownames(cd45neg_s2b[["RNA"]])
]

if (length(genes_present) == 0L) {
  stop("No hay genes históricos de S2B presentes en el assay RNA.")
}

dot_plot_s2b <- DotPlot(
  cd45neg_s2b,
  features = genes_present,
  group.by = "Cell_Type",
  assay = "RNA",
  scale = FALSE
)

dot_table_s2b <- dot_plot_s2b$data %>%
  rename(
    Gene = features.plot,
    Cell_Type = id,
    AvgExpression = avg.exp,
    PctExpressed = pct.exp
  )

dot_table_s2b$Gene <- factor(
  dot_table_s2b$Gene,
  levels = genes_present
)

dot_table_s2b$Cell_Type <- factor(
  dot_table_s2b$Cell_Type,
  levels = rev(cd45neg_clusters)
)

color_cap_s2b <- quantile(
  dot_table_s2b$AvgExpression,
  probs = 0.85,
  na.rm = TRUE
)

dot_table_s2b$AvgExpression_plot <- pmin(
  dot_table_s2b$AvgExpression,
  color_cap_s2b
)

write_tsv(
  dot_table_s2b,
  file.path(table_dir, "SuppFigure2B_source_data.tsv")
)

bubble_cols_s2b <- c(
  "#4575B4",
  "#91BFDB",
  "#FFFFBF",
  "#FDAE61",
  "#D73027"
)

bubble_values_s2b <- c(0.00, 0.18, 0.40, 0.65, 1.00)

figure2b <- ggplot(
  dot_table_s2b,
  aes(x = Gene, y = Cell_Type)
) +
  geom_point(
    aes(size = PctExpressed, color = AvgExpression_plot),
    alpha = 1
  ) +
  scale_size(
    range = c(0.2, 8),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    name = "Percentage Expressed"
  ) +
  scale_color_gradientn(
    colors = bubble_cols_s2b,
    values = bubble_values_s2b,
    limits = c(0, color_cap_s2b),
    breaks = round(pretty(c(0, color_cap_s2b), n = 5), 1),
    oob = scales::squish,
    name = "Average Expression"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.45),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9
    ),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  labs(x = "Gene", y = "Cluster")

ggsave(
  file.path(plot_dir, "SuppFigure2B_CD45neg_Cell_Type_Bubbleplot.pdf"),
  figure2b,
  width = max(12, length(genes_present) * 0.38),
  height = 6,
  units = "in",
  device = grDevices::cairo_pdf
)

# ============================================================
# Supplementary Figure 2C - CD45-negative population composition
# ============================================================

S2C_genotype_order <- c(
  "WT",
  "CD169DTR",
  "CCR2KO"
)

composition_s2c <- cd45neg[[]] %>%
  transmute(
    genotype = factor(
      as.character(orig.ident),
      levels = S2C_genotype_order
    ),
    condition = factor(
      as.character(group),
      levels = c("Naive", "Tumor")
    ),
    annotation = factor(
      as.character(Clustering.Mesenchymal),
      levels = cd45neg_clusters
    )
  ) %>%
  filter(!is.na(genotype), !is.na(condition), !is.na(annotation)) %>%
  count(genotype, condition, annotation, name = "n")

composition_s2c_complete <- expand.grid(
  genotype = factor(
    S2C_genotype_order,
    levels = S2C_genotype_order
  ),
  condition = factor(
    c("Naive", "Tumor"),
    levels = c("Naive", "Tumor")
  ),
  annotation = factor(
    cd45neg_clusters,
    levels = cd45neg_clusters
  )
) %>%
  left_join(
    composition_s2c,
    by = c("genotype", "condition", "annotation")
  ) %>%
  mutate(n = ifelse(is.na(n), 0L, n)) %>%
  group_by(genotype, condition) %>%
  mutate(
    proportion = n / sum(n)
  ) %>%
  ungroup()

figure2c <- ggplot(
  composition_s2c_complete,
  aes(x = condition, y = proportion, fill = annotation)
) +
  geom_col(width = 0.9) +
  facet_grid(. ~ genotype) +
  scale_fill_manual(
    values = S2A_colors,
    breaks = cd45neg_clusters,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = function(x) paste0(round(x * 100), "%")
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  labs(
    x = "Condition",
    y = "Proportion",
    fill = "Annotation"
  )

ggsave(
  file.path(plot_dir, "SuppFigure2C_CD45neg_Composition_by_Genotype.pdf"),
  figure2c,
  width = 8,
  height = 5.5,
  units = "in",
  device = grDevices::cairo_pdf
)

write_tsv(
  composition_s2c_complete,
  file.path(table_dir, "SuppFigure2C_source_data.tsv")
)

# ============================================================
# External validation helpers for Supplementary Figures 2D-F
# ============================================================

cardoso_file <- file.path(
  "data", "SuppFigure2", "Fibros.External",
  "Cardoso", "Total_Integration_fibroblasts.rds"
)

steady_file <- file.path(
  "data", "SuppFigure2", "Fibros.External",
  "Shannon", "Mouse_SS_Fibro.RDS"
)

perturbed_file <- file.path(
  "data", "SuppFigure2", "Fibros.External",
  "Shannon", "Mouse_PS_Fibro.RDS"
)

human_perturbed_file <- file.path(
  "data", "SuppFigure2", "Fibros.External",
  "Shannon", "Human_PS_Fibro.RDS"
)

tsukui_file <- file.path(
  "data", "SuppFigure2", "Fibros.External",
  "Sheppard", "GSM6428697_seurat_v4_obj.rds"
)

fib_univ_mouse <- c(
  "Col1a2", "Clec3b", "Serping1", "Dpep1", "Slit3",
  "Fn1", "Scara5", "Cd248", "Pi16", "Abca8a",
  "C7", "Pcolce", "Serpinf1", "Svep1", "Col15a1",
  "Ebf2", "Axl", "Pdpn"
)

fib_act1_mouse <- c(
  "Hhip", "Grem2", "Bmp5", "Mustn1", "Angpt2",
  "Wif1", "Gper1", "Kcnma1", "Lgr6", "P3h2",
  "Slc8a3", "Tagln", "Acta2", "Tnc", "Myh11"
)

fib_act2_mouse <- c(
  "Ddit4", "Nampt", "Per1", "Klf15", "Npr3",
  "Cdo1", "Fmo2", "Map1lc3b", "Glul", "Trp53inp1",
  "Per2", "Mt1"
)

fib_univ_human <- c(
  "COL1A2", "CLEC3B", "SERPING1", "DPEP1", "SLIT3",
  "FN1", "SCARA5", "CD248", "PI16", "ABCA8",
  "C7", "PCOLCE", "SERPINF1", "SVEP1", "COL15A1",
  "EBF2", "AXL", "PDPN"
)

fib_act1_human <- c(
  "HHIP", "GREM2", "BMP5", "MUSTN1", "ANGPT2",
  "WIF1", "GPER1", "KCNMA1", "LGR6", "P3H2",
  "SLC8A3", "TAGLN", "ACTA2", "TNC", "MYH11"
)

fib_act2_human <- c(
  "DDIT4", "NAMPT", "PER1", "KLF15", "NPR3",
  "CDO1", "FMO2", "MAP1LC3B", "GLUL", "TP53INP1",
  "PER2", "MT1X"
)

fib_signatures_mouse <- list(
  Fib_Univ = fib_univ_mouse,
  Fib_Act1 = fib_act1_mouse,
  Fib_Act2 = fib_act2_mouse
)

fib_signatures_human <- list(
  Fib_Univ = fib_univ_human,
  Fib_Act1 = fib_act1_human,
  Fib_Act2 = fib_act2_human
)

match_signature_to_object <- function(signature_genes, object_genes) {
  signature_genes <- unique(signature_genes)
  exact_found <- intersect(signature_genes, object_genes)
  remaining <- setdiff(signature_genes, exact_found)

  if (length(remaining) > 0L) {
    idx <- match(toupper(remaining), toupper(object_genes))
    ci_found <- object_genes[idx[!is.na(idx)]]
  } else {
    ci_found <- character()
  }

  unique(c(exact_found, ci_found))
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)

  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }

  (x - mean(x, na.rm = TRUE)) / s
}

add_signature_scores <- function(object, signatures, assay_use = "RNA") {
  DefaultAssay(object) <- assay_use
  set.seed(123)

  for (signature_name in names(signatures)) {
    genes_found <- match_signature_to_object(
      signatures[[signature_name]],
      rownames(object)
    )

    if (length(genes_found) == 0L) {
      stop("No se encontraron genes para la firma ", signature_name)
    }

    object <- AddModuleScore(
      object,
      features = list(genes_found),
      assay = assay_use,
      name = paste0(signature_name, "_AMS"),
      ctrl = 50
    )

    object[[paste0(signature_name, "_score")]] <- object[[
      paste0(signature_name, "_AMS1")
    ]]
    object[[paste0(signature_name, "_z")]] <- zscore(
      object[[paste0(signature_name, "_score")]][, 1]
    )
  }

  object
}

get_umap_reduction <- function(obj) {
  if ("umap" %in% Reductions(obj)) return("umap")
  if ("umap.harmony" %in% Reductions(obj)) return("umap.harmony")
  Reductions(obj)[1]
}

geom_umap_point <- function(...) {
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::geom_point_rast(..., raster.dpi = 300)
  } else {
    geom_point(...)
  }
}

gradient_cols <- c(
  "#313695",
  "#74ADD1",
  "#F7F7F7",
  "#FDAE61",
  "#A50026"
)

make_signature_umap <- function(
  obj,
  title,
  first_signature_includes_dataset = FALSE,
  show_signature_titles = TRUE
) {
  reduction_use <- get_umap_reduction(obj)
  umap_df <- as.data.frame(Embeddings(obj, reduction = reduction_use))
  colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
  umap_df$cell <- rownames(umap_df)

  signature_labels <- c(
    Fib_Univ_z = "Fib. Univ",
    Fib_Act1_z = "Fib. Act1",
    Fib_Act2_z = "Fib. Act2"
  )

  plot_list <- lapply(names(signature_labels), function(signature_name) {
    plot_df <- umap_df
    plot_df$z_score <- obj[[]][plot_df$cell, signature_name]

    signature_title <- if (!show_signature_titles) {
      NULL
    } else if (
      first_signature_includes_dataset && signature_name == "Fib_Univ_z"
    ) {
      paste(title, signature_labels[[signature_name]], sep = "\n")
    } else {
      signature_labels[[signature_name]]
    }

    ggplot(plot_df, aes(UMAP_1, UMAP_2, color = z_score)) +
      geom_umap_point(size = 0.22, alpha = 0.9) +
      scale_color_gradientn(
        colors = gradient_cols,
        values = scales::rescale(c(-2, -1, 0, 1, 2)),
        limits = c(-2, 2),
        oob = scales::squish,
        name = "z-score"
      ) +
      coord_equal() +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        plot.margin = margin(2, 4, 2, 2)
      ) +
      labs(
        title = signature_title,
        x = "UMAP 1",
        y = "UMAP 2"
      )
  })

  names(plot_list) <- names(signature_labels)
  plot_list
}

make_signature_source <- function(obj, dataset, group_column = NULL) {
  reduction_use <- get_umap_reduction(obj)
  umap_df <- as.data.frame(Embeddings(obj, reduction = reduction_use))
  colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
  umap_df$cell <- rownames(umap_df)
  meta <- obj[[]][umap_df$cell, , drop = FALSE]

  source <- data.frame(
    dataset = dataset,
    cell = umap_df$cell,
    UMAP_1 = umap_df$UMAP_1,
    UMAP_2 = umap_df$UMAP_2,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  if (!is.null(group_column)) {
    source[[group_column]] <- as.character(meta[[group_column]])
  }

  source$Fib_Univ_z <- as.numeric(meta$Fib_Univ_z)
  source$Fib_Act1_z <- as.numeric(meta$Fib_Act1_z)
  source$Fib_Act2_z <- as.numeric(meta$Fib_Act2_z)
  source
}

shared_cluster_cols <- c(
  "Pi16" = "#CC79A7",
  "Ccl19" = "#D55E00",
  "Col15a1" = "#009E73",
  "Comp" = "#00A087",
  "Cxcl12" = "#0072B2",
  "Hhip" = "#3B4CC0",
  "Npnt" = "#984EA3",
  "Bmp4" = "#E64B35",
  "Coch" = "#A6A600",
  "Fbln1" = "#56B4E9",
  "Adamdec1" = "#E64B35",
  "Cxcl5" = "#F39B7F",
  "Lrrc15" = "#DC0000",
  "PI16" = "#CC79A7",
  "CCL19" = "#D55E00",
  "COL15A1" = "#009E73",
  "COL3A1" = "#66A61E",
  "LRRC15" = "#DC0000",
  "NPNT" = "#984EA3",
  "ADAMDEC1" = "#E64B35"
)

make_clustname_umap <- function(obj, title) {
  reduction_use <- get_umap_reduction(obj)
  umap_df <- as.data.frame(Embeddings(obj, reduction = reduction_use))
  colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
  umap_df$ClustName <- as.character(
    obj[[]][rownames(umap_df), "ClustName"]
  )

  cols_use <- shared_cluster_cols
  missing_cols <- setdiff(unique(umap_df$ClustName), names(cols_use))
  if (length(missing_cols) > 0L) {
    cols_use[missing_cols] <- "grey70"
  }

  label_df <- umap_df %>%
    group_by(ClustName) %>%
    summarise(
      UMAP_1 = median(UMAP_1, na.rm = TRUE),
      UMAP_2 = median(UMAP_2, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot(umap_df, aes(UMAP_1, UMAP_2, color = ClustName)) +
    geom_umap_point(size = 0.20, alpha = 0.9) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(UMAP_1, UMAP_2, label = ClustName),
      inherit.aes = FALSE,
      size = 2.2,
      color = "black",
      fontface = "bold",
      box.padding = 0.22,
      point.padding = 0.12,
      max.overlaps = Inf,
      segment.color = NA
    ) +
    scale_color_manual(values = cols_use, drop = FALSE) +
    coord_equal() +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 6.5),
      axis.text = element_text(size = 5.8),
      plot.title = element_text(size = 8.2, hjust = 0.5),
      plot.margin = margin(1, 1, 1, 1)
    ) +
    labs(title = title, color = NULL)
}

# ============================================================
# Supplementary Figure 2D - Buechler et al.
# ============================================================

steady <- readRDS(steady_file)
steady <- add_signature_scores(
  steady,
  fib_signatures_mouse,
  assay_use = "RNA"
)
steady_source <- make_signature_source(steady, "Mouse steady", "ClustName")
steady_original <- make_clustname_umap(steady, "Mouse steady")
steady_signatures <- make_signature_umap(
  steady,
  "Mouse steady"
)
rm(steady)
gc()

perturbed <- readRDS(perturbed_file)
perturbed <- add_signature_scores(
  perturbed,
  fib_signatures_mouse,
  assay_use = "RNA"
)
perturbed_source <- make_signature_source(perturbed, "Mouse perturbed", "ClustName")
perturbed_original <- make_clustname_umap(perturbed, "Mouse perturbed")
perturbed_signatures <- make_signature_umap(
  perturbed,
  "",
  show_signature_titles = FALSE
)
rm(perturbed)
gc()

human_perturbed <- readRDS(human_perturbed_file)
human_perturbed <- add_signature_scores(
  human_perturbed,
  fib_signatures_human,
  assay_use = "RNA"
)
human_perturbed_source <- make_signature_source(
  human_perturbed,
  "Human perturbed",
  "ClustName"
)
human_perturbed_original <- make_clustname_umap(
  human_perturbed,
  "Human perturbed"
)
human_perturbed_signatures <- make_signature_umap(
  human_perturbed,
  "",
  show_signature_titles = FALSE
)
rm(human_perturbed)
gc()

figure2d <- wrap_plots(
  steady_original,
  perturbed_original,
  human_perturbed_original,
  steady_signatures$Fib_Univ_z,
  steady_signatures$Fib_Act1_z,
  steady_signatures$Fib_Act2_z,
  perturbed_signatures$Fib_Univ_z,
  perturbed_signatures$Fib_Act1_z,
  perturbed_signatures$Fib_Act2_z,
  human_perturbed_signatures$Fib_Univ_z,
  human_perturbed_signatures$Fib_Act1_z,
  human_perturbed_signatures$Fib_Act2_z,
  ncol = 3,
  byrow = TRUE,
  guides = "collect"
) +
  plot_annotation(
    title = "Buechler et al. external validation"
  ) &
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.tag = element_text(face = "bold")
  )

ggsave(
  file.path(plot_dir, "SuppFigure2D_Buechler_Validation.pdf"),
  figure2d,
  width = 13,
  height = 12,
  units = "in",
  device = grDevices::cairo_pdf
)

write_tsv(
  bind_rows(steady_source, perturbed_source, human_perturbed_source),
  file.path(table_dir, "SuppFigure2D_source_data.tsv")
)

# ============================================================
# Supplementary Figure 2E - Cardoso et al.
# ============================================================

cardoso <- readRDS(cardoso_file)
cardoso <- add_signature_scores(
  cardoso,
  fib_signatures_mouse,
  assay_use = "RNA"
)

adventitial_clusters <- c("11", "8", "2", "9", "3")
reprogrammed_clusters <- c("6")
alveolar_clusters <- c("4", "5", "0", "7", "1", "10")

cardoso$cardoso_state <- dplyr::case_when(
  as.character(cardoso$seurat_clusters) %in% adventitial_clusters ~
    "Adventitial fibroblasts",
  as.character(cardoso$seurat_clusters) %in% reprogrammed_clusters ~
    "Reprogrammed fibroblasts",
  as.character(cardoso$seurat_clusters) %in% alveolar_clusters ~
    "Alveolar fibroblasts",
  TRUE ~ "Other"
)

fibro_cols <- c(
  "Alveolar fibroblasts" = "#dd4e4e",
  "Adventitial fibroblasts" = "#c9b821",
  "Reprogrammed fibroblasts" = "#3d8eda",
  "Other" = "grey70"
)

cardoso_reduction <- get_umap_reduction(cardoso)
cardoso_umap <- as.data.frame(Embeddings(cardoso, reduction = cardoso_reduction))
colnames(cardoso_umap)[1:2] <- c("UMAP_1", "UMAP_2")
cardoso_umap$cell <- rownames(cardoso_umap)
cardoso_umap$cell.identity <- as.character(cardoso$cell.identity[cardoso_umap$cell])
cardoso_umap$seurat_clusters <- as.character(cardoso$seurat_clusters[cardoso_umap$cell])
cardoso_umap$fibro_annotation_paper <- as.character(
  cardoso$cardoso_state[cardoso_umap$cell]
)

cardoso_label_df <- cardoso_umap %>%
  filter(fibro_annotation_paper != "Other") %>%
  group_by(fibro_annotation_paper) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  )

cardoso_original <- ggplot(
  cardoso_umap,
  aes(UMAP_1, UMAP_2, color = fibro_annotation_paper)
) +
  geom_umap_point(size = 0.25, alpha = 0.9) +
  scale_color_manual(values = fibro_cols, drop = FALSE) +
  coord_equal() +
  theme_classic() +
  theme(legend.position = "right") +
  guides(
    color = guide_legend(
      override.aes = list(size = 3.5, alpha = 1)
    )
  ) +
  labs(title = NULL, color = "Fibroblast state")

cardoso_signatures <- make_signature_umap(cardoso, "Cardoso")
cardoso_source <- make_signature_source(cardoso, "Cardoso", "cell.identity") %>%
  mutate(
    seurat_clusters = as.character(cardoso$seurat_clusters[cell]),
    fibro_annotation_paper = as.character(cardoso$cardoso_state[cell])
  )

cardoso_signature_block <- wrap_plots(
  cardoso_signatures$Fib_Univ_z,
  cardoso_signatures$Fib_Act1_z,
  cardoso_signatures$Fib_Act2_z,
  ncol = 3,
  guides = "collect"
)

figure2e <- cardoso_original + cardoso_signature_block +
  plot_layout(widths = c(1, 3)) +
  plot_annotation(
    title = "Cardoso et al. external validation"
  )

ggsave(
  file.path(plot_dir, "SuppFigure2E_Cardoso_Validation.pdf"),
  figure2e,
  width = 12,
  height = 5,
  units = "in",
  device = grDevices::cairo_pdf
)

write_tsv(
  cardoso_source,
  file.path(table_dir, "SuppFigure2E_source_data.tsv")
)

rm(cardoso)
gc()

# ============================================================
# Supplementary Figure 2F - Tsukui et al.
# ============================================================

tsukui <- readRDS(tsukui_file)
tsukui <- add_signature_scores(
  tsukui,
  fib_signatures_mouse,
  assay_use = "RNA"
)

cluster_to_tsukui_state <- c(
  "0" = "Alveolar",
  "1" = "Fibrotic",
  "2" = "Inflammatory",
  "3" = "Stress-activated",
  "4" = "Adventitial",
  "5" = "Smooth muscle",
  "6" = "Peribronchial",
  "7" = "Stress-activated",
  "8" = "Proliferating",
  "9" = "Mesothelial",
  "10" = "Pericyte"
)

tsukui$tsukui_state <- unname(
  cluster_to_tsukui_state[
    as.character(tsukui$RNA_snn_res.0.3)
  ]
)

tsukui_state_order <- c(
  "Alveolar",
  "Inflammatory",
  "Stress-activated",
  "Fibrotic",
  "Adventitial",
  "Peribronchial",
  "Smooth muscle",
  "Pericyte",
  "Mesothelial",
  "Proliferating"
)

tsukui_state_cols <- c(
  "Alveolar" = "#E69F00",
  "Inflammatory" = "#009E73",
  "Stress-activated" = "#CC79A7",
  "Fibrotic" = "#56B4E9",
  "Adventitial" = "#F0E442",
  "Peribronchial" = "#D55E00",
  "Smooth muscle" = "#0072B2",
  "Pericyte" = "#882255",
  "Mesothelial" = "grey70",
  "Proliferating" = "grey35"
)

tsukui$tsukui_state <- factor(
  tsukui$tsukui_state,
  levels = tsukui_state_order
)

tsukui_reduction <- get_umap_reduction(tsukui)
tsukui_umap <- as.data.frame(Embeddings(tsukui, reduction = tsukui_reduction))
colnames(tsukui_umap)[1:2] <- c("UMAP_1", "UMAP_2")
tsukui_umap$cell <- rownames(tsukui_umap)
tsukui_umap$cluster_res03 <- as.character(tsukui$RNA_snn_res.0.3[tsukui_umap$cell])
tsukui_umap$tsukui_state <- as.character(tsukui$tsukui_state[tsukui_umap$cell])

tsukui_label_df <- tsukui_umap %>%
  filter(!is.na(tsukui_state)) %>%
  group_by(tsukui_state) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  )

tsukui_original <- ggplot(
  tsukui_umap,
  aes(UMAP_1, UMAP_2, color = tsukui_state)
) +
  geom_umap_point(size = 0.22, alpha = 0.9) +
  scale_color_manual(values = tsukui_state_cols, drop = FALSE) +
  coord_equal() +
  theme_classic() +
  theme(legend.position = "right") +
  guides(
    color = guide_legend(
      override.aes = list(size = 3.5, alpha = 1)
    )
  ) +
  labs(title = NULL, color = "Cell state")

tsukui_signatures <- make_signature_umap(tsukui, "")
tsukui_source <- make_signature_source(tsukui, "Tsukui") %>%
  mutate(
    cluster_res03 = as.character(tsukui$RNA_snn_res.0.3[cell]),
    tsukui_state = as.character(tsukui$tsukui_state[cell])
  )

tsukui_signature_block <- wrap_plots(
  tsukui_signatures$Fib_Univ_z,
  tsukui_signatures$Fib_Act1_z,
  tsukui_signatures$Fib_Act2_z,
  ncol = 3,
  guides = "collect"
)

figure2f <- tsukui_original + tsukui_signature_block +
  plot_layout(widths = c(1, 3)) +
  plot_annotation(
    title = "Tsukui et al. external validation"
  )

ggsave(
  file.path(plot_dir, "SuppFigure2F_Tsukui_Validation.pdf"),
  figure2f,
  width = 12,
  height = 5,
  units = "in",
  device = grDevices::cairo_pdf
)

write_tsv(
  tsukui_source,
  file.path(table_dir, "SuppFigure2F_source_data.tsv")
)

rm(tsukui)
gc()
