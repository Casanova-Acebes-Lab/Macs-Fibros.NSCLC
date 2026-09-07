# ============================================================
# 10_Figure3.R
# Figure 3E-F: T / lymphoid UMAPs
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

set.seed(1234)

# ============================================================
# Configuration
# ============================================================

config <- read.delim(
  "config.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

values <- setNames(
  config$value,
  config$key
)

rsdir <- values[["rsdir"]]
outdir <- values[["outdir"]]

objects_dir <- file.path(
  rsdir,
  "objects"
)

figure_dir <- file.path(
  outdir,
  "Figure3"
)

table_dir <- file.path(
  rsdir,
  "tables",
  "Figure3"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

downsample_n <- 300L

# ============================================================
# Load final object
# ============================================================

data <- readRDS(
  file.path(
    objects_dir,
    "fibroblasts.Not.Integrated.annotated.rds"
  )
)

if (!"Clustering" %in% colnames(data[[]])) {
  stop("The final object does not contain the metadata column 'Clustering'.")
}

if (!"umap" %in% Reductions(data)) {
  stop("The final object does not contain the UMAP reduction.")
}

# ============================================================
# T / lymphoid order and palette
# ============================================================

tcell_order <- c(
  "Cd4 Early Activated",
  "Cd4 Naive",
  "Cd4 IL7R Memory Like",
  "Treg",
  "Inmature T cells",
  "Stem like T cells",
  "Tpex Transitional",
  "Cd8 Memory",
  "Cd8 Exhausted",
  "Tgd",
  "NKT",
  "NK",
  "ILC"
)

tcell_colors <- c(
  "Cd4 Early Activated" = "#4B7591",
  "Cd4 Naive" = "#2ECC71",
  "Cd4 IL7R Memory Like" = "#008F5A",
  "Treg" = "#FFFF00",
  "Inmature T cells" = "#FFC300",
  "Stem like T cells" = "#00CED1",
  "Tpex Transitional" = "#FF5A36",
  "Cd8 Memory" = "#C95F82",
  "Cd8 Exhausted" = "#6A0DAD",
  "Tgd" = "#80ED99",
  "NKT" = "#003566",
  "NK" = "#1F77B4",
  "ILC" = "#C77D00"
)

# ============================================================
# Harmonize labels from Clustering
# ============================================================

labels <- trimws(
  as.character(data$Clustering)
)

label_aliases <- c(
  "CD4 Early Activated" = "Cd4 Early Activated",
  "CD4 Naive" = "Cd4 Naive",
  "CD4 IL7R Memory Like" = "Cd4 IL7R Memory Like",
  "Immature T Cells" = "Inmature T cells",
  "Immature T cells" = "Inmature T cells",
  "Inmature T Cells" = "Inmature T cells",
  "Stem like T Cells" = "Stem like T cells",
  "CD8 Memory" = "Cd8 Memory",
  "CD8 Exhausted" = "Cd8 Exhausted"
)

alias_index <- match(
  labels,
  names(label_aliases)
)

replace_cells <- !is.na(alias_index)

labels[replace_cells] <- unname(
  label_aliases[
    alias_index[replace_cells]
  ]
)

data$Figure3_annotation <- labels

# ============================================================
# Subset T / lymphoid cells
# ============================================================

cells_tcells <- colnames(data)[
  data$Figure3_annotation %in% tcell_order
]

Tcells <- subset(
  data,
  cells = cells_tcells
)

cat(
  "\nT / lymphoid cells before UMAP filtering: ",
  ncol(Tcells),
  "\n",
  sep = ""
)

# ============================================================
# Remove cells with UMAP_2 > 0
# ============================================================

umap_coordinates <- Embeddings(
  Tcells,
  reduction = "umap"
)

cells_keep <- rownames(umap_coordinates)[
  umap_coordinates[, 2] <= 0
]

Tcells <- subset(
  Tcells,
  cells = cells_keep
)

tcell_order_plot <- tcell_order[
  tcell_order %in%
    unique(as.character(Tcells$Figure3_annotation))
]

Tcells$Figure3_annotation <- factor(
  Tcells$Figure3_annotation,
  levels = tcell_order_plot
)

Tcells$orig.ident <- factor(
  as.character(Tcells$orig.ident),
  levels = c(
    "WT",
    "CD169DTR",
    "CCR2KO"
  )
)

tcell_colors_plot <- tcell_colors[
  tcell_order_plot
]

# ============================================================
# Downsampling
# Balance the total number of cells across WT, CD169DTR
# and CCR2KO. The smallest total is rounded down to the
# nearest multiple of 10.
# ============================================================

groups_plot <- c(
  "WT",
  "CD169DTR",
  "CCR2KO"
)

metadata <- Tcells[[]]
metadata$Cell <- rownames(metadata)

total_cells_by_group <- table(
  factor(
    metadata$orig.ident,
    levels = groups_plot
  )
)

if (any(total_cells_by_group == 0)) {
  stop(
    "At least one orig.ident group contains no cells after UMAP filtering: ",
    paste(
      names(total_cells_by_group)[total_cells_by_group == 0],
      collapse = ", "
    )
  )
}

minimum_total <- min(total_cells_by_group)

downsample_n <- floor(
  minimum_total / 10
) * 10L

if (downsample_n < 10) {
  stop(
    "The smallest group contains fewer than 10 cells after UMAP filtering."
  )
}

cat(
  "\nTotal cells before downsampling:\n"
)

print(total_cells_by_group)

cat(
  "\nDownsampling each group to: ",
  downsample_n,
  " cells\n",
  sep = ""
)

set.seed(1234)

cells_downsampled <- unlist(
  lapply(
    groups_plot,
    function(current_group) {

      current_cells <- metadata$Cell[
        as.character(metadata$orig.ident) == current_group
      ]

      sample(
        current_cells,
        size = downsample_n,
        replace = FALSE
      )
    }
  ),
  use.names = FALSE
)

Tcells_plot <- subset(
  Tcells,
  cells = cells_downsampled
)

tcell_order_plot <- tcell_order[
  tcell_order %in%
    unique(
      as.character(
        Tcells_plot$Figure3_annotation
      )
    )
]

Tcells_plot$Figure3_annotation <- factor(
  as.character(Tcells_plot$Figure3_annotation),
  levels = tcell_order_plot
)

Tcells_plot$orig.ident <- factor(
  as.character(Tcells_plot$orig.ident),
  levels = groups_plot
)

tcell_colors_plot <- tcell_colors[
  tcell_order_plot
]

Idents(Tcells_plot) <- "Figure3_annotation"

cat(
  "\nCells after UMAP filtering: ",
  ncol(Tcells),
  "\nCells after total-group downsampling: ",
  ncol(Tcells_plot),
  "\n",
  sep = ""
)

cat(
  "\nTotal cells after downsampling:\n"
)

print(
  table(
    Tcells_plot$orig.ident
  )
)

cat(
  "\nCell-type composition after downsampling:\n"
)

print(
  table(
    Tcells_plot$orig.ident,
    Tcells_plot$Figure3_annotation
  )
)


# ============================================================
# Figure 3E
# Global T / lymphoid UMAP
# ============================================================

Figure3E <- DimPlot(
  object = Tcells_plot,
  reduction = "umap",
  group.by = "Figure3_annotation",
  label = FALSE,
  pt.size = 0.6,
  raster = FALSE
) +
  scale_color_manual(
    values = tcell_colors_plot,
    breaks = tcell_order_plot,
    drop = FALSE
  ) +
  labs(
    title = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 4),
      ncol = 1
    )
  )

ggsave(
  filename = file.path(
    figure_dir,
    "Figure3E_Tcells_annotation.pdf"
  ),
  plot = Figure3E,
  width = 14,
  height = 9,
  units = "in",
  useDingbats = FALSE
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure3E_Tcells_annotation.png"
  ),
  plot = Figure3E,
  width = 14,
  height = 9,
  units = "in",
  dpi = 300
)

# ============================================================
# Figure 3F
# T / lymphoid UMAP split by orig.ident
# ============================================================

plot_list <- DimPlot(
  object = Tcells_plot,
  reduction = "umap",
  group.by = "Figure3_annotation",
  split.by = "orig.ident",
  label = FALSE,
  pt.size = 0.6,
  combine = FALSE,
  raster = FALSE
)

plot_list <- lapply(
  plot_list,
  function(p) {
    p +
      scale_color_manual(
        values = tcell_colors_plot,
        breaks = tcell_order_plot,
        drop = FALSE
      ) +
      labs(
        title = NULL,
        color = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        legend.title = element_blank(),
        legend.text = element_text(size = 9)
      ) +
      guides(
        color = guide_legend(
          override.aes = list(size = 4),
          ncol = 1
        )
      )
  }
)

Figure3F <- wrap_plots(
  plot_list,
  nrow = 1,
  guides = "collect"
) &
  theme(
    legend.position = "right",
    legend.title = element_blank()
  )

ggsave(
  filename = file.path(
    figure_dir,
    "Figure3F_Tcells_annotation_split_orig.ident.pdf"
  ),
  plot = Figure3F,
  width = 28,
  height = 8,
  units = "in",
  useDingbats = FALSE,
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure3F_Tcells_annotation_split_orig.ident.png"
  ),
  plot = Figure3F,
  width = 28,
  height = 8,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

# ============================================================
# Source table
# ============================================================

umap_plot <- Embeddings(
  Tcells_plot,
  reduction = "umap"
)

metadata_plot <- Tcells_plot[[]]

source_table <- data.frame(
  Cell = rownames(umap_plot),
  UMAP_1 = umap_plot[, 1],
  UMAP_2 = umap_plot[, 2],
  orig.ident = as.character(
    metadata_plot[
      rownames(umap_plot),
      "orig.ident"
    ]
  ),
  group = as.character(
    metadata_plot[
      rownames(umap_plot),
      "group"
    ]
  ),
  Clustering_original = as.character(
    metadata_plot[
      rownames(umap_plot),
      "Clustering"
    ]
  ),
  Figure3_annotation = as.character(
    metadata_plot[
      rownames(umap_plot),
      "Figure3_annotation"
    ]
  ),
  stringsAsFactors = FALSE
)

write.table(
  source_table,
  file = file.path(
    table_dir,
    "Figure3EF_Tcells_umap_downsampled.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Downsampling summary
# ============================================================

counts_before <- as.data.frame(
  table(
    orig.ident = Tcells$orig.ident,
    Cell_Type = Tcells$Figure3_annotation
  )
)

counts_after <- as.data.frame(
  table(
    orig.ident = Tcells_plot$orig.ident,
    Cell_Type = Tcells_plot$Figure3_annotation
  )
)

colnames(counts_before)[3] <- "Cells_before"
colnames(counts_after)[3] <- "Cells_after"

downsampling_summary <- merge(
  counts_before,
  counts_after,
  by = c(
    "orig.ident",
    "Cell_Type"
  ),
  all = TRUE
)

downsampling_summary[
  is.na(downsampling_summary)
] <- 0

write.table(
  downsampling_summary,
  file = file.path(
    table_dir,
    "Figure3EF_downsampling_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    table_dir,
    "Figure3_sessionInfo.txt"
  )
)

cat(
  "\nFigure 3E-F completed.\n",
  "Figures: ",
  figure_dir,
  "\nTables: ",
  table_dir,
  "\n",
  sep = ""
)


# ============================================================
# Figure 3G
# T / lymphoid composition
# Same downsampled cells used in Figure 3E-F
# ============================================================

obj_quant <- Tcells_plot

groups_quant <- c(
  "WT",
  "CD169DTR",
  "CCR2KO"
)

types_quant <- c(
  "Naive",
  "Tumor"
)

cluster_order_quant <- tcell_order_plot[
  tcell_order_plot %in%
    unique(as.character(obj_quant$Figure3_annotation))
]

# Same palette used in the UMAPs
nora_colors <- tcell_colors_plot[
  cluster_order_quant
]

missing_colors <- setdiff(
  cluster_order_quant,
  names(nora_colors)
)

if (length(missing_colors) > 0) {
  stop(
    "Missing colors for: ",
    paste(missing_colors, collapse = ", ")
  )
}

# ============================================================
# Frequency and proportion table
# ============================================================

table_tcells <- obj_quant[[]] |>
  dplyr::mutate(
    orig.ident = as.character(orig.ident),
    group = as.character(group),
    Figure3_annotation = trimws(
      as.character(Figure3_annotation)
    )
  ) |>
  dplyr::filter(
    orig.ident %in% groups_quant,
    group %in% types_quant,
    Figure3_annotation %in% cluster_order_quant
  ) |>
  dplyr::count(
    orig.ident,
    group,
    Figure3_annotation,
    name = "Freq"
  ) |>
  tidyr::complete(
    orig.ident = groups_quant,
    group = types_quant,
    Figure3_annotation = cluster_order_quant,
    fill = list(
      Freq = 0
    )
  ) |>
  dplyr::group_by(
    orig.ident,
    group
  ) |>
  dplyr::mutate(
    Total_cells = sum(Freq),
    Proportion = ifelse(
      Total_cells > 0,
      Freq / Total_cells,
      0
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::rename(
    Group = orig.ident,
    Type = group,
    Clusters = Figure3_annotation
  )

table_tcells$Group <- factor(
  table_tcells$Group,
  levels = groups_quant
)

table_tcells$Type <- factor(
  table_tcells$Type,
  levels = types_quant
)

table_tcells$Clusters <- factor(
  table_tcells$Clusters,
  levels = cluster_order_quant
)

cat(
  "\nT-cell composition calculated using the ",
  "same downsampled cells as Figure 3E-F:\n",
  sep = ""
)

print(table_tcells)

# ============================================================
# Save source-data table
# ============================================================

write.table(
  table_tcells,
  file = file.path(
    table_dir,
    "Figure3G_Tcells_composition_downsampled.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Stacked proportion plot
# ============================================================

Figure3G <- ggplot(
  table_tcells,
  aes(
    x = Type,
    y = Proportion,
    fill = Clusters
  )
) +
  geom_col(
    width = 0.75
  ) +
  facet_wrap(
    ~Group,
    nrow = 1
  ) +
  scale_fill_manual(
    values = nora_colors,
    breaks = cluster_order_quant,
    drop = FALSE,
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    x = NULL,
    y = "Proportion",
    fill = NULL
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    strip.background = element_rect(
      fill = "white",
      color = "black"
    ),
    axis.text.x = element_text(
      size = 10
    ),
    axis.text.y = element_text(
      size = 10
    ),
    axis.title.y = element_text(
      size = 11
    ),
    legend.title = element_blank(),
    legend.text = element_text(
      size = 8
    ),
    legend.position = "right"
  ) +
  guides(
    fill = guide_legend(
      title = NULL,
      ncol = 1
    )
  )

# ============================================================
# Save Figure 3G
# ============================================================

ggsave(
  filename = file.path(
    figure_dir,
    "Figure3G_Tcells_composition_downsampled.pdf"
  ),
  plot = Figure3G,
  width = 12,
  height = 8,
  units = "in",
  useDingbats = FALSE
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure3G_Tcells_composition_downsampled.png"
  ),
  plot = Figure3G,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)

cat(
  "\nFigure 3G completed.\n",
  "Figure: ",
  file.path(
    figure_dir,
    "Figure3G_Tcells_composition_downsampled.pdf"
  ),
  "\nTable: ",
  file.path(
    table_dir,
    "Figure3G_Tcells_composition_downsampled.tsv"
  ),
  "\n",
  sep = ""
)