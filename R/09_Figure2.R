# ============================================================
# Figure 2E-F
# Fibroblast trajectory and pseudotime heatmap
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(dendextend)
  library(grid)
})

set.seed(1234)

# ============================================================
# Configuration
# ============================================================

read_config <- function(path = "config.tsv") {
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x) & !startsWith(x, "#")]

  parts <- strsplit(x, "[\t=]+")

  keys <- vapply(
    parts,
    `[`,
    character(1),
    1
  )

  values <- vapply(
    parts,
    function(z) paste(z[-1], collapse = "="),
    character(1)
  )

  setNames(
    trimws(values),
    trimws(keys)
  )
}

config <- read_config()

rsdir <- unname(config["rsdir"])
outdir <- unname(config["outdir"])

objects_dir <- file.path(
  rsdir,
  "objects"
)

figure_dir <- file.path(
  outdir,
  "Figure2"
)

table_dir <- file.path(
  rsdir,
  "tables",
  "Figure2"
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

n_cores <- 4L

# ============================================================
# Load final annotated object
# ============================================================

fibroblasts <- readRDS(
  file.path(
    objects_dir,
    "fibroblasts.Not.Integrated.annotated.rds"
  )
)

DefaultAssay(fibroblasts) <- "RNA"

# ============================================================
# Fibroblast subset
# Mesothelial cells are excluded before Monocle
# ============================================================

fibro_states <- c(
  "Fibroblasts Universal",
  "Fibroblasts Activated 1",
  "Fibroblasts Activated 2"
)

metadata <- fibroblasts[[]]

cells_keep <- rownames(metadata)[
  as.character(
    metadata$Clustering.Mesenchymal
  ) == "Mesenchymal" &
    as.character(
      metadata$Clustering
    ) %in% fibro_states
]

fibro_subset <- subset(
  fibroblasts,
  cells = cells_keep
)

fibro_subset$Clustering <- factor(
  as.character(
    fibro_subset$Clustering
  ),
  levels = fibro_states
)

Idents(fibro_subset) <- "Clustering"
DefaultAssay(fibro_subset) <- "RNA"

cat("\nCélulas fibroblásticas utilizadas:\n")
print(table(fibro_subset$Clustering))

# ============================================================
# Existing RNA layers
# No Seurat renormalization
# ============================================================

fibro_subset <- JoinLayers(
  fibro_subset,
  assay = "RNA"
)

counts_matrix <- GetAssayData(
  fibro_subset,
  assay = "RNA",
  layer = "counts"
)

normalized_matrix <- GetAssayData(
  fibro_subset,
  assay = "RNA",
  layer = "data"
)

# ============================================================
# Create Monocle 3 object
# ============================================================

cell_metadata <- fibro_subset[[]]

cell_metadata <- cell_metadata[
  colnames(counts_matrix),
  ,
  drop = FALSE
]

gene_metadata <- data.frame(
  gene_short_name = rownames(counts_matrix),
  row.names = rownames(counts_matrix),
  stringsAsFactors = FALSE
)

state_labels <- c(
  "Fibroblasts Universal" = "Universal",
  "Fibroblasts Activated 1" = "Activated fibroblast 1",
  "Fibroblasts Activated 2" = "Activated fibroblast 2"
)

cell_metadata$Figure2_state <- factor(
  unname(
    state_labels[
      as.character(
        cell_metadata$Clustering
      )
    ]
  ),
  levels = unname(
    state_labels[
      fibro_states
    ]
  )
)

cds <- new_cell_data_set(
  expression_data = counts_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

cds <- estimate_size_factors(cds)

# ============================================================
# Monocle-specific PCA and UMAP
# ============================================================

cat("\nEjecutando preprocess_cds...\n")

cds <- preprocess_cds(
  cds,
  num_dim = 50
)

set.seed(1234)

cat("Calculando UMAP de Monocle...\n")

cds <- reduce_dimension(
  cds,
  reduction_method = "UMAP",
  umap.min_dist = 0.4,
  umap.n_neighbors = 30,
  cores = n_cores
)

set.seed(1234)

cds <- cluster_cells(
  cds,
  reduction_method = "UMAP",
  resolution = 1e-2
)

cat("Aprendiendo principal graph...\n")

cds <- learn_graph(
  cds,
  use_partition = TRUE,
  learn_graph_control = list(
    minimal_branch_len = 40,
    prune_graph = FALSE
  )
)

# ============================================================
# Pi16-high Universal fibroblasts as root candidates
# ============================================================

pi16_values <- as.numeric(
  normalized_matrix[
    "Pi16",
    colnames(cds),
    drop = TRUE
  ]
)

names(pi16_values) <- colnames(cds)

colData(cds)$Pi16_normalized <- pi16_values

universal_cells <- colnames(cds)[
  as.character(
    colData(cds)$Clustering
  ) == "Fibroblasts Universal"
]

n_root_candidates <- max(
  20L,
  ceiling(
    0.10 * length(universal_cells)
  )
)

n_root_candidates <- min(
  n_root_candidates,
  length(universal_cells)
)

universal_pi16 <- sort(
  pi16_values[
    universal_cells
  ],
  decreasing = TRUE,
  na.last = NA
)

root_candidate_cells <- names(
  head(
    universal_pi16,
    n_root_candidates
  )
)

colData(cds)$Root_candidate <- factor(
  ifelse(
    colnames(cds) %in% root_candidate_cells,
    "Pi16-high Universal",
    "Other cells"
  ),
  levels = c(
    "Other cells",
    "Pi16-high Universal"
  )
)

cat(
  "\nUniversal fibroblasts: ",
  length(universal_cells),
  "\nPi16-high root candidates: ",
  length(root_candidate_cells),
  "\n",
  sep = ""
)

# ============================================================
# Root-selection diagnostic
# ============================================================

fibro_colors <- c(
  "Fibroblasts Universal" = "#FBCDC7",
  "Fibroblasts Activated 1" = "#F37B79",
  "Fibroblasts Activated 2" = "#6D1A1A"
)

p_root_states <- plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = "Clustering",
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  show_trajectory_graph = TRUE,
  label_branch_points = TRUE,
  label_leaves = TRUE,
  label_roots = FALSE,
  graph_label_size = 2,
  cell_size = 0.6
) +
  scale_color_manual(
    values = fibro_colors
  ) +
  labs(
    title = "Fibroblast states",
    x = "UMAP1",
    y = "UMAP2"
  )

p_root_pi16 <- plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = "Pi16_normalized",
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  show_trajectory_graph = TRUE,
  label_branch_points = TRUE,
  label_leaves = TRUE,
  label_roots = FALSE,
  graph_label_size = 2,
  cell_size = 0.6
) +
  scale_color_viridis_c(
    option = "magma",
    name = "Pi16"
  ) +
  labs(
    title = "Pi16 expression",
    x = "UMAP1",
    y = "UMAP2"
  )

p_root_candidates <- plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = "Root_candidate",
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  show_trajectory_graph = TRUE,
  label_branch_points = TRUE,
  label_leaves = TRUE,
  label_roots = FALSE,
  graph_label_size = 2,
  cell_size = 0.6
) +
  scale_color_manual(
    values = c(
      "Other cells" = "grey85",
      "Pi16-high Universal" = "#0072B2"
    )
  ) +
  labs(
    title = "Proposed root cells",
    x = "UMAP1",
    y = "UMAP2"
  )

Figure2E_root_diagnostic <-
  p_root_states +
  p_root_pi16 +
  p_root_candidates +
  plot_layout(
    ncol = 3
  )

ggsave(
  filename = file.path(
    figure_dir,
    "Figure2E_root_selection_diagnostic.pdf"
  ),
  plot = Figure2E_root_diagnostic,
  width = 12,
  height = 4,
  units = "in"
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure2E_root_selection_diagnostic.png"
  ),
  plot = Figure2E_root_diagnostic,
  width = 12,
  height = 4,
  units = "in",
  dpi = 300
)

# ============================================================
# Map Pi16-high cells to the principal graph
# ============================================================

closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$
  pr_graph_cell_proj_closest_vertex

closest_vertex <- as.matrix(
  closest_vertex
)

if (is.null(rownames(closest_vertex))) {
  rownames(closest_vertex) <- colnames(cds)
}

closest_vertex <- closest_vertex[
  colnames(cds),
  ,
  drop = FALSE
]

candidate_vertices <- as.character(
  closest_vertex[
    root_candidate_cells,
    1
  ]
)

dominant_vertex <- names(
  sort(
    table(candidate_vertices),
    decreasing = TRUE
  )
)[1]

graph_vertex_names <- igraph::V(
  principal_graph(cds)[["UMAP"]]
)$name

if (dominant_vertex %in% graph_vertex_names) {
  root_principal_node <- dominant_vertex
} else {
  root_principal_node <- graph_vertex_names[
    as.integer(dominant_vertex)
  ]
}

cat(
  "\nNodo raíz seleccionado: ",
  root_principal_node,
  "\n",
  sep = ""
)

# ============================================================
# Calculate pseudotime
# ============================================================

cds <- order_cells(
  cds,
  reduction_method = "UMAP",
  root_pr_nodes = root_principal_node
)

pseudotime_values <- pseudotime(
  cds,
  reduction_method = "UMAP"
)

names(pseudotime_values) <- colnames(cds)

colData(cds)$monocle3_pseudotime <-
  pseudotime_values

cat(
  "Células con pseudotime infinito: ",
  sum(
    !is.finite(
      pseudotime_values
    )
  ),
  "\n",
  sep = ""
)

suppressWarnings(
  saveRDS(
    cds,
    file.path(
      objects_dir,
      "Figure2_fibroblast_monocle3_cds.rds"
    )
  )
)



# ============================================================
# Figure 2E
# Fibroblast states and pseudotime
# ============================================================

fibro_colors <- c(
  "Fibroblasts Universal" = "#FBCDC7",
  "Fibroblasts Activated 1" = "#F37B79",
  "Fibroblasts Activated 2" = "#6D1A1A"
)

# ============================================================
# Copy Monocle UMAP to Seurat
# ============================================================

monocle_umap <- SingleCellExperiment::reducedDim(
  cds,
  "UMAP"
)

rownames(monocle_umap) <- colnames(cds)

monocle_umap <- monocle_umap[
  colnames(fibro_subset),
  ,
  drop = FALSE
]

colnames(monocle_umap) <- c(
  "mUMAP_1",
  "mUMAP_2"
)

fibro_subset[["monocle_umap"]] <- CreateDimReducObject(
  embeddings = monocle_umap,
  key = "mUMAP_",
  assay = "RNA"
)

# ============================================================
# Left panel: clustering without trajectories
# ============================================================

Figure2E_clusters <- DimPlot(
  object = fibro_subset,
  reduction = "monocle_umap",
  group.by = "Clustering",
  cols = fibro_colors,
  pt.size = 0.7,
  raster = FALSE
) +
  labs(
    title = "Fibroblast states",
    x = "UMAP1",
    y = "UMAP2",
    color = NULL
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    legend.title = element_blank()
  )

# ============================================================
# Right panel: pseudotime with trajectories
# ============================================================

Figure2E_pseudotime <- plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "black",
  trajectory_graph_segment_size = 0.6,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_leaves = FALSE,
  cell_size = 0.7
) +
  scale_color_viridis_c(
    option = "plasma",
    name = "Pseudotime"
  ) +
  labs(
    title = "Pseudotime",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1
  )

# ============================================================
# Combined figure
# ============================================================

Figure2E <- Figure2E_clusters + Figure2E_pseudotime +
  plot_layout(
    ncol = 2,
    widths = c(1, 1)
  )

ggsave(
  filename = file.path(
    figure_dir,
    "Figure2E_fibroblast_states_pseudotime.pdf"
  ),
  plot = Figure2E,
  width = 9,
  height = 4.2,
  units = "in",
  useDingbats = FALSE
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure2E_fibroblast_states_pseudotime.png"
  ),
  plot = Figure2E,
  width = 9,
  height = 4.2,
  units = "in",
  dpi = 300
)

# ============================================================
# Figure 2E source table
# ============================================================

figure2e_cells <- colnames(cds)

Figure2E_table <- data.frame(
  Cell = figure2e_cells,
  UMAP1 = monocle_umap[figure2e_cells, 1],
  UMAP2 = monocle_umap[figure2e_cells, 2],
  Clustering = as.character(
    colData(cds)$Clustering
  ),
  Pseudotime = monocle3::pseudotime(cds)[
    figure2e_cells
  ],
  stringsAsFactors = FALSE
)

write.table(
  Figure2E_table,
  file = file.path(
    table_dir,
    "Figure2E_fibroblast_states_pseudotime.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# Figure 2F: graph_test
# ============================================================

cat("\nEjecutando graph_test...\n")

modulated_genes <- graph_test(
  cds,
  neighbor_graph = "principal_graph",
  cores = n_cores
)

graph_test_table <- as.data.frame(
  modulated_genes
)

graph_test_table$Gene <- rownames(
  graph_test_table
)

write.table(
  graph_test_table,
  file.path(
    table_dir,
    "Figure2F_graph_test_all.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sig_genes <- graph_test_table[
  graph_test_table$status == "OK" &
    graph_test_table$q_value < 0.05 &
    is.finite(
      graph_test_table$morans_I
    ),
  ,
  drop = FALSE
]

sig_genes <- sig_genes[
  order(
    sig_genes$morans_I,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

top_genes <- head(
  sig_genes,
  60
)

valid_genes <- intersect(
  top_genes$Gene,
  rownames(
    rowData(cds)
  )
)

valid_genes <- setdiff(
  valid_genes,
  c(
    "Mapkapk3",
    "Dpep1"
  )
)

cat(
  "Genes seleccionados: ",
  length(valid_genes),
  "\n",
  sep = ""
)

# ============================================================
# Expression ordered by pseudotime
# ============================================================

finite_cells <- names(
  pseudotime_values
)[
  is.finite(
    pseudotime_values
  )
]

ordered_cells <- finite_cells[
  order(
    pseudotime_values[
      finite_cells
    ]
  )
]

ordered_pseudotime <- pseudotime_values[
  ordered_cells
]

pt_matrix <- exprs(cds)[
  valid_genes,
  ordered_cells,
  drop = FALSE
]

# ============================================================
# Smooth expression along ordered cell position
# ============================================================

cat("Aplicando smooth.spline(df = 3)...\n")

pt_matrix <- t(
  apply(
    pt_matrix,
    1,
    function(x) {
      x <- as.numeric(x)

      if (
        length(unique(x)) < 2 ||
          sd(x) == 0
      ) {
        return(
          rep(
            0,
            length(x)
          )
        )
      }

      smooth.spline(
        x = seq_along(x),
        y = x,
        df = 3
      )$y
    }
  )
)

rownames(pt_matrix) <- valid_genes
colnames(pt_matrix) <- ordered_cells

# ============================================================
# Row Z-score
# ============================================================

pt_matrix <- t(
  apply(
    pt_matrix,
    1,
    function(x) {
      current_sd <- sd(x)

      if (
        !is.finite(current_sd) ||
          current_sd == 0
      ) {
        return(
          rep(
            0,
            length(x)
          )
        )
      }

      (
        x - mean(x)
      ) / current_sd
    }
  )
)

pt_matrix <- pt_matrix[
  complete.cases(pt_matrix),
  ,
  drop = FALSE
]

cat(
  "Genes finales en el heatmap: ",
  nrow(pt_matrix),
  "\n",
  sep = ""
)

# ============================================================
# Heatmap colors and column annotations
# ============================================================

expression_color <- colorRamp2(
  seq(
    -2,
    2,
    length.out = 11
  ),
  rev(
    brewer.pal(
      11,
      "Spectral"
    )
  )
)

group_colors <- c(
  "Naive" = "#FA5858",
  "Tumor" = "#04B486"
)

sample_colors <- c(
  "WT" = "#A6CEE3",
  "CD169DTR" = "#1F78B4",
  "CCR2KO" = "#08306B"
)

heatmap_metadata <- as.data.frame(
  colData(cds)
)[
  ordered_cells,
  ,
  drop = FALSE
]

cell_clusters <- factor(
  as.character(
    heatmap_metadata$Clustering
  ),
  levels = fibro_states
)

cell_group <- factor(
  as.character(
    heatmap_metadata$group
  ),
  levels = c(
    "Naive",
    "Tumor"
  )
)

cell_sample <- factor(
  as.character(
    heatmap_metadata$orig.ident
  ),
  levels = c(
    "WT",
    "CD169DTR",
    "CCR2KO"
  )
)

pseudotime_color <- colorRamp2(
  range(
    ordered_pseudotime,
    finite = TRUE
  ),
  c(
    "#F7FCF5",
    "#238B45"
  )
)

column_annotation <- HeatmapAnnotation(
  Cluster = cell_clusters,
  Group = cell_group,
  Sample = cell_sample,
  Pseudotime = ordered_pseudotime,
  col = list(
    Cluster = fibro_colors,
    Group = group_colors,
    Sample = sample_colors,
    Pseudotime = pseudotime_color
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  simple_anno_size = unit(
    4,
    "mm"
  )
)
# ============================================================
# Hierarchical clustering and temporal module order
# ============================================================

calculate_activation_center <- function(
  expression,
  pseudotime
) {
  weights <- expression - min(
    expression,
    na.rm = TRUE
  )

  if (sum(weights, na.rm = TRUE) == 0) {
    return(
      mean(
        pseudotime,
        na.rm = TRUE
      )
    )
  }

  weighted.mean(
    pseudotime,
    weights,
    na.rm = TRUE
  )
}

gene_activation_center <- apply(
  pt_matrix,
  1,
  calculate_activation_center,
  pseudotime = ordered_pseudotime
)

# Clustering jerárquico global para definir tres módulos
row_hc <- hclust(
  dist(pt_matrix),
  method = "ward.D2"
)

row_clusters <- cutree(
  row_hc,
  k = 3
)

# Centro temporal de cada módulo
cluster_activation_center <- tapply(
  gene_activation_center,
  row_clusters,
  median,
  na.rm = TRUE
)

# Módulos ordenados naturalmente:
# menor pseudotime -> mayor pseudotime
temporal_cluster_order <- names(
  sort(
    cluster_activation_center,
    decreasing = FALSE
  )
)

cluster_to_phase <- setNames(
  c(
    "Early",
    "Intermediate",
    "Late"
  ),
  temporal_cluster_order
)

gene_phase <- cluster_to_phase[
  as.character(
    row_clusters
  )
]

names(gene_phase) <- names(
  row_clusters
)

cat("\nAsignación temporal de módulos:\n")
print(
  data.frame(
    Cluster = temporal_cluster_order,
    Phase = c(
      "Early",
      "Intermediate",
      "Late"
    ),
    Activation_center = cluster_activation_center[
      temporal_cluster_order
    ]
  )
)

# ============================================================
# Hierarchical order within each module
# ============================================================

get_hierarchical_order <- function(genes) {
  if (length(genes) <= 1) {
    return(genes)
  }

  module_hc <- hclust(
    dist(
      pt_matrix[
        genes,
        ,
        drop = FALSE
      ]
    ),
    method = "ward.D2"
  )

  module_hc$labels[
    module_hc$order
  ]
}

late_genes <- names(
  gene_phase
)[
  gene_phase == "Late"
]

intermediate_genes <- names(
  gene_phase
)[
  gene_phase == "Intermediate"
]

early_genes <- names(
  gene_phase
)[
  gene_phase == "Early"
]

# Orden jerárquico dentro de cada módulo
late_order <- get_hierarchical_order(
  late_genes
)

intermediate_order <- get_hierarchical_order(
  intermediate_genes
)

early_order <- get_hierarchical_order(
  early_genes
)

# Invertir solamente el módulo superior
late_order <- rev(
  late_order
)

# ============================================================
# Final visual order
# Top: Late
# Middle: Intermediate
# Bottom: Early
# ============================================================

final_gene_order <- c(
  early_order,
  intermediate_order,
  late_order
)

# Reordenar físicamente la matriz
pt_matrix_final <- pt_matrix[
  final_gene_order,
  ,
  drop = FALSE
]

gene_phase_final <- gene_phase[
  final_gene_order
]

row_clusters_final <- row_clusters[
  final_gene_order
]

gene_activation_center_final <-
  gene_activation_center[
    final_gene_order
  ]

cat("\nOrden visual final:\n")

print(
  table(
    factor(
      gene_phase_final,
      levels = c(
        "Late",
        "Intermediate",
        "Early"
      )
    )
  )
)

# ============================================================
# Figure 2F heatmap
# Hierarchical clustering used within modules
# No spaces between modules
# ============================================================

Figure2F <- Heatmap(
  pt_matrix_final,
  name = "Scaled Expression",
  col = expression_color,

  # El orden jerárquico ya está incorporado
  # físicamente en pt_matrix_final
  cluster_rows = FALSE,
  cluster_columns = FALSE,

  show_row_names = TRUE,
  show_column_names = FALSE,

  row_names_gp = gpar(
    fontsize = 9
  ),

  top_annotation = column_annotation,

  column_title = "Pseudotime",
  column_title_gp = gpar(
    fontsize = 16,
    fontface = "bold"
  ),

  use_raster = TRUE,
  raster_quality = 3
)

# ============================================================
# Save Figure 2F
# ============================================================

pdf(
  file = file.path(
    figure_dir,
    "Figure2F_pseudotime_heatmap.pdf"
  ),
  width = 13,
  height = 8,
  useDingbats = FALSE
)

draw(
  Figure2F,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

png(
  filename = file.path(
    figure_dir,
    "Figure2F_pseudotime_heatmap.png"
  ),
  width = 3600,
  height = 2200,
  res = 250
)

draw(
  Figure2F,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

# ============================================================
# Figure 2F source tables
# ============================================================

selected_gene_table <- top_genes[
  match(
    final_gene_order,
    top_genes$Gene
  ),
  ,
  drop = FALSE
]

selected_gene_table$Heatmap_cluster <- unname(
  row_clusters_final[
    selected_gene_table$Gene
  ]
)

selected_gene_table$Activation_center <- unname(
  gene_activation_center_final[
    selected_gene_table$Gene
  ]
)

selected_gene_table$Temporal_phase <- unname(
  gene_phase_final[
    selected_gene_table$Gene
  ]
)

selected_gene_table$Heatmap_order <- seq_len(
  nrow(selected_gene_table)
)

write.table(
  selected_gene_table,
  file = file.path(
    table_dir,
    "Figure2F_selected_genes.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

scaled_expression_table <- data.frame(
  Gene = rownames(
    pt_matrix_final
  ),
  pt_matrix_final,
  check.names = FALSE
)

write.table(
  scaled_expression_table,
  file = file.path(
    table_dir,
    "Figure2F_scaled_expression.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    table_dir,
    "Figure2_sessionInfo.txt"
  )
)

cat(
  "\nFigure 2F finalizada.\n",
  "Orden superior: Late, invertido internamente\n",
  "Orden central: Intermediate\n",
  "Orden inferior: Early\n",
  "PDF: ",
  file.path(
    figure_dir,
    "Figure2F_pseudotime_heatmap.pdf"
  ),
  "\nPNG: ",
  file.path(
    figure_dir,
    "Figure2F_pseudotime_heatmap.png"
  ),
  "\n",
  sep = ""
)





# ============================================================
# Figure 2G
# Cytokines, collagens and ECM heatmaps
# ============================================================

figure2g_states <- c(
  "Fibroblasts Universal",
  "Fibroblasts Activated 1",
  "Fibroblasts Activated 2",
  "Mesothelial"
)

meta <- fibroblasts[[]]

figure2g_cells <- rownames(meta)[
  as.character(meta$Clustering.Mesenchymal) == "Mesenchymal" &
    as.character(meta$Clustering) %in% figure2g_states
]

fibro_meso <- subset(
  fibroblasts,
  cells = figure2g_cells
)

fibro_meso$Clustering <- factor(
  as.character(fibro_meso$Clustering),
  levels = figure2g_states
)

DefaultAssay(fibro_meso) <- "RNA"

fibro_meso <- JoinLayers(
  fibro_meso,
  assay = "RNA"
)

cat("\nCélulas para los heatmaps de Figure 2G:\n")
print(table(fibro_meso$Clustering))

# ============================================================
# Gene sets
# ============================================================

cytokine_genes <- c(
  "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl7", "Ccl8", "Ccl11",
  "Ccl12", "Ccl17", "Ccl19", "Ccl20", "Ccl21", "Ccl22",
  "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl5", "Cxcl9", "Cxcl10",
  "Cxcl11", "Cxcl12", "Cxcl13", "Cxcl14", "Cxcl16",
  "Il1a", "Il1b", "Il6", "Il7", "Il11", "Il15", "Il18",
  "Il33", "Il34",
  "Tnf", "Tnfsf10", "Tnfsf11", "Tnfsf12", "Tnfsf13b",
  "Tnfsf14", "Lta",
  "Tgfb1", "Tgfb2", "Tgfb3",
  "Lif", "Osm", "Csf1", "Csf2", "Csf3", "Kitl", "Flt3l",
  "Areg", "Ereg", "Gdf15", "Vegfa", "Ptgs2"
)

collagen_genes <- c(
  "Col1a1", "Col1a2", "Col2a1", "Col3a1",
  "Col4a1", "Col4a2", "Col4a3", "Col4a4", "Col4a5", "Col4a6",
  "Col5a1", "Col5a2", "Col5a3",
  "Col6a1", "Col6a2", "Col6a3", "Col6a4", "Col6a5", "Col6a6",
  "Col7a1", "Col8a1", "Col8a2",
  "Col9a1", "Col9a2", "Col9a3",
  "Col10a1", "Col11a1", "Col11a2",
  "Col12a1", "Col13a1", "Col14a1", "Col15a1", "Col16a1",
  "Col17a1", "Col18a1", "Col19a1", "Col20a1", "Col21a1",
  "Col22a1", "Col23a1", "Col24a1", "Col25a1", "Col26a1",
  "Col27a1", "Col28a1"
)

ecm_genes <- c(
  "Fn1", "Fbn1", "Fbn2", "Postn", "Tnc", "Thbs1", "Thbs2",
  "Sparc", "Sparcl1", "Vtn", "Nid1", "Nid2", "Emilin1", "Emilin2",
  "Mfap2", "Mfap4", "Mfap5", "Eln", "Efemp1", "Efemp2",
  "Fbln1", "Fbln2", "Fbln5", "Fbln7",
  "Dcn", "Lum", "Bgn", "Fmod", "Aspn", "Ogn", "Prelp",
  "Hspg2", "Sdc1", "Sdc2", "Sdc3", "Sdc4",
  "Gpc1", "Gpc3", "Gpc4", "Gpc6",
  "Ccn1", "Ccn2", "Ccn3", "Ccn4", "Ccn5",
  "Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4",
  "Serpine1", "Serpinh1", "Pcolce", "Pcolce2",
  "Plod1", "Plod2", "Plod3",
  "P4ha1", "P4ha2", "P4ha3",
  "Timp1", "Timp2", "Timp3", "Timp4",
  "Acta2", "Tagln", "Myl9", "Myh11", "Cnn1"
)

# ============================================================
# Heatmap function
# ============================================================

make_figure2g_heatmap <- function(
  object,
  genes,
  title,
  filename
) {
  genes <- intersect(
    genes,
    rownames(object)
  )

  if (length(genes) == 0) {
    stop("No hay genes presentes para: ", title)
  }

  avg_mat <- AverageExpression(
    object = object,
    assays = "RNA",
    features = genes,
    group.by = "Clustering",
    layer = "data",
    verbose = FALSE
  )[["RNA"]]

  avg_mat <- as.matrix(avg_mat)
  storage.mode(avg_mat) <- "double"

  # Ordenar los estados
  state_keys <- make.names(figure2g_states)
  column_keys <- make.names(colnames(avg_mat))

  column_index <- match(
    state_keys,
    column_keys
  )

  keep_states <- !is.na(column_index)

  avg_mat <- avg_mat[
    ,
    column_index[keep_states],
    drop = FALSE
  ]

  colnames(avg_mat) <- figure2g_states[
    keep_states
  ]

  # Eliminar genes sin expresión
  avg_mat <- avg_mat[
    rowSums(avg_mat > 0, na.rm = TRUE) > 0,
    ,
    drop = FALSE
  ]

  # Z-score por gen
  row_mean <- rowMeans(
    avg_mat,
    na.rm = TRUE
  )

  row_sd <- apply(
    avg_mat,
    1,
    sd,
    na.rm = TRUE
  )

  z_mat <- sweep(
    avg_mat,
    1,
    row_mean,
    "-"
  )

  z_mat <- sweep(
    z_mat,
    1,
    row_sd,
    "/"
  )

  z_mat[!is.finite(z_mat)] <- 0
  z_mat[z_mat > 2.5] <- 2.5
  z_mat[z_mat < -2.5] <- -2.5

  heatmap_matrix <- t(z_mat)

  ph <- pheatmap::pheatmap(
    heatmap_matrix,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    border_color = NA,
    fontsize_row = 10,
    fontsize_col = 8,
    angle_col = 45,
    color = colorRampPalette(
      c(
        "navy",
        "white",
        "firebrick3"
      )
    )(100),
    breaks = seq(
      -2.5,
      2.5,
      length.out = 101
    ),
    main = title,
    silent = TRUE
  )

  plot_width <- max(
    14,
    ncol(heatmap_matrix) * 0.32
  )

  # PDF
  pdf(
    file.path(
      figure_dir,
      paste0(filename, ".pdf")
    ),
    width = plot_width,
    height = 4.5,
    useDingbats = FALSE
  )

  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()

  # PNG
  png(
    file.path(
      figure_dir,
      paste0(filename, ".png")
    ),
    width = plot_width,
    height = 4.5,
    units = "in",
    res = 300
  )

  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()

  # Source table
  write.table(
    data.frame(
      Gene = rownames(z_mat),
      z_mat,
      check.names = FALSE
    ),
    file.path(
      table_dir,
      paste0(filename, "_zscore.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat(
    title,
    ": ",
    nrow(z_mat),
    " genes\n",
    sep = ""
  )
}

# ============================================================
# Generate the three heatmaps
# ============================================================

make_figure2g_heatmap(
  object = fibro_meso,
  genes = cytokine_genes,
  title = "Cytokines / chemokines",
  filename = "Figure2G_cytokines_chemokines"
)

make_figure2g_heatmap(
  object = fibro_meso,
  genes = collagen_genes,
  title = "Collagens",
  filename = "Figure2G_collagens"
)

make_figure2g_heatmap(
  object = fibro_meso,
  genes = ecm_genes,
  title = "ECM / matrix-associated genes",
  filename = "Figure2G_ECM"
)

cat("\nFigure 2G finalizada.\n")