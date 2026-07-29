# ============================================================
# Figure 2E-G
# Fibroblast trajectory, pseudotime heatmap and expression heatmaps
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
  x <- trimws(readLines(path, warn = FALSE))
  x <- x[nzchar(x) & !startsWith(x, "#")]
  parts <- strsplit(x, "[\t=]+")

  keys <- vapply(parts, `[`, character(1), 1)
  values <- vapply(
    parts,
    function(z) paste(z[-1], collapse = "="),
    character(1)
  )

  setNames(trimws(values), trimws(keys))
}

config <- read_config()
rsdir <- unname(config["rsdir"])
outdir <- unname(config["outdir"])

objects_dir <- file.path(rsdir, "objects")
figure_dir <- file.path(outdir, "Figure2")
table_dir <- file.path(rsdir, "tables", "Figure2")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

n_cores <- 4L

preorder_cds_file <- file.path(
  objects_dir,
  "Figure2_fibroblast_mesothelial_before_ordering.rds"
)

manual_root_nodes_file <- file.path(
  table_dir,
  "Figure2E_manual_root_pr_nodes.txt"
)

root_selection_summary_file <- file.path(
  table_dir,
  "Figure2E_root_selection_used.tsv"
)

# ============================================================
# Load annotated object
# ============================================================

fibroblasts <- readRDS(
  file.path(
    objects_dir,
    "fibroblasts.Not.Integrated.annotated.rds"
  )
)

DefaultAssay(fibroblasts) <- "RNA"

fibro_states <- c(
  "Fibroblasts Universal",
  "Fibroblasts Activated 1",
  "Fibroblasts Activated 2"
)

trajectory_states <- c(
  fibro_states,
  "Mesothelial"
)

trajectory_colors <- c(
  "Fibroblasts Universal" = "#FBCDC7",
  "Fibroblasts Activated 1" = "#F37B79",
  "Fibroblasts Activated 2" = "#6D1A1A",
  "Mesothelial" = "#FED579"
)

fibro_colors <- trajectory_colors[fibro_states]

# ============================================================
# Build trajectory input with fibroblasts + mesothelial cells
# ============================================================

meta <- fibroblasts[[]]
trajectory_cells <- rownames(meta)[
  as.character(meta$Clustering) %in% trajectory_states
]

trajectory_subset <- subset(
  fibroblasts,
  cells = trajectory_cells
)

trajectory_subset$Clustering <- factor(
  as.character(trajectory_subset$Clustering),
  levels = trajectory_states
)

Idents(trajectory_subset) <- "Clustering"
DefaultAssay(trajectory_subset) <- "RNA"

cat("\nCélulas utilizadas para construir la trayectoria:\n")
print(table(trajectory_subset$Clustering, useNA = "ifany"))

if (sum(as.character(trajectory_subset$Clustering) == "Mesothelial") == 0) {
  stop(
    "No se han incluido células Mesothelial. ",
    "Revisa la columna Clustering del objeto anotado."
  )
}

trajectory_subset <- JoinLayers(
  trajectory_subset,
  assay = "RNA"
)

counts_matrix <- GetAssayData(
  trajectory_subset,
  assay = "RNA",
  layer = "counts"
)

cell_metadata <- trajectory_subset[[]]
cell_metadata <- cell_metadata[colnames(counts_matrix), , drop = FALSE]

gene_metadata <- data.frame(
  gene_short_name = rownames(counts_matrix),
  row.names = rownames(counts_matrix),
  stringsAsFactors = FALSE
)

if (file.exists(preorder_cds_file)) {
  cat(
    "\nReutilizando el CDS previo a order_cells:\n",
    preorder_cds_file,
    "\n",
    sep = ""
  )

  cds_all <- readRDS(preorder_cds_file)

} else {
  cds_all <- new_cell_data_set(
    expression_data = counts_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )

  cds_all <- estimate_size_factors(cds_all)

  # ==========================================================
  # Monocle PCA, UMAP, clustering and principal graph
  # Fibroblasts and mesothelial cells remain present
  # ==========================================================

  cat("\nEjecutando preprocess_cds...\n")
  cds_all <- preprocess_cds(
    cds_all,
    num_dim = 50
  )

  set.seed(1234)
  cat("Calculando UMAP de Monocle...\n")
  cds_all <- reduce_dimension(
    cds_all,
    reduction_method = "UMAP",
    umap.min_dist = 0.4,
    umap.n_neighbors = 30,
    cores = 1L
  )

  set.seed(1234)
  cds_all <- cluster_cells(
    cds_all,
    reduction_method = "UMAP",
    resolution = 1e-2
  )

  cat("Aprendiendo principal graph...\n")
  cds_all <- learn_graph(
    cds_all,
    use_partition = TRUE,
    learn_graph_control = list(
      minimal_branch_len = 40,
      prune_graph = FALSE
    )
  )

  saveRDS(cds_all, preorder_cds_file)
  cat(
    "CDS previo a order_cells guardado en: ",
    preorder_cds_file,
    "\n",
    sep = ""
  )
}

cat("\nPoblaciones presentes en el CDS usado para seleccionar el root:\n")
print(table(colData(cds_all)$Clustering, useNA = "ifany"))

if (sum(as.character(colData(cds_all)$Clustering) == "Mesothelial") == 0) {
  stop(
    "El CDS previo a order_cells no contiene células Mesothelial. ",
    "Borra el RDS previo y vuelve a ejecutar el script."
  )
}

graph_vertex_names <- igraph::V(
  principal_graph(cds_all)[["UMAP"]]
)$name

# ============================================================
# Root selection and pseudotime
# ============================================================

if (!file.exists(manual_root_nodes_file)) {
  cat(
    "\nEl grafo con fibroblastos y mesoteliales ya está preparado.\n",
    "No existe todavía una selección manual del nodo raíz.\n\n",
    "Genera ahora el mapa estático de nodos:\n",
    "  Rscript R/09_Figure2_select_root_static.R\n\n",
    "Después de mirar el PNG, guarda el root con, por ejemplo:\n",
    "  Rscript R/09_Figure2_select_root_static.R Y_17\n\n",
    "Finalmente vuelve a ejecutar:\n",
    "  Rscript R/09_Figure2.R\n",
    sep = ""
  )

  quit(
    save = "no",
    status = 0
  )
}

root_nodes_used <- unique(
  trimws(
    readLines(
      manual_root_nodes_file,
      warn = FALSE
    )
  )
)
root_nodes_used <- root_nodes_used[nzchar(root_nodes_used)]

if (length(root_nodes_used) == 0) {
  stop(
    "El archivo de nodos raíz está vacío: ",
    manual_root_nodes_file
  )
}

invalid_nodes <- setdiff(root_nodes_used, graph_vertex_names)
if (length(invalid_nodes) > 0) {
  stop(
    "Los nodos manuales no existen en el grafo actual: ",
    paste(invalid_nodes, collapse = ", "),
    ". Borra la selección anterior y ejecuta de nuevo ",
    "R/09_Figure2_select_root_static.R."
  )
}

root_selection_method <- "Saved manual principal nodes"

cds_all <- order_cells(
  cds_all,
  reduction_method = "UMAP",
  root_pr_nodes = root_nodes_used
)

pseudotime_all <- pseudotime(
  cds_all,
  reduction_method = "UMAP"
)
names(pseudotime_all) <- colnames(cds_all)
colData(cds_all)$monocle3_pseudotime <- pseudotime_all

cat(
  "\nMétodo de selección del root: ", root_selection_method,
  "\nRoot node(s): ", paste(root_nodes_used, collapse = ", "),
  "\nCélulas con pseudotime infinito: ",
  sum(!is.finite(pseudotime_all)),
  "\n",
  sep = ""
)

write.table(
  data.frame(
    Selection_method = root_selection_method,
    Root_pr_node = root_nodes_used,
    stringsAsFactors = FALSE
  ),
  root_selection_summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  cds_all,
  file.path(
    objects_dir,
    "Figure2_fibroblast_mesothelial_monocle3_cds.rds"
  )
)

# ============================================================
# Remove mesothelial cells only after root selection, order_cells()
# and pseudotime calculation have been completed with all four states.
# No dimensional reduction, graph learning or ordering is repeated.
# ============================================================

# cds_all retains the complete trajectory used for:
# - PCA and UMAP
# - principal graph
# - manual root-node selection
# - order_cells()
# - pseudotime calculation

fibro_cells <- colnames(cds_all)[
  as.character(colData(cds_all)$Clustering) %in% fibro_states
]

# Monocle object used only for the final fibroblast figures and analyses.
# It inherits the UMAP, principal graph and pseudotime from cds_all.
cds_fibro <- cds_all[, fibro_cells]

pseudotime_fibro_raw <- pseudotime_all[fibro_cells]
names(pseudotime_fibro_raw) <- fibro_cells

# Rescale fibroblast pseudotime to a fixed 0-20 range.
# This changes only the numerical display scale, not cell ordering.
finite_pseudotime <- is.finite(pseudotime_fibro_raw)

if (!any(finite_pseudotime)) {
  stop("No hay células de fibroblastos con pseudotime finito.")
}

pseudotime_limits_raw <- range(
  pseudotime_fibro_raw[finite_pseudotime],
  finite = TRUE
)

if (diff(pseudotime_limits_raw) == 0) {
  stop("El pseudotime de fibroblastos no presenta variación.")
}

pseudotime_fibro <- rep(
  Inf,
  length(pseudotime_fibro_raw)
)

names(pseudotime_fibro) <- fibro_cells

pseudotime_fibro[finite_pseudotime] <- (
  pseudotime_fibro_raw[finite_pseudotime] -
    pseudotime_limits_raw[1]
) / diff(pseudotime_limits_raw) * 20

colData(cds_fibro)$monocle3_pseudotime_raw <-
  pseudotime_fibro_raw

colData(cds_fibro)$monocle3_pseudotime <-
  pseudotime_fibro

saveRDS(
  cds_fibro,
  file.path(
    objects_dir,
    "Figure2_fibroblast_monocle3_cds.rds"
  )
)

fibro_subset_plot <- subset(
  trajectory_subset,
  cells = fibro_cells
)

fibro_subset_plot$Clustering <- factor(
  as.character(fibro_subset_plot$Clustering),
  levels = fibro_states
)

Idents(fibro_subset_plot) <- "Clustering"

cat("\nCélulas usadas para calcular el pseudotime:\n")
print(table(colData(cds_all)$Clustering, useNA = "ifany"))

cat("\nCélulas conservadas para Figure 2E-F:\n")
print(table(colData(cds_fibro)$Clustering, useNA = "ifany"))

if (any(as.character(colData(cds_fibro)$Clustering) == "Mesothelial")) {
  stop("Mesothelial sigue presente en el CDS final de Figure 2E-F.")
}

# ============================================================
# Figure 2E: fibroblast states and fibroblast pseudotime
#
# The UMAP, principal graph, manual root selection and pseudotime
# were calculated with fibroblasts + mesothelial cells.
# For the final figure, only fibroblast cells and the graph branches
# needed to connect fibroblast-associated principal points are shown.
# ============================================================

# Adjustable dimensions for Figure 2E.
# Because coord_cartesian() is used instead of coord_fixed(), increasing
# Figure 2E: tamaño final y anchura relativa de los dos paneles.
# coord_fixed(ratio = 1) conserva la geometría original del UMAP.
figure2e_width <- 11.5
figure2e_height <- 5.2
figure2e_panel_widths <- c(1, 1)

# ------------------------------------------------------------
# Monocle UMAP coordinates and fibroblast plotting table
# ------------------------------------------------------------

monocle_umap_all <- SingleCellExperiment::reducedDim(
  cds_all,
  "UMAP"
)

rownames(monocle_umap_all) <- colnames(cds_all)
colnames(monocle_umap_all) <- c("UMAP1", "UMAP2")

fibro_pseudotime_table <- data.frame(
  Cell = fibro_cells,
  UMAP1 = monocle_umap_all[fibro_cells, "UMAP1"],
  UMAP2 = monocle_umap_all[fibro_cells, "UMAP2"],
  Clustering = factor(
    as.character(colData(cds_all)[fibro_cells, "Clustering"]),
    levels = fibro_states
  ),
  Pseudotime_raw = as.numeric(
    pseudotime_fibro_raw[fibro_cells]
  ),
  Pseudotime = as.numeric(
    pseudotime_fibro[fibro_cells]
  ),
  stringsAsFactors = FALSE
)

fibro_pseudotime_table$Pseudotime_raw[
  !is.finite(fibro_pseudotime_table$Pseudotime_raw)
] <- NA_real_

fibro_pseudotime_table$Pseudotime[
  !is.finite(fibro_pseudotime_table$Pseudotime)
] <- NA_real_

# ------------------------------------------------------------
# Build a fibroblast-only graph for plotting
#
# No UMAP, principal graph, root selection or pseudotime is recalculated.
# We only remove graph branches supported exclusively by mesothelial cells.
# ------------------------------------------------------------

principal_graph_all <- principal_graph(cds_all)[["UMAP"]]
principal_graph_aux_all <- principal_graph_aux(cds_all)[["UMAP"]]

graph_vertex_names <- igraph::V(principal_graph_all)$name

closest_vertex_raw <-
  principal_graph_aux_all$pr_graph_cell_proj_closest_vertex

if (is.matrix(closest_vertex_raw) ||
    is.data.frame(closest_vertex_raw)) {

  closest_vertex <- closest_vertex_raw[, 1]

  if (!is.null(rownames(closest_vertex_raw))) {
    names(closest_vertex) <- rownames(closest_vertex_raw)
  }

} else {
  closest_vertex <- closest_vertex_raw
}

if (is.null(names(closest_vertex))) {
  stop(
    "pr_graph_cell_proj_closest_vertex no contiene nombres de células."
  )
}

closest_vertex_names <- names(closest_vertex)
closest_vertex_character <- as.character(closest_vertex)
numeric_vertex_indices <- suppressWarnings(
  as.integer(closest_vertex_character)
)

# Depending on the Monocle3 version, vertices may be stored as Y_* names
# or as numeric indices.
if (
  all(
    is.na(closest_vertex_character) |
      !is.na(numeric_vertex_indices)
  )
) {

  valid_index <- !is.na(numeric_vertex_indices) &
    numeric_vertex_indices >= 1 &
    numeric_vertex_indices <= length(graph_vertex_names)

  converted_vertices <- rep(
    NA_character_,
    length(numeric_vertex_indices)
  )

  converted_vertices[valid_index] <-
    graph_vertex_names[numeric_vertex_indices[valid_index]]

  closest_vertex <- converted_vertices

} else {
  closest_vertex <- closest_vertex_character
}

names(closest_vertex) <- closest_vertex_names

mapped_fibro_cells <- intersect(
  fibro_cells,
  names(closest_vertex)
)

fibro_principal_vertices <- unique(
  closest_vertex[mapped_fibro_cells]
)

fibro_principal_vertices <- fibro_principal_vertices[
  !is.na(fibro_principal_vertices) &
    fibro_principal_vertices %in% graph_vertex_names
]

if (length(fibro_principal_vertices) < 2) {
  stop(
    "No se identificaron suficientes principal points asociados ",
    "a fibroblastos para dibujar la trayectoria."
  )
}

# Keep the minimal set of paths connecting fibroblast-associated nodes
# within each connected component. This removes the mesothelial-only arm.
graph_components <- igraph::components(
  principal_graph_all
)$membership
names(graph_components) <- graph_vertex_names

fibro_component_ids <- unique(
  graph_components[fibro_principal_vertices]
)

fibro_graph_vertices <- character(0)

for (current_component in fibro_component_ids) {

  current_targets <- fibro_principal_vertices[
    graph_components[fibro_principal_vertices] ==
      current_component
  ]

  current_targets <- unique(current_targets)

  if (length(current_targets) == 1) {
    fibro_graph_vertices <- union(
      fibro_graph_vertices,
      current_targets
    )
  } else {
    reference_vertex <- current_targets[1]

    fibro_graph_vertices <- union(
      fibro_graph_vertices,
      reference_vertex
    )

    for (target_vertex in current_targets[-1]) {
      current_path <- igraph::shortest_paths(
        principal_graph_all,
        from = reference_vertex,
        to = target_vertex,
        output = "vpath"
      )$vpath[[1]]

      fibro_graph_vertices <- union(
        fibro_graph_vertices,
        igraph::as_ids(current_path)
      )
    }
  }
}

principal_graph_fibro_plot <- igraph::induced_subgraph(
  principal_graph_all,
  vids = unique(fibro_graph_vertices)
)

# ------------------------------------------------------------
# Principal-point coordinates and graph segments
# ------------------------------------------------------------

dp_mst <- principal_graph_aux_all$dp_mst

if (
  nrow(dp_mst) == 2 &&
    ncol(dp_mst) == length(graph_vertex_names)
) {
  principal_point_coordinates <- t(dp_mst)

  if (!is.null(colnames(dp_mst))) {
    rownames(principal_point_coordinates) <- colnames(dp_mst)
  }

} else if (
  ncol(dp_mst) == 2 &&
    nrow(dp_mst) == length(graph_vertex_names)
) {
  principal_point_coordinates <- dp_mst

} else {
  stop(
    "No se reconoce la estructura de coordenadas de dp_mst."
  )
}

if (is.null(rownames(principal_point_coordinates))) {
  rownames(principal_point_coordinates) <- graph_vertex_names
}

colnames(principal_point_coordinates)[1:2] <- c(
  "UMAP1",
  "UMAP2"
)

missing_graph_coordinates <- setdiff(
  igraph::V(principal_graph_fibro_plot)$name,
  rownames(principal_point_coordinates)
)

if (length(missing_graph_coordinates) > 0) {
  stop(
    "Faltan coordenadas para estos principal points: ",
    paste(missing_graph_coordinates, collapse = ", ")
  )
}

fibro_graph_edges <- igraph::as_data_frame(
  principal_graph_fibro_plot,
  what = "edges"
)

fibro_graph_segments <- data.frame(
  x = principal_point_coordinates[
    fibro_graph_edges$from,
    "UMAP1"
  ],
  y = principal_point_coordinates[
    fibro_graph_edges$from,
    "UMAP2"
  ],
  xend = principal_point_coordinates[
    fibro_graph_edges$to,
    "UMAP1"
  ],
  yend = principal_point_coordinates[
    fibro_graph_edges$to,
    "UMAP2"
  ],
  stringsAsFactors = FALSE
)

fibro_graph_segments <- fibro_graph_segments[
  complete.cases(fibro_graph_segments),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# Common coordinate limits for both panels
# ------------------------------------------------------------

common_xlim <- range(
  c(
    fibro_pseudotime_table$UMAP1,
    fibro_graph_segments$x,
    fibro_graph_segments$xend
  ),
  finite = TRUE
)

common_ylim <- range(
  c(
    fibro_pseudotime_table$UMAP2,
    fibro_graph_segments$y,
    fibro_graph_segments$yend
  ),
  finite = TRUE
)

add_plot_padding <- function(x, fraction = 0.035) {
  span <- diff(x)

  if (!is.finite(span) || span == 0) {
    span <- 1
  }

  x + c(-1, 1) * span * fraction
}

common_xlim <- add_plot_padding(common_xlim)
common_ylim <- add_plot_padding(common_ylim)

# ------------------------------------------------------------
# Left panel: fibroblast states
# ------------------------------------------------------------

Figure2E_clusters <- ggplot(
  fibro_pseudotime_table,
  aes(
    x = UMAP1,
    y = UMAP2,
    color = Clustering
  )
) +
  geom_point(
    size = 0.8,
    alpha = 1
  ) +
  scale_color_manual(
    values = fibro_colors,
    breaks = fibro_states,
    drop = FALSE,
    name = NULL
  ) +
  guides(
    color = guide_legend(
      direction = "vertical",
      ncol = 1,
      byrow = TRUE,
      override.aes = list(
        size = 2.2,
        alpha = 1
      )
    )
  ) +
  labs(
    title = "Fibroblast states",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  coord_fixed(
    xlim = common_xlim,
    ylim = common_ylim,
    ratio = 1,
    clip = "off"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 4
    ),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

# ------------------------------------------------------------
# Right panel: fibroblast pseudotime
# ------------------------------------------------------------

Figure2E_pseudotime <- ggplot(
  fibro_pseudotime_table,
  aes(
    x = UMAP1,
    y = UMAP2
  )
) +
  geom_segment(
    data = fibro_graph_segments,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.6,
    lineend = "round"
  ) +
  geom_point(
    aes(color = Pseudotime),
    size = 0.8,
    alpha = 1
  ) +
  scale_color_viridis_c(
    option = "plasma",
    name = "Pseudotime",
    limits = c(0, 20),
    breaks = seq(0, 20, by = 5),
    na.value = "grey80"
  ) +
  guides(
    color = guide_colorbar(
      direction = "vertical",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(3, "mm"),
      barheight = unit(38, "mm")
    )
  ) +
  labs(
    title = "Pseudotime",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  coord_fixed(
    xlim = common_xlim,
    ylim = common_ylim,
    ratio = 1,
    clip = "off"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "right",
    legend.justification = "center",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 4
    ),
    plot.margin = margin(5.5, 5.5, 5.5, 7)
  )

# ------------------------------------------------------------
# Combine and save
# ------------------------------------------------------------

Figure2E <- Figure2E_clusters + Figure2E_pseudotime +
  plot_layout(
    ncol = 2,
    widths = figure2e_panel_widths,
    guides = "keep"
  )

ggsave(
  filename = file.path(
    figure_dir,
    "Figure2E_fibroblast_states_pseudotime.pdf"
  ),
  plot = Figure2E,
  width = figure2e_width,
  height = figure2e_height,
  units = "in",
  useDingbats = FALSE,
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure2E_fibroblast_states_pseudotime.png"
  ),
  plot = Figure2E,
  width = figure2e_width,
  height = figure2e_height,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

Figure2E_table <- fibro_pseudotime_table[
  ,
  c(
    "Cell",
    "UMAP1",
    "UMAP2",
    "Clustering",
    "Pseudotime_raw",
    "Pseudotime"
  )
]

write.table(
  Figure2E_table,
  file.path(
    table_dir,
    "Figure2E_fibroblast_states_pseudotime.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nFigure 2E finalizada.\n",
  "PDF: ",
  file.path(
    figure_dir,
    "Figure2E_fibroblast_states_pseudotime.pdf"
  ),
  "\nPNG: ",
  file.path(
    figure_dir,
    "Figure2E_fibroblast_states_pseudotime.png"
  ),
  "\n",
  sep = ""
)

# ============================================================
# Figure 2F: dynamic genes along fibroblast pseudotime
# Mesothelial cells are absent only after the complete trajectory
# and pseudotime have already been calculated in cds_all.
# ============================================================

cat("\nEjecutando graph_test en fibroblastos...\n")

modulated_genes <- graph_test(
  cds_fibro,
  neighbor_graph = "principal_graph",
  cores = n_cores
)

graph_test_table <- as.data.frame(modulated_genes)
graph_test_table$Gene <- rownames(graph_test_table)

write.table(
  graph_test_table,
  file.path(table_dir, "Figure2F_graph_test_all.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sig_genes <- graph_test_table[
  graph_test_table$status == "OK" &
    graph_test_table$q_value < 0.05 &
    is.finite(graph_test_table$morans_I),
  ,
  drop = FALSE
]

sig_genes <- sig_genes[
  order(sig_genes$morans_I, decreasing = TRUE),
  ,
  drop = FALSE
]

if (nrow(sig_genes) == 0) {
  stop("No se encontraron genes significativos en graph_test().")
}

# ============================================================
# Select 60 genes for the pseudotime heatmap
# Highest Moran's I genes plus mandatory manually selected genes
# ============================================================

manual_genes <- c(
  "Axl",
  "Pdgfra",
  "Pdgfrb",
  "Dcn",
  "Pdpn",
  "Notch3",
  "Lipg",
  "Per1",
  "Col24a1",
  "Macf1"
)

excluded_genes <- c(
  "Mapkapk3",
  "Dpep1",
  "Mgp",
  "Cfh",
  "Spon2",
  "Tgfbi",
  "Txnip"
)

manual_genes <- unique(
  setdiff(
    manual_genes,
    excluded_genes
  )
)

missing_manual_genes <- setdiff(
  manual_genes,
  rownames(rowData(cds_fibro))
)

if (length(missing_manual_genes) > 0) {
  stop(
    "Estos genes manuales no están presentes en cds_fibro: ",
    paste(missing_manual_genes, collapse = ", ")
  )
}

final_gene_number <- 60L
n_automatic_genes <- final_gene_number -
  length(manual_genes)

if (n_automatic_genes < 0) {
  stop(
    "La lista manual contiene más de ",
    final_gene_number,
    " genes."
  )
}

automatic_candidates <- sig_genes[
  !sig_genes$Gene %in% c(
    manual_genes,
    excluded_genes
  ),
  ,
  drop = FALSE
]

if (nrow(automatic_candidates) < n_automatic_genes) {
  stop(
    "No hay suficientes genes automáticos después de las exclusiones. ",
    "Disponibles: ",
    nrow(automatic_candidates),
    "; necesarios: ",
    n_automatic_genes,
    "."
  )
}

automatic_genes <- head(
  automatic_candidates$Gene,
  n_automatic_genes
)

valid_genes <- unique(
  c(
    automatic_genes,
    manual_genes
  )
)

valid_genes <- intersect(
  valid_genes,
  rownames(rowData(cds_fibro))
)

valid_genes <- setdiff(
  valid_genes,
  excluded_genes
)

if (length(valid_genes) != final_gene_number) {
  stop(
    "La selección final contiene ",
    length(valid_genes),
    " genes en lugar de ",
    final_gene_number,
    "."
  )
}

if (!all(manual_genes %in% valid_genes)) {
  stop(
    "No se conservaron todos los genes manuales en la selección final."
  )
}

top_genes <- graph_test_table[
  match(
    valid_genes,
    graph_test_table$Gene
  ),
  ,
  drop = FALSE
]

top_genes$Gene <- valid_genes

top_genes$Selection <- ifelse(
  top_genes$Gene %in% manual_genes,
  "Manual",
  "Moran_I"
)

rownames(top_genes) <- top_genes$Gene

cat(
  "\nGenes seleccionados para Figure 2F:\n",
  "  Automáticos por Moran's I: ",
  length(automatic_genes),
  "\n  Incluidos manualmente: ",
  length(manual_genes),
  "\n  Total: ",
  length(valid_genes),
  "\n  Genes manuales: ",
  paste(manual_genes, collapse = ", "),
  "\n  Genes excluidos: ",
  paste(excluded_genes, collapse = ", "),
  "\n",
  sep = ""
)

finite_cells <- names(pseudotime_fibro)[is.finite(pseudotime_fibro)]
ordered_cells <- finite_cells[
  order(pseudotime_fibro[finite_cells])
]
ordered_pseudotime <- pseudotime_fibro[ordered_cells]

pt_matrix <- exprs(cds_fibro)[
  valid_genes,
  ordered_cells,
  drop = FALSE
]

cat("Aplicando smooth.spline(df = 3)...\n")

pt_matrix <- t(
  apply(
    pt_matrix,
    1,
    function(x) {
      x <- as.numeric(x)
      if (length(unique(x)) < 2 || sd(x) == 0) {
        return(rep(0, length(x)))
      }
      smooth.spline(x = seq_along(x), y = x, df = 3)$y
    }
  )
)

rownames(pt_matrix) <- valid_genes
colnames(pt_matrix) <- ordered_cells

pt_matrix <- t(
  apply(
    pt_matrix,
    1,
    function(x) {
      current_sd <- sd(x)
      if (!is.finite(current_sd) || current_sd == 0) {
        return(rep(0, length(x)))
      }
      (x - mean(x)) / current_sd
    }
  )
)

pt_matrix <- pt_matrix[complete.cases(pt_matrix), , drop = FALSE]

if (nrow(pt_matrix) < 3) {
  stop("Hay menos de tres genes después del suavizado y escalado.")
}

expression_color <- colorRamp2(
  seq(-2, 2, length.out = 11),
  rev(brewer.pal(11, "Spectral"))
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

heatmap_metadata <- as.data.frame(colData(cds_fibro))[
  ordered_cells,
  ,
  drop = FALSE
]

cell_clusters <- factor(
  as.character(heatmap_metadata$Clustering),
  levels = fibro_states
)

cell_group <- factor(
  as.character(heatmap_metadata$group),
  levels = c("Naive", "Tumor")
)

cell_sample <- factor(
  as.character(heatmap_metadata$orig.ident),
  levels = c("WT", "CD169DTR", "CCR2KO")
)

pseudotime_color <- colorRamp2(
  c(0, 20),
  c("#F7FCF5", "#238B45")
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
  annotation_legend_param = list(
    Pseudotime = list(
      at = seq(0, 20, by = 5)
    )
  ),
  simple_anno_size = unit(4, "mm")
)
# ============================================================
# Hierarchical clustering of genes
# Ward.D2 dendrogram with manual branch orientation
# ============================================================

row_hc <- hclust(
  dist(pt_matrix),
  method = "ward.D2"
)

row_dend <- as.dendrogram(
  row_hc
)

# ------------------------------------------------------------
# Function to exchange the two immediate branches
# ------------------------------------------------------------

swap_dendrogram_children <- function(dend) {

  if (!inherits(dend, "dendrogram")) {
    stop("El objeto proporcionado no es un dendrograma.")
  }

  if (
    is.leaf(dend) ||
      length(dend) != 2L
  ) {
    stop(
      "La rama seleccionada no contiene exactamente dos subramas."
    )
  }

  first_branch <- dend[[1]]
  second_branch <- dend[[2]]

  dend[[1]] <- second_branch
  dend[[2]] <- first_branch

  dend
}

# ------------------------------------------------------------
# Exchange only the two main clusters
# ------------------------------------------------------------

row_dend <- swap_dendrogram_children(
  row_dend
)

# IMPORTANTE:
# No intercambiamos row_dend[[1]], porque corresponde al
# cluster superior, cuya orientación ya es correcta.
#
# Tampoco intercambiamos row_dend[[2]] aquí. Al eliminar el
# intercambio que tenía el script anterior, sus dos subclusters
# quedan invertidos respecto al heatmap actual.

# ------------------------------------------------------------
# Final visual gene order
# ------------------------------------------------------------

heatmap_gene_order <- labels(
  row_dend
)

final_gene_order <- heatmap_gene_order

cat(
  "\nGenes incluidos en el clustering jerárquico: ",
  length(final_gene_order),
  "\n",
  sep = ""
)

Figure2F <- Heatmap(
  pt_matrix,
  name = "Scaled Expression",
  col = expression_color,

  cluster_rows = row_dend,
  cluster_columns = FALSE,
  show_row_dend = TRUE,

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

pdf(
  file = file.path(
    figure_dir,
    "Figure2F_pseudotime_heatmap.pdf"
  ),
  width = 14,
  height = 8.5,
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
  width = 4000,
  height = 2400,
  res = 250
)

draw(
  Figure2F,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

selected_gene_table <- top_genes[
  match(heatmap_gene_order, top_genes$Gene),
  ,
  drop = FALSE
]

selected_gene_table$Heatmap_order <- match(
  selected_gene_table$Gene,
  final_gene_order
)

selected_gene_table <- selected_gene_table[
  order(selected_gene_table$Heatmap_order),
  ,
  drop = FALSE
]

write.table(
  selected_gene_table,
  file.path(table_dir, "Figure2F_selected_genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

scaled_expression_table <- data.frame(
  Gene = final_gene_order,
  pt_matrix[final_gene_order, , drop = FALSE],
  check.names = FALSE
)

write.table(
  scaled_expression_table,
  file.path(table_dir, "Figure2F_scaled_expression.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nFigure 2F finalizada con clustering jerárquico Ward.D2.\n",
  "PDF: ", file.path(figure_dir, "Figure2F_pseudotime_heatmap.pdf"),
  "\nPNG: ", file.path(figure_dir, "Figure2F_pseudotime_heatmap.png"),
  "\n",
  sep = ""
)

# ============================================================
# Figure 2G
# Cytokines, collagens and ECM heatmaps
# ============================================================

figure2g_states <- trajectory_states

figure2g_cells <- rownames(meta)[
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
fibro_meso <- JoinLayers(fibro_meso, assay = "RNA")

cat("\nCélulas para los heatmaps de Figure 2G:\n")
print(table(fibro_meso$Clustering, useNA = "ifany"))

cytokine_genes <- c(
  "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl7", "Ccl8", "Ccl11",
  "Ccl12", "Ccl17", "Ccl19", "Ccl20", "Ccl21", "Ccl22",
  "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl5", "Cxcl9", "Cxcl10",
  "Cxcl11", "Cxcl12", "Cxcl13", "Cxcl14", "Cxcl16",
  "Il1a", "Il1b", "Il6", "Il7", "Il11", "Il15", "Il18",
  "Il33", "Il34", "Tnf", "Tnfsf10", "Tnfsf11", "Tnfsf12",
  "Tnfsf13b", "Tnfsf14", "Lta", "Tgfb1", "Tgfb2", "Tgfb3",
  "Lif", "Osm", "Csf1", "Csf2", "Csf3", "Kitl", "Flt3l",
  "Areg", "Ereg", "Gdf15", "Vegfa", "Ptgs2"
)

collagen_genes <- c(
  "Col1a1", "Col1a2", "Col2a1", "Col3a1",
  "Col4a1", "Col4a2", "Col4a3", "Col4a4", "Col4a5", "Col4a6",
  "Col5a1", "Col5a2", "Col5a3",
  "Col6a1", "Col6a2", "Col6a3", "Col6a4", "Col6a5", "Col6a6",
  "Col7a1", "Col8a1", "Col8a2", "Col9a1", "Col9a2", "Col9a3",
  "Col10a1", "Col11a1", "Col11a2", "Col12a1", "Col13a1",
  "Col14a1", "Col15a1", "Col16a1", "Col17a1", "Col18a1",
  "Col19a1", "Col20a1", "Col21a1", "Col22a1", "Col23a1",
  "Col24a1", "Col25a1", "Col26a1", "Col27a1", "Col28a1"
)

ecm_genes <- c(
  "Fn1", "Fbn1", "Fbn2", "Postn", "Tnc", "Thbs1", "Thbs2",
  "Sparc", "Sparcl1", "Vtn", "Nid1", "Nid2", "Emilin1", "Emilin2",
  "Mfap2", "Mfap4", "Mfap5", "Eln", "Efemp1", "Efemp2",
  "Fbln1", "Fbln2", "Fbln5", "Fbln7", "Dcn", "Lum", "Bgn",
  "Fmod", "Aspn", "Ogn", "Prelp", "Hspg2", "Sdc1", "Sdc2",
  "Sdc3", "Sdc4", "Gpc1", "Gpc3", "Gpc4", "Gpc6", "Ccn1",
  "Ccn2", "Ccn3", "Ccn4", "Ccn5", "Lox", "Loxl1", "Loxl2",
  "Loxl3", "Loxl4", "Serpine1", "Serpinh1", "Pcolce", "Pcolce2",
  "Plod1", "Plod2", "Plod3", "P4ha1", "P4ha2", "P4ha3",
  "Timp1", "Timp2", "Timp3", "Timp4", "Acta2", "Tagln", "Myl9",
  "Myh11", "Cnn1"
)

make_figure2g_heatmap <- function(object, genes, title, filename) {
  genes <- intersect(genes, rownames(object))
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

  state_keys <- make.names(figure2g_states)
  column_keys <- make.names(colnames(avg_mat))
  column_index <- match(state_keys, column_keys)
  keep_states <- !is.na(column_index)

  avg_mat <- avg_mat[, column_index[keep_states], drop = FALSE]
  colnames(avg_mat) <- figure2g_states[keep_states]
  avg_mat <- avg_mat[rowSums(avg_mat > 0, na.rm = TRUE) > 0, , drop = FALSE]

  row_mean <- rowMeans(avg_mat, na.rm = TRUE)
  row_sd <- apply(avg_mat, 1, sd, na.rm = TRUE)
  z_mat <- sweep(avg_mat, 1, row_mean, "-")
  z_mat <- sweep(z_mat, 1, row_sd, "/")
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
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    breaks = seq(-2.5, 2.5, length.out = 101),
    main = title,
    silent = TRUE
  )

  plot_width <- max(14, ncol(heatmap_matrix) * 0.32)

  pdf(
    file.path(figure_dir, paste0(filename, ".pdf")),
    width = plot_width,
    height = 4.5,
    useDingbats = FALSE
  )
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()

  png(
    file.path(figure_dir, paste0(filename, ".png")),
    width = plot_width,
    height = 4.5,
    units = "in",
    res = 300
  )
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()

  write.table(
    data.frame(Gene = rownames(z_mat), z_mat, check.names = FALSE),
    file.path(table_dir, paste0(filename, "_zscore.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat(title, ": ", nrow(z_mat), " genes\n", sep = "")
}

make_figure2g_heatmap(
  fibro_meso,
  cytokine_genes,
  "Cytokines / chemokines",
  "Figure2G_cytokines_chemokines"
)

make_figure2g_heatmap(
  fibro_meso,
  collagen_genes,
  "Collagens",
  "Figure2G_collagens"
)

make_figure2g_heatmap(
  fibro_meso,
  ecm_genes,
  "ECM / matrix-associated genes",
  "Figure2G_ECM"
)

capture.output(
  sessionInfo(),
  file = file.path(table_dir, "Figure2_sessionInfo.txt")
)

cat("\nFigure 2 finalizada.\n")