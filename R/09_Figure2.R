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
  "Figure2_fibroblast_before_ordering.rds"
)

manual_root_nodes_file <- file.path(
  table_dir,
  "Figure2E_manual_root_pr_nodes_fibro_only.txt"
)

root_selection_summary_file <- file.path(
  table_dir,
  "Figure2E_root_selection_used_fibro_only.tsv"
)


# ============================================================
# Load annotated fibroblast object
# ============================================================

input_object_file <- file.path(
  objects_dir,
  "fibroblasts.Fibros.annotated.rds"
)

if (!file.exists(input_object_file)) {
  stop(
    "No se encuentra el objeto de fibroblastos anotados:\n",
    input_object_file
  )
}

fibroblasts <- readRDS(
  input_object_file
)

DefaultAssay(fibroblasts) <- "RNA"

required_metadata <- c(
  "fibro_annotation",
  "group",
  "orig.ident"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(fibroblasts[[]])
)

if (length(missing_metadata) > 0) {
  stop(
    "Faltan columnas de metadatos en el objeto de entrada: ",
    paste(missing_metadata, collapse = ", ")
  )
}

fibro_states <- c(
  "Fibroblasts Universal",
  "Fibroblasts Activated 1",
  "Fibroblasts Activated 2"
)

# Estados permitidos en el objeto de origen.
# Mesothelial puede estar presente en el RDS, pero se elimina
# antes de construir el objeto de Monocle.
source_states <- c(
  fibro_states,
  "Mesothelial"
)

trajectory_states <- fibro_states

observed_states <- sort(
  unique(
    na.omit(
      as.character(
        fibroblasts$fibro_annotation
      )
    )
  )
)

unexpected_states <- setdiff(
  observed_states,
  source_states
)

if (length(unexpected_states) > 0) {
  stop(
    "Se han encontrado anotaciones inesperadas en ",
    "'fibro_annotation':\n",
    paste(unexpected_states, collapse = "\n")
  )
}

if (anyNA(fibroblasts$fibro_annotation)) {
  stop(
    "El objeto contiene células con fibro_annotation = NA."
  )
}

fibroblasts$Clustering <- factor(
  as.character(
    fibroblasts$fibro_annotation
  ),
  levels = source_states
)

# Metadatos que permiten comprobar el origen y el tipo de análisis
# cuando se reutiliza un CDS ya calculado.
fibroblasts$Figure2_source_object <- basename(
  input_object_file
)

fibroblasts$Figure2_analysis <- "fibroblast_only_seurat_umap"
fibroblasts$Figure2_embedding <- "Seurat_umap"

trajectory_colors <- c(
  "Fibroblasts Universal" = "#FBCDC7",
  "Fibroblasts Activated 1" = "#F37B79",
  "Fibroblasts Activated 2" = "#6D1A1A"
)

fibro_colors <- trajectory_colors[
  fibro_states
]

cat("\nObjeto de entrada de Figure 2:\n")
cat(input_object_file, "\n")

cat("\nAnotaciones presentes en el objeto de origen:\n")
print(
  table(
    fibroblasts$Clustering,
    useNA = "ifany"
  )
)

cat(
  "\nCélulas totales del objeto de origen: ",
  ncol(fibroblasts),
  "\n",
  sep = ""
)

# ============================================================
# Build trajectory input with fibroblasts only
# ============================================================

meta <- fibroblasts[[]]

trajectory_cells <- rownames(meta)[
  as.character(meta$Clustering) %in% fibro_states
]

trajectory_subset <- subset(
  fibroblasts,
  cells = trajectory_cells
)

trajectory_subset$Clustering <- factor(
  as.character(
    trajectory_subset$Clustering
  ),
  levels = fibro_states
)

Idents(trajectory_subset) <- "Clustering"
DefaultAssay(trajectory_subset) <- "RNA"

cat("\nCélulas utilizadas para construir la trayectoria:\n")
print(
  table(
    trajectory_subset$Clustering,
    useNA = "ifany"
  )
)

missing_fibro_states <- setdiff(
  fibro_states,
  unique(
    as.character(
      trajectory_subset$Clustering
    )
  )
)

if (length(missing_fibro_states) > 0) {
  stop(
    "Faltan estados de fibroblastos: ",
    paste(missing_fibro_states, collapse = ", ")
  )
}

if (
  any(
    as.character(
      trajectory_subset$Clustering
    ) == "Mesothelial"
  )
) {
  stop(
    "Se han incorporado células Mesothelial por error."
  )
}

cat(
  "Número total de fibroblastos incorporados: ",
  ncol(trajectory_subset),
  "\n",
  sep = ""
)

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

# ============================================================
# Recover the original Seurat UMAP
# No PCA, UMAP or clustering is recalculated in Monocle
# ============================================================

if (!"umap" %in% Reductions(trajectory_subset)) {
  stop(
    "El objeto de fibroblastos no contiene la reducción 'umap'."
  )
}

seurat_umap <- Embeddings(
  trajectory_subset,
  reduction = "umap"
)

if (ncol(seurat_umap) < 2) {
  stop("La reducción UMAP de Seurat tiene menos de dos dimensiones.")
}

missing_umap_cells <- setdiff(
  colnames(counts_matrix),
  rownames(seurat_umap)
)

if (length(missing_umap_cells) > 0) {
  stop(
    "Faltan células en el UMAP original de Seurat: ",
    paste(head(missing_umap_cells, 10), collapse = ", ")
  )
}

seurat_umap <- seurat_umap[
  colnames(counts_matrix),
  1:2,
  drop = FALSE
]

colnames(seurat_umap) <- c(
  "UMAP1",
  "UMAP2"
)

expected_analysis <- "fibroblast_only_seurat_umap"
expected_embedding <- "Seurat_umap"

if (file.exists(preorder_cds_file)) {
  cat(
    "\nReutilizando el CDS previo a order_cells:\n",
    preorder_cds_file,
    "\n",
    sep = ""
  )

  cds_all <- readRDS(
    preorder_cds_file
  )

  required_cds_metadata <- c(
    "Figure2_source_object",
    "Figure2_analysis",
    "Figure2_embedding"
  )

  missing_cds_metadata <- setdiff(
    required_cds_metadata,
    colnames(colData(cds_all))
  )

  if (length(missing_cds_metadata) > 0) {
    stop(
      "El CDS existente no contiene estos metadatos de control: ",
      paste(missing_cds_metadata, collapse = ", "),
      "\nElimina este archivo y vuelve a ejecutar el script:\n",
      preorder_cds_file
    )
  }

  cds_source <- unique(
    as.character(
      colData(cds_all)$Figure2_source_object
    )
  )

  cds_source <- cds_source[
    !is.na(cds_source) & nzchar(cds_source)
  ]

  expected_source <- basename(
    input_object_file
  )

  if (
    length(cds_source) != 1L ||
      cds_source != expected_source
  ) {
    stop(
      "El CDS existente fue generado desde otro objeto.\n",
      "Objeto esperado: ",
      expected_source,
      "\nObjeto registrado en el CDS: ",
      paste(cds_source, collapse = ", "),
      "\n\nElimina este archivo y vuelve a ejecutar:\n",
      preorder_cds_file
    )
  }

  cds_analysis <- unique(
    as.character(
      colData(cds_all)$Figure2_analysis
    )
  )

  cds_analysis <- cds_analysis[
    !is.na(cds_analysis) & nzchar(cds_analysis)
  ]

  if (
    length(cds_analysis) != 1L ||
      cds_analysis != expected_analysis
  ) {
    stop(
      "El CDS existente no corresponde al análisis con UMAP de Seurat.\n",
      "Tipo esperado: ", expected_analysis,
      "\nTipo registrado: ",
      paste(cds_analysis, collapse = ", "),
      "\n\nElimina este archivo y vuelve a ejecutar:\n",
      preorder_cds_file
    )
  }

  cds_embedding <- unique(
    as.character(
      colData(cds_all)$Figure2_embedding
    )
  )

  cds_embedding <- cds_embedding[
    !is.na(cds_embedding) & nzchar(cds_embedding)
  ]

  if (
    length(cds_embedding) != 1L ||
      cds_embedding != expected_embedding
  ) {
    stop(
      "El CDS existente no contiene el UMAP original de Seurat.\n",
      "Embedding esperado: ", expected_embedding,
      "\nEmbedding registrado: ",
      paste(cds_embedding, collapse = ", "),
      "\n\nElimina este archivo y vuelve a ejecutar:\n",
      preorder_cds_file
    )
  }

  if (!"UMAP" %in% SingleCellExperiment::reducedDimNames(cds_all)) {
    stop(
      "El CDS existente no contiene la reducción UMAP.\n",
      "Elimina este archivo y vuelve a ejecutar:\n",
      preorder_cds_file
    )
  }

  stored_umap <- SingleCellExperiment::reducedDim(
    cds_all,
    "UMAP"
  )

  missing_stored_cells <- setdiff(
    rownames(seurat_umap),
    rownames(stored_umap)
  )

  if (length(missing_stored_cells) > 0) {
    stop(
      "El CDS existente no contiene todas las células del UMAP original.\n",
      "Elimina este archivo y vuelve a ejecutar:\n",
      preorder_cds_file
    )
  }

  stored_umap <- stored_umap[
    rownames(seurat_umap),
    1:2,
    drop = FALSE
  ]

  umap_max_difference <- max(
    abs(
      unname(stored_umap) -
        unname(seurat_umap)
    )
  )

  if (!isTRUE(all.equal(umap_max_difference, 0, tolerance = 0))) {
    stop(
      "El UMAP del CDS no es idéntico al UMAP original de Seurat.\n",
      "Diferencia máxima absoluta: ",
      umap_max_difference,
      "\nElimina este archivo y vuelve a ejecutar:\n",
      preorder_cds_file
    )
  }

  cat(
    "UMAP verificado: diferencia máxima absoluta = 0.\n"
  )

} else {
  cds_all <- new_cell_data_set(
    expression_data = counts_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )

  cds_all <- estimate_size_factors(
    cds_all
  )

  # Import the exact Seurat coordinates.
  SingleCellExperiment::reducedDim(
    cds_all,
    "UMAP"
  ) <- seurat_umap[
    colnames(cds_all),
    ,
    drop = FALSE
  ]

  # Register the existing fibroblast annotations in Monocle's
  # internal structure. This does not calculate new clusters.
  seurat_fibro_states <- factor(
    as.character(
      colData(cds_all)$Clustering
    ),
    levels = fibro_states
  )

  names(seurat_fibro_states) <- colnames(
    cds_all
  )

  single_partition <- factor(
    rep(
      "1",
      ncol(cds_all)
    ),
    levels = "1"
  )

  names(single_partition) <- colnames(
    cds_all
  )

  cds_all@clusters[["UMAP"]] <- list(
    clusters = seurat_fibro_states,
    partitions = single_partition
  )

  imported_umap <- SingleCellExperiment::reducedDim(
    cds_all,
    "UMAP"
  )

  imported_umap <- imported_umap[
    rownames(seurat_umap),
    1:2,
    drop = FALSE
  ]

  umap_max_difference <- max(
    abs(
      unname(imported_umap) -
        unname(seurat_umap)
    )
  )

  if (!isTRUE(all.equal(umap_max_difference, 0, tolerance = 0))) {
    stop(
      "La importación del UMAP de Seurat ha modificado las coordenadas.\n",
      "Diferencia máxima absoluta: ",
      umap_max_difference
    )
  }

  cat(
    "\nUMAP original de Seurat incorporado sin recalcularlo.\n",
    "Diferencia máxima absoluta respecto a Seurat: ",
    umap_max_difference,
    "\n",
    sep = ""
  )

  cat(
    "Anotaciones originales registradas desde fibro_annotation.\n",
    "No se ha calculado un clustering nuevo.\n"
  )

  set.seed(1234)
  cat(
    "Aprendiendo principal graph sobre el UMAP original de Seurat...\n"
  )

  cds_all <- learn_graph(
    cds_all,
    use_partition = FALSE,
    learn_graph_control = list(
      minimal_branch_len = 40,
      prune_graph = FALSE
    )
  )

  saveRDS(
    cds_all,
    preorder_cds_file
  )

  cat(
    "CDS previo a order_cells guardado en: ",
    preorder_cds_file,
    "\n",
    sep = ""
  )
}

cat("\nFibroblastos presentes en el CDS usado para seleccionar el root:\n")
print(
  table(
    colData(cds_all)$Clustering,
    useNA = "ifany"
  )
)

if (
  any(
    as.character(
      colData(cds_all)$Clustering
    ) == "Mesothelial"
  )
) {
  stop(
    "El CDS contiene células Mesothelial. ",
    "Borra el RDS previo y vuelve a ejecutar el script."
  )
}

cds_states <- unique(
  as.character(
    colData(cds_all)$Clustering
  )
)

cds_states <- cds_states[
  !is.na(cds_states)
]

unexpected_cds_states <- setdiff(
  cds_states,
  fibro_states
)

if (length(unexpected_cds_states) > 0) {
  stop(
    "El CDS contiene poblaciones distintas de los tres ",
    "estados de fibroblastos: ",
    paste(unexpected_cds_states, collapse = ", ")
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
    "\nEl grafo con fibroblastos ya está preparado.\n",
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

# ============================================================
# Fibroblast-only CDS and pseudotime scale
# ============================================================

# The CDS already contains only fibroblasts. No dimensional reduction,
# graph learning or cell ordering is repeated.
fibro_cells <- colnames(cds_all)
cds_fibro <- cds_all

pseudotime_fibro_raw <- pseudotime_all[
  fibro_cells
]

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

fibro_subset_plot <- trajectory_subset

fibro_subset_plot$Clustering <- factor(
  as.character(
    fibro_subset_plot$Clustering
  ),
  levels = fibro_states
)

Idents(fibro_subset_plot) <- "Clustering"

cat("\nCélulas usadas para calcular el pseudotime:\n")
print(
  table(
    colData(cds_fibro)$Clustering,
    useNA = "ifany"
  )
)

if (
  any(
    as.character(
      colData(cds_fibro)$Clustering
    ) == "Mesothelial"
  )
) {
  stop("Mesothelial está presente en el CDS fibroblast-only.")
}

# ============================================================
# Figure 2E: fibroblast states and fibroblast pseudotime
#
# The UMAP, principal graph, manual root selection and pseudotime
# are calculated using only fibroblast cells.
# ============================================================

# Figure 2E: final dimensions and relative panel widths.
# coord_fixed(ratio = 1) preserves the UMAP geometry.
figure2e_width <- 14
figure2e_height <- 5.2
figure2e_panel_widths <- c(1.05, 1)

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
# Principal graph used for plotting
#
# The graph was learned using only fibroblast cells.
# ------------------------------------------------------------

principal_graph_all <- principal_graph(
  cds_all
)[["UMAP"]]

principal_graph_aux_all <- principal_graph_aux(
  cds_all
)[["UMAP"]]

if (is.null(principal_graph_all)) {
  stop("El CDS no contiene principal graph para UMAP.")
}

principal_graph_fibro_plot <- principal_graph_all

graph_vertex_names <- igraph::V(
  principal_graph_fibro_plot
)$name

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
    ratio = 0.75,
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
    ratio = 0.75,
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
# The trajectory and pseudotime contain fibroblasts only.
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
  "Macf1",
  "Pnoc",
  "Cxcl13",
  "Styk1",
  "Ptn",
  "Tnc"
)

excluded_genes <- c(
  "Mapkapk3",
  "Dpep1",
  "Mgp",
  "Cfh",
  "Spon2",
  "C3",
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
# Order the three gene blocks by the pseudotime of maximum
# expression:
#
#   1. Early pseudotime
#   2. Intermediate pseudotime
#   3. Late pseudotime
#
# The hierarchical clustering is preserved. Only the visual
# orientation of compatible dendrogram branches is changed.
# ------------------------------------------------------------

row_clusters <- cutree(
  row_hc,
  k = 3
)

if (length(unique(row_clusters)) != 3L) {
  stop(
    "No se identificaron exactamente tres bloques de genes."
  )
}

# Pseudotime at which each smoothed gene profile reaches
# its maximum expression.
gene_peak_index <- apply(
  pt_matrix,
  1,
  which.max
)

gene_peak_pseudotime <- ordered_pseudotime[
  gene_peak_index
]

names(gene_peak_pseudotime) <- rownames(
  pt_matrix
)

# Representative peak pseudotime for each hierarchical block.
cluster_peak_pseudotime <- tapply(
  gene_peak_pseudotime[names(row_clusters)],
  row_clusters,
  median,
  na.rm = TRUE
)

if (
  length(cluster_peak_pseudotime) != 3L ||
    any(!is.finite(cluster_peak_pseudotime))
) {
  stop(
    "No se pudo calcular el pseudotime de pico de los tres bloques."
  )
}

# Early -> intermediate -> late.
desired_cluster_order <- as.integer(
  names(
    sort(
      cluster_peak_pseudotime,
      decreasing = FALSE
    )
  )
)

current_gene_order <- labels(
  row_dend
)

# Preserve the internal dendrogram order of genes within
# each block.
desired_gene_order <- unlist(
  lapply(
    desired_cluster_order,
    function(cluster_id) {
      current_gene_order[
        row_clusters[current_gene_order] == cluster_id
      ]
    }
  ),
  use.names = FALSE
)

if (
  length(desired_gene_order) != nrow(pt_matrix) ||
    anyDuplicated(desired_gene_order) ||
    !setequal(desired_gene_order, rownames(pt_matrix))
) {
  stop(
    "El nuevo orden no contiene exactamente todos los genes."
  )
}

# Safe branch rotation: dendrogram attributes are preserved.
row_dend <- dendextend::rotate(
  row_dend,
  order = desired_gene_order
)

# Confirm that the three major blocks were ordered correctly.
major_block_gene_order <- labels(
  row_dend
)

if (!identical(
  as.character(major_block_gene_order),
  as.character(desired_gene_order)
)) {
  stop(
    "No se pudo obtener el orden general ",
    "early-intermediate-late solicitado."
  )
}

# ------------------------------------------------------------
# Local rotations within the hierarchical dendrogram
#
# These rotations preserve:
#   - Ward.D2 clustering
#   - branch heights
#   - gene membership
#   - distances
#
# Only the visual orientation of sibling branches is changed.
# ------------------------------------------------------------

place_sibling_block_after <- function(
  dend,
  block_to_move,
  reference_block,
  description
) {

  dend_labels <- labels(
    dend
  )

  missing_genes <- setdiff(
    c(
      block_to_move,
      reference_block
    ),
    dend_labels
  )

  if (length(missing_genes) > 0) {
    stop(
      "No se encuentran estos genes para la rotación '",
      description,
      "': ",
      paste(missing_genes, collapse = ", ")
    )
  }

  find_sibling_split <- function(node) {

    if (is.leaf(node)) {
      return(NULL)
    }

    left_genes <- labels(
      node[[1]]
    )

    right_genes <- labels(
      node[[2]]
    )

    move_in_left <- all(
      block_to_move %in% left_genes
    )

    move_in_right <- all(
      block_to_move %in% right_genes
    )

    reference_in_left <- all(
      reference_block %in% left_genes
    )

    reference_in_right <- all(
      reference_block %in% right_genes
    )

    # The two requested gene groups are located in
    # opposite immediate branches of this node.
    if (
      move_in_left &&
        reference_in_right
    ) {
      return(
        list(
          node_genes = labels(node),
          block_to_move_genes = left_genes,
          reference_branch_genes = right_genes
        )
      )
    }

    if (
      move_in_right &&
        reference_in_left
    ) {
      return(
        list(
          node_genes = labels(node),
          block_to_move_genes = right_genes,
          reference_branch_genes = left_genes
        )
      )
    }

    target_genes <- unique(
      c(
        block_to_move,
        reference_block
      )
    )

    # Continue down the only branch containing both groups.
    if (all(target_genes %in% left_genes)) {
      return(
        find_sibling_split(
          node[[1]]
        )
      )
    }

    if (all(target_genes %in% right_genes)) {
      return(
        find_sibling_split(
          node[[2]]
        )
      )
    }

    NULL
  }

  sibling_split <- find_sibling_split(
    dend
  )

  if (is.null(sibling_split)) {
    stop(
      "No se encontraron ramas hermanas compatibles para: ",
      description
    )
  }

  current_order <- labels(
    dend
  )

  node_positions <- which(
    current_order %in%
      sibling_split$node_genes
  )

  if (
    length(node_positions) == 0L ||
      !identical(
        current_order[node_positions],
        sibling_split$node_genes
      )
  ) {
    stop(
      "No se pudo localizar de forma continua la rama para: ",
      description
    )
  }

  # Desired local order:
  # reference branch first, moved branch afterwards.
  desired_local_order <- c(
    sibling_split$reference_branch_genes,
    sibling_split$block_to_move_genes
  )

  desired_full_order <- current_order
  desired_full_order[node_positions] <-
    desired_local_order

  # Nothing to do if it is already correctly oriented.
  if (identical(
    current_order,
    desired_full_order
  )) {
    cat(
      "La rotación ya estaba aplicada: ",
      description,
      "\n",
      sep = ""
    )

    return(dend)
  }

  rotated_dend <- dendextend::rotate(
    dend,
    order = desired_full_order
  )

  if (!identical(
    as.character(labels(rotated_dend)),
    as.character(desired_full_order)
  )) {
    stop(
      "No se pudo aplicar la rotación local: ",
      description
    )
  }

  cat(
    "Rotación local aplicada: ",
    description,
    "\n",
    sep = ""
  )

  rotated_dend
}

# ------------------------------------------------------------
# 1. Move the Inmt/Fn1/Col15a1 branch below the
#    neighbouring early-pseudotime branch ending in Fbln2.
# ------------------------------------------------------------

row_dend <- place_sibling_block_after(
  dend = row_dend,
  block_to_move = c(
    "Inmt",
    "Fn1",
    "Col15a1"
  ),
  reference_block = c(
    "Fbln2"
  ),
  description = paste0(
    "rama Inmt/Fn1/Col15a1 después ",
    "de la rama que termina en Fbln2"
  )
)

# ------------------------------------------------------------
# 2. Move the Nampt/Tnc branch below the
#    Fmo2/Ddit4/Per1 branch.
# ------------------------------------------------------------

row_dend <- place_sibling_block_after(
  dend = row_dend,
  block_to_move = c(
    "Nampt",
    "Tnc"
  ),
  reference_block = c(
    "Fmo2",
    "Ddit4",
    "Per1"
  ),
  description = paste0(
    "rama Nampt/Tnc después ",
    "de la rama Fmo2/Ddit4/Per1"
  )
)

# ------------------------------------------------------------
# 3. Move the Notch3/Ptn branch below the
#    Pdgfrb/Macf1 branch.
# ------------------------------------------------------------

row_dend <- place_sibling_block_after(
  dend = row_dend,
  block_to_move = c(
    "Notch3",
    "Ptn"
  ),
  reference_block = c(
    "Pdgfrb",
    "Macf1"
  ),
  description = paste0(
    "rama Notch3/Ptn después ",
    "de la rama Pdgfrb/Macf1"
  )
)

# Final visual order after the two local rotations.
heatmap_gene_order <- labels(
  row_dend
)

final_gene_order <- heatmap_gene_order

cat(
  "\nOrden final del dendrograma:\n",
  "  - Bloques principales: early, intermediate, late.\n",
  "  - Rama Inmt/Fn1/Col15a1 situada después de Fbln2.\n",
  "  - Rama Nampt/Tnc situada después de Fmo2/Ddit4/Per1.\n",
  sep = ""
)

cat(
  "\nOrden de los bloques de genes:\n",
  "  Arriba: pseudotime bajo\n",
  "  Centro: pseudotime intermedio\n",
  "  Abajo: pseudotime alto\n\n",
  "Pseudotime mediano del pico por bloque:\n",
  sep = ""
)

print(
  cluster_peak_pseudotime[
    as.character(desired_cluster_order)
  ]
)

cat(
  "\nOrientación manual del dendrograma aplicada:\n",
  "  - Intercambiados los dos clusters principales.\n",
  "  - Intercambiados los dos subclusters del nuevo cluster superior.\n",
  sep = ""
)
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
# Cytokines, collagens and ECM heatmaps in fibroblast states
# ============================================================

figure2g_states <- fibro_states

figure2g_cells <- rownames(meta)[
  as.character(meta$Clustering) %in% figure2g_states
]

fibro_only <- subset(
  fibroblasts,
  cells = figure2g_cells
)

fibro_only$Clustering <- factor(
  as.character(
    fibro_only$Clustering
  ),
  levels = figure2g_states
)

DefaultAssay(fibro_only) <- "RNA"
fibro_only <- JoinLayers(
  fibro_only,
  assay = "RNA"
)

cat("\nCélulas para los heatmaps de Figure 2G:\n")
print(
  table(
    fibro_only$Clustering,
    useNA = "ifany"
  )
)

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
  fibro_only,
  cytokine_genes,
  "Cytokines / chemokines",
  "Figure2G_cytokines_chemokines"
)

make_figure2g_heatmap(
  fibro_only,
  collagen_genes,
  "Collagens",
  "Figure2G_collagens"
)

make_figure2g_heatmap(
  fibro_only,
  ecm_genes,
  "ECM / matrix-associated genes",
  "Figure2G_ECM"
)

capture.output(
  sessionInfo(),
  file = file.path(table_dir, "Figure2_sessionInfo.txt")
)

cat("\nFigure 2 finalizada.\n")