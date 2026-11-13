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
  library(dplyr)
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
  values <- vapply( parts, function(z) paste(z[-1], collapse = "="), character(1) )
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
preorder_cds_file <- file.path( objects_dir, "Figure2_fibroblast_before_ordering.rds" )
manual_root_nodes_file <- file.path( table_dir, "Figure2E_manual_root_pr_nodes_fibro_only.txt" )
root_selection_summary_file <- file.path( table_dir, "Figure2E_root_selection_used_fibro_only.tsv" )
# ============================================================
# Load annotated fibroblast object
# ============================================================
input_object_file <- file.path( objects_dir, "fibroblasts.Fibros.annotated.rds" )
if (!file.exists(input_object_file)) {
  stop( "No se encuentra el objeto de fibroblastos anotados:\n", input_object_file )
}
fibroblasts <- readRDS( input_object_file )
DefaultAssay(fibroblasts) <- "RNA"
required_metadata <- c( "fibro_annotation", "group", "orig.ident" )
missing_metadata <- setdiff( required_metadata, colnames(fibroblasts[[]]) )
if (length(missing_metadata) > 0) {
  stop( "Faltan columnas de metadatos en el objeto de entrada: ", paste(missing_metadata, collapse = ", ") )
}
fibro_states <- c( "Fibroblasts Universal", "Fibroblasts Activated 1", "Fibroblasts Activated 2" )
# Estados permitidos en el objeto de origen.
# Mesothelial puede estar presente en el RDS, pero se elimina
# antes de construir el objeto de Monocle.
source_states <- c( fibro_states, "Mesothelial" )
trajectory_states <- fibro_states
observed_states <- sort( unique( na.omit( as.character( fibroblasts$fibro_annotation ) ) ) )
unexpected_states <- setdiff( observed_states, source_states )
if (length(unexpected_states) > 0) {
  stop(
    "Se han encontrado anotaciones inesperadas en ",
    "'fibro_annotation':\n",
    paste(unexpected_states, collapse = "\n")
  )
}
if (anyNA(fibroblasts$fibro_annotation)) {
  stop( "El objeto contiene células con fibro_annotation = NA." )
}
fibroblasts$Clustering <- factor( as.character( fibroblasts$fibro_annotation ), levels = source_states )
# Metadatos que permiten comprobar el origen y el tipo de análisis
# cuando se reutiliza un CDS ya calculado.
fibroblasts$Figure2_source_object <- basename( input_object_file )
fibroblasts$Figure2_analysis <- "fibroblast_only_seurat_umap"
fibroblasts$Figure2_embedding <- "Seurat_umap"
trajectory_colors <- c(
  "Fibroblasts Universal" = "#FBCDC7",
  "Fibroblasts Activated 1" = "#F37B79",
  "Fibroblasts Activated 2" = "#6D1A1A"
)
fibro_colors <- trajectory_colors[ fibro_states ]
cat("\nObjeto de entrada de Figure 2:\n")
cat(input_object_file, "\n")
cat("\nAnotaciones presentes en el objeto de origen:\n")
print( table( fibroblasts$Clustering, useNA = "ifany" ) )
cat( "\nCélulas totales del objeto de origen: ", ncol(fibroblasts), "\n", sep = "" )
# ============================================================
# Build trajectory input with fibroblasts only
# ============================================================
meta <- fibroblasts[[]]
trajectory_cells <- rownames(meta)[ as.character(meta$Clustering) %in% fibro_states ]
trajectory_subset <- subset( fibroblasts, cells = trajectory_cells )
trajectory_subset$Clustering <- factor( as.character( trajectory_subset$Clustering ), levels = fibro_states )
Idents(trajectory_subset) <- "Clustering"
DefaultAssay(trajectory_subset) <- "RNA"
cat("\nCélulas utilizadas para construir la trayectoria:\n")
print( table( trajectory_subset$Clustering, useNA = "ifany" ) )
missing_fibro_states <- setdiff( fibro_states, unique( as.character( trajectory_subset$Clustering ) ) )
if (length(missing_fibro_states) > 0) {
  stop( "Faltan estados de fibroblastos: ", paste(missing_fibro_states, collapse = ", ") )
}
if (
  any( as.character( trajectory_subset$Clustering ) == "Mesothelial" )
) {
  stop( "Se han incorporado células Mesothelial por error." )
}
cat( "Número total de fibroblastos incorporados: ", ncol(trajectory_subset), "\n", sep = "" )
trajectory_subset <- JoinLayers( trajectory_subset, assay = "RNA" )
counts_matrix <- GetAssayData( trajectory_subset, assay = "RNA", layer = "counts" )
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
  stop( "El objeto de fibroblastos no contiene la reducción 'umap'." )
}
seurat_umap <- Embeddings( trajectory_subset, reduction = "umap" )
if (ncol(seurat_umap) < 2) {
  stop("La reducción UMAP de Seurat tiene menos de dos dimensiones.")
}
missing_umap_cells <- setdiff( colnames(counts_matrix), rownames(seurat_umap) )
if (length(missing_umap_cells) > 0) {
  stop( "Faltan células en el UMAP original de Seurat: ", paste(head(missing_umap_cells, 10), collapse = ", ") )
}
seurat_umap <- seurat_umap[ colnames(counts_matrix), 1:2, drop = FALSE ]
colnames(seurat_umap) <- c( "UMAP1", "UMAP2" )
expected_analysis <- "fibroblast_only_seurat_umap"
expected_embedding <- "Seurat_umap"
if (file.exists(preorder_cds_file)) {
  cat( "\nReutilizando el CDS previo a order_cells:\n", preorder_cds_file, "\n", sep = "" )
  cds_all <- readRDS( preorder_cds_file )
  required_cds_metadata <- c( "Figure2_source_object", "Figure2_analysis", "Figure2_embedding" )
  missing_cds_metadata <- setdiff( required_cds_metadata, colnames(colData(cds_all)) )
  if (length(missing_cds_metadata) > 0) {
    stop(
      "El CDS existente no contiene estos metadatos de control: ",
      paste(missing_cds_metadata, collapse = ", "),
      "\nElimina este archivo y vuelve a ejecutar el script:\n",
      preorder_cds_file
    )
  }
  cds_source <- unique( as.character( colData(cds_all)$Figure2_source_object ) )
  cds_source <- cds_source[ !is.na(cds_source) & nzchar(cds_source) ]
  expected_source <- basename( input_object_file )
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
  cds_analysis <- unique( as.character( colData(cds_all)$Figure2_analysis ) )
  cds_analysis <- cds_analysis[ !is.na(cds_analysis) & nzchar(cds_analysis) ]
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
  cds_embedding <- unique( as.character( colData(cds_all)$Figure2_embedding ) )
  cds_embedding <- cds_embedding[ !is.na(cds_embedding) & nzchar(cds_embedding) ]
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
  stored_umap <- SingleCellExperiment::reducedDim( cds_all, "UMAP" )
  missing_stored_cells <- setdiff( rownames(seurat_umap), rownames(stored_umap) )
  if (length(missing_stored_cells) > 0) {
    stop(
      "El CDS existente no contiene todas las células del UMAP original.\n",
      "Elimina este archivo y vuelve a ejecutar:\n",
      preorder_cds_file
    )
  }
  stored_umap <- stored_umap[ rownames(seurat_umap), 1:2, drop = FALSE ]
  umap_max_difference <- max( abs( unname(stored_umap) - unname(seurat_umap) ) )
  if (!isTRUE(all.equal(umap_max_difference, 0, tolerance = 0))) {
    stop(
      "El UMAP del CDS no es idéntico al UMAP original de Seurat.\n",
      "Diferencia máxima absoluta: ",
      umap_max_difference,
      "\nElimina este archivo y vuelve a ejecutar:\n",
      preorder_cds_file
    )
  }
  cat( "UMAP verificado: diferencia máxima absoluta = 0.\n" )
} else {
  cds_all <- new_cell_data_set(
expression_data = counts_matrix,
cell_metadata = cell_metadata,
gene_metadata = gene_metadata
  )
  cds_all <- estimate_size_factors( cds_all )
  # Import the exact Seurat coordinates.
  SingleCellExperiment::reducedDim( cds_all, "UMAP" ) <- seurat_umap[ colnames(cds_all), , drop = FALSE ]
  # Register the existing fibroblast annotations in Monocle's
  # internal structure. This does not calculate new clusters.
  seurat_fibro_states <- factor( as.character( colData(cds_all)$Clustering ), levels = fibro_states )
  names(seurat_fibro_states) <- colnames( cds_all )
  single_partition <- factor( rep( "1", ncol(cds_all) ), levels = "1" )
  names(single_partition) <- colnames( cds_all )
  cds_all@clusters[["UMAP"]] <- list( clusters = seurat_fibro_states, partitions = single_partition )
  imported_umap <- SingleCellExperiment::reducedDim( cds_all, "UMAP" )
  imported_umap <- imported_umap[ rownames(seurat_umap), 1:2, drop = FALSE ]
  umap_max_difference <- max( abs( unname(imported_umap) - unname(seurat_umap) ) )
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
  cat( "Anotaciones originales registradas desde fibro_annotation.\n", "No se ha calculado un clustering nuevo.\n" )
  set.seed(1234)
  cat( "Aprendiendo principal graph sobre el UMAP original de Seurat...\n" )
  cds_all <- learn_graph(
    cds_all,
use_partition = FALSE,
learn_graph_control = list( minimal_branch_len = 40, prune_graph = FALSE )
  )
  saveRDS( cds_all, preorder_cds_file )
  cat( "CDS previo a order_cells guardado en: ", preorder_cds_file, "\n", sep = "" )
}
cat("\nFibroblastos presentes en el CDS usado para seleccionar el root:\n")
print( table( colData(cds_all)$Clustering, useNA = "ifany" ) )
if (
  any( as.character( colData(cds_all)$Clustering ) == "Mesothelial" )
) {
  stop( "El CDS contiene células Mesothelial. ", "Borra el RDS previo y vuelve a ejecutar el script." )
}
cds_states <- unique( as.character( colData(cds_all)$Clustering ) )
cds_states <- cds_states[ !is.na(cds_states) ]
unexpected_cds_states <- setdiff( cds_states, fibro_states )
if (length(unexpected_cds_states) > 0) {
  stop(
    "El CDS contiene poblaciones distintas de los tres ",
    "estados de fibroblastos: ",
    paste(unexpected_cds_states, collapse = ", ")
  )
}
graph_vertex_names <- igraph::V( principal_graph(cds_all)[["UMAP"]] )$name
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
  quit( save = "no", status = 0 )
}
root_nodes_used <- unique( trimws( readLines( manual_root_nodes_file, warn = FALSE ) ) )
root_nodes_used <- root_nodes_used[nzchar(root_nodes_used)]
if (length(root_nodes_used) == 0) {
  stop( "El archivo de nodos raíz está vacío: ", manual_root_nodes_file )
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
cds_all <- order_cells( cds_all, reduction_method = "UMAP", root_pr_nodes = root_nodes_used )
pseudotime_all <- pseudotime( cds_all, reduction_method = "UMAP" )
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
  data.frame( Selection_method = root_selection_method, Root_pr_node = root_nodes_used, stringsAsFactors = FALSE ),
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
pseudotime_fibro_raw <- pseudotime_all[ fibro_cells ]
names(pseudotime_fibro_raw) <- fibro_cells
# Rescale fibroblast pseudotime to a fixed 0-20 range.
# This changes only the numerical display scale, not cell ordering.
finite_pseudotime <- is.finite(pseudotime_fibro_raw)
if (!any(finite_pseudotime)) {
  stop("No hay células de fibroblastos con pseudotime finito.")
}
pseudotime_limits_raw <- range( pseudotime_fibro_raw[finite_pseudotime], finite = TRUE )
if (diff(pseudotime_limits_raw) == 0) {
  stop("El pseudotime de fibroblastos no presenta variación.")
}
pseudotime_fibro <- rep( Inf, length(pseudotime_fibro_raw) )
names(pseudotime_fibro) <- fibro_cells
pseudotime_fibro[finite_pseudotime] <- (
  pseudotime_fibro_raw[finite_pseudotime] -
    pseudotime_limits_raw[1]
) / diff(pseudotime_limits_raw) * 20
colData(cds_fibro)$monocle3_pseudotime_raw <-
  pseudotime_fibro_raw
colData(cds_fibro)$monocle3_pseudotime <-
  pseudotime_fibro
saveRDS( cds_fibro, file.path( objects_dir, "Figure2_fibroblast_monocle3_cds.rds" ) )
fibro_subset_plot <- trajectory_subset
fibro_subset_plot$Clustering <- factor( as.character( fibro_subset_plot$Clustering ), levels = fibro_states )
Idents(fibro_subset_plot) <- "Clustering"
cat("\nCélulas usadas para calcular el pseudotime:\n")
print( table( colData(cds_fibro)$Clustering, useNA = "ifany" ) )
if (
  any( as.character( colData(cds_fibro)$Clustering ) == "Mesothelial" )
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
monocle_umap_all <- SingleCellExperiment::reducedDim( cds_all, "UMAP" )
rownames(monocle_umap_all) <- colnames(cds_all)
colnames(monocle_umap_all) <- c("UMAP1", "UMAP2")
fibro_pseudotime_table <- data.frame(
Cell = fibro_cells,
UMAP1 = monocle_umap_all[fibro_cells, "UMAP1"],
UMAP2 = monocle_umap_all[fibro_cells, "UMAP2"],
Clustering = factor( as.character(colData(cds_all)[fibro_cells, "Clustering"]), levels = fibro_states ),
Pseudotime_raw = as.numeric( pseudotime_fibro_raw[fibro_cells] ),
Pseudotime = as.numeric( pseudotime_fibro[fibro_cells] ),
stringsAsFactors = FALSE
)
fibro_pseudotime_table$Pseudotime_raw[ !is.finite(fibro_pseudotime_table$Pseudotime_raw) ] <- NA_real_
fibro_pseudotime_table$Pseudotime[ !is.finite(fibro_pseudotime_table$Pseudotime) ] <- NA_real_
# ------------------------------------------------------------
# Principal graph used for plotting
#
# The graph was learned using only fibroblast cells.
# ------------------------------------------------------------
principal_graph_all <- principal_graph( cds_all )[["UMAP"]]
principal_graph_aux_all <- principal_graph_aux( cds_all )[["UMAP"]]
if (is.null(principal_graph_all)) {
  stop("El CDS no contiene principal graph para UMAP.")
}
principal_graph_fibro_plot <- principal_graph_all
graph_vertex_names <- igraph::V( principal_graph_fibro_plot )$name
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
  stop( "No se reconoce la estructura de coordenadas de dp_mst." )
}
if (is.null(rownames(principal_point_coordinates))) {
  rownames(principal_point_coordinates) <- graph_vertex_names
}
colnames(principal_point_coordinates)[1:2] <- c( "UMAP1", "UMAP2" )
missing_graph_coordinates <- setdiff(
  igraph::V(principal_graph_fibro_plot)$name,
  rownames(principal_point_coordinates)
)
if (length(missing_graph_coordinates) > 0) {
  stop( "Faltan coordenadas para estos principal points: ", paste(missing_graph_coordinates, collapse = ", ") )
}
fibro_graph_edges <- igraph::as_data_frame( principal_graph_fibro_plot, what = "edges" )
fibro_graph_segments <- data.frame(
x = principal_point_coordinates[ fibro_graph_edges$from, "UMAP1" ],
y = principal_point_coordinates[ fibro_graph_edges$from, "UMAP2" ],
xend = principal_point_coordinates[ fibro_graph_edges$to, "UMAP1" ],
yend = principal_point_coordinates[ fibro_graph_edges$to, "UMAP2" ],
stringsAsFactors = FALSE
)
fibro_graph_segments <- fibro_graph_segments[ complete.cases(fibro_graph_segments), , drop = FALSE ]
# ------------------------------------------------------------
# Common coordinate limits for both panels
# ------------------------------------------------------------
common_xlim <- range(
  c( fibro_pseudotime_table$UMAP1, fibro_graph_segments$x, fibro_graph_segments$xend ),
finite = TRUE
)
common_ylim <- range(
  c( fibro_pseudotime_table$UMAP2, fibro_graph_segments$y, fibro_graph_segments$yend ),
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
Figure2E_clusters <- ggplot( fibro_pseudotime_table, aes( x = UMAP1, y = UMAP2, color = Clustering ) ) +
  geom_point( size = 0.8, alpha = 1 ) +
  scale_color_manual( values = fibro_colors, breaks = fibro_states, drop = FALSE, name = NULL ) +
  guides(
color = guide_legend( direction = "vertical", ncol = 1, byrow = TRUE, override.aes = list( size = 2.2, alpha = 1 ) )
  ) +
  labs( title = "Fibroblast states", x = "UMAP1", y = "UMAP2" ) +
  coord_fixed( xlim = common_xlim, ylim = common_ylim, ratio = 0.75, clip = "off" ) +
  theme_classic(base_size = 12) +
  theme(
plot.title = element_text( hjust = 0.5, face = "bold" ),
legend.position = "right",
legend.justification = "center",
legend.direction = "vertical",
legend.box = "vertical",
legend.margin = margin( t = 0, r = 0, b = 0, l = 4 ),
plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )
# ------------------------------------------------------------
# Right panel: fibroblast pseudotime
# ------------------------------------------------------------
Figure2E_pseudotime <- ggplot( fibro_pseudotime_table, aes( x = UMAP1, y = UMAP2 ) ) +
  geom_segment(
data = fibro_graph_segments,
    aes( x = x, y = y, xend = xend, yend = yend ),
inherit.aes = FALSE,
color = "black",
linewidth = 0.6,
lineend = "round"
  ) +
  geom_point( aes(color = Pseudotime), size = 0.8, alpha = 1 ) +
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
  labs( title = "Pseudotime", x = "UMAP1", y = "UMAP2" ) +
  coord_fixed( xlim = common_xlim, ylim = common_ylim, ratio = 0.75, clip = "off" ) +
  theme_classic(base_size = 12) +
  theme(
plot.title = element_text( hjust = 0.5, face = "bold" ),
legend.position = "right",
legend.justification = "center",
legend.direction = "vertical",
legend.box = "vertical",
legend.margin = margin( t = 0, r = 0, b = 0, l = 4 ),
plot.margin = margin(5.5, 5.5, 5.5, 7)
  )
# ------------------------------------------------------------
# Combine and save
# ------------------------------------------------------------
Figure2E <- Figure2E_clusters + Figure2E_pseudotime +
  plot_layout( ncol = 2, widths = figure2e_panel_widths, guides = "keep" )
ggsave(
filename = file.path( figure_dir, "Figure2E_fibroblast_states_pseudotime.pdf" ),
plot = Figure2E,
width = figure2e_width,
height = figure2e_height,
units = "in",
useDingbats = FALSE,
limitsize = FALSE
)
ggsave(
filename = file.path( figure_dir, "Figure2E_fibroblast_states_pseudotime.png" ),
plot = Figure2E,
width = figure2e_width,
height = figure2e_height,
units = "in",
dpi = 300,
limitsize = FALSE
)
Figure2E_table <- fibro_pseudotime_table[
  ,
  c( "Cell", "UMAP1", "UMAP2", "Clustering", "Pseudotime_raw", "Pseudotime" )
]
write.table(
  Figure2E_table,
  file.path( table_dir, "Figure2E_fibroblast_states_pseudotime.tsv" ),
sep = "\t",
quote = FALSE,
row.names = FALSE
)
cat(
  "\nFigure 2E finalizada.\n",
  "PDF: ",
  file.path( figure_dir, "Figure2E_fibroblast_states_pseudotime.pdf" ),
  "\nPNG: ",
  file.path( figure_dir, "Figure2E_fibroblast_states_pseudotime.png" ),
  "\n",
sep = ""
)
# ============================================================
# Figure 2F: dynamic genes along fibroblast pseudotime
# The trajectory and pseudotime contain fibroblasts only.
# ============================================================
# ============================================================
# Figure 2F: graph_test
#
# graph_test() is the slowest recurrent step in this script.
# Reuse the complete previous result when it is compatible with
# the current CDS and newer than the objects from which it was
# generated. Otherwise calculate it again.
# ============================================================
graph_test_table_file <- file.path( table_dir, "Figure2F_graph_test_all.tsv" )
graph_test_table <- NULL
reuse_graph_test <- FALSE
if (file.exists(graph_test_table_file)) {
  graph_test_is_current <-
    file.mtime(graph_test_table_file) >=
      file.mtime(input_object_file) &&
    file.mtime(graph_test_table_file) >=
      file.mtime(preorder_cds_file)
  if (isTRUE(graph_test_is_current)) {
    previous_graph_test <- tryCatch(
      read.delim( graph_test_table_file, check.names = FALSE, stringsAsFactors = FALSE ),
error = function(e) NULL
    )
    required_graph_test_columns <- c( "Gene", "status", "q_value", "morans_I" )
    if (
      !is.null(previous_graph_test) &&
        all( required_graph_test_columns %in% colnames(previous_graph_test) ) &&
        nrow(previous_graph_test) ==
          nrow(cds_fibro) &&
        setequal( previous_graph_test$Gene, rownames(cds_fibro) )
    ) {
      graph_test_table <- previous_graph_test
      reuse_graph_test <- TRUE
      cat( "\nReutilizando graph_test previo: ", graph_test_table_file, "\n", sep = "" )
    }
  }
}
if (!reuse_graph_test) {
  cat( "\nEjecutando graph_test en fibroblastos...\n" )
  modulated_genes <- graph_test( cds_fibro, neighbor_graph = "principal_graph", cores = n_cores )
  graph_test_table <- as.data.frame( modulated_genes )
  graph_test_table$Gene <- rownames( graph_test_table )
  write.table( graph_test_table, graph_test_table_file, sep = "\t", quote = FALSE, row.names = FALSE )
  cat( "graph_test actualizado y guardado en: ", graph_test_table_file, "\n", sep = "" )
}
sig_genes <- graph_test_table[
  graph_test_table$status == "OK" &
    graph_test_table$q_value < 0.05 &
    is.finite(graph_test_table$morans_I),
  ,
  drop = FALSE
]
sig_genes <- sig_genes[ order(sig_genes$morans_I, decreasing = TRUE), , drop = FALSE ]
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
excluded_genes <- c( "Mapkapk3", "Dpep1", "Mgp", "Cfh", "Spon2", "C3", "Txnip" )
manual_genes <- unique( setdiff( manual_genes, excluded_genes ) )
missing_manual_genes <- setdiff( manual_genes, rownames(rowData(cds_fibro)) )
if (length(missing_manual_genes) > 0) {
  stop( "Estos genes manuales no están presentes en cds_fibro: ", paste(missing_manual_genes, collapse = ", ") )
}
final_gene_number <- 60L
n_automatic_genes <- final_gene_number -
  length(manual_genes)
if (n_automatic_genes < 0) {
  stop( "La lista manual contiene más de ", final_gene_number, " genes." )
}
automatic_candidates <- sig_genes[ !sig_genes$Gene %in% c( manual_genes, excluded_genes ), , drop = FALSE ]
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
automatic_genes <- head( automatic_candidates$Gene, n_automatic_genes )
valid_genes <- unique( c( automatic_genes, manual_genes ) )
valid_genes <- intersect( valid_genes, rownames(rowData(cds_fibro)) )
valid_genes <- setdiff( valid_genes, excluded_genes )
if (length(valid_genes) != final_gene_number) {
  stop( "La selección final contiene ", length(valid_genes), " genes en lugar de ", final_gene_number, "." )
}
if (!all(manual_genes %in% valid_genes)) {
  stop( "No se conservaron todos los genes manuales en la selección final." )
}
top_genes <- graph_test_table[ match( valid_genes, graph_test_table$Gene ), , drop = FALSE ]
top_genes$Gene <- valid_genes
top_genes$Selection <- ifelse( top_genes$Gene %in% manual_genes, "Manual", "Moran_I" )
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
ordered_cells <- finite_cells[ order(pseudotime_fibro[finite_cells]) ]
ordered_pseudotime <- pseudotime_fibro[ordered_cells]
pt_matrix <- exprs(cds_fibro)[ valid_genes, ordered_cells, drop = FALSE ]
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
expression_color <- colorRamp2( seq(-2, 2, length.out = 11), rev(brewer.pal(11, "Spectral")) )
group_colors <- c( "Naive" = "#FA5858", "Tumor" = "#04B486" )
sample_colors <- c( "WT" = "#A6CEE3", "CD169DTR" = "#1F78B4", "CCR2KO" = "#08306B" )
heatmap_metadata <- as.data.frame(colData(cds_fibro))[ ordered_cells, , drop = FALSE ]
cell_clusters <- factor( as.character(heatmap_metadata$Clustering), levels = fibro_states )
cell_group <- factor( as.character(heatmap_metadata$group), levels = c("Naive", "Tumor") )
cell_sample <- factor( as.character(heatmap_metadata$orig.ident), levels = c("WT", "CD169DTR", "CCR2KO") )
pseudotime_color <- colorRamp2( c(0, 20), c("#F7FCF5", "#238B45") )
column_annotation <- HeatmapAnnotation(
Cluster = cell_clusters,
Group = cell_group,
Sample = cell_sample,
Pseudotime = ordered_pseudotime,
col = list( Cluster = fibro_colors, Group = group_colors, Sample = sample_colors, Pseudotime = pseudotime_color ),
show_annotation_name = TRUE,
annotation_name_side = "left",
annotation_legend_param = list( Pseudotime = list( at = seq(0, 20, by = 5) ) ),
simple_anno_size = unit(4, "mm")
)
# ============================================================
# Hierarchical clustering of genes
# Ward.D2 dendrogram with manual branch orientation
# ============================================================
row_hc <- hclust( dist(pt_matrix), method = "ward.D2" )
row_dend <- as.dendrogram( row_hc )
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
row_clusters <- cutree( row_hc, k = 3 )
if (length(unique(row_clusters)) != 3L) {
  stop( "No se identificaron exactamente tres bloques de genes." )
}
# Pseudotime at which each smoothed gene profile reaches
# its maximum expression.
gene_peak_index <- apply( pt_matrix, 1, which.max )
gene_peak_pseudotime <- ordered_pseudotime[ gene_peak_index ]
names(gene_peak_pseudotime) <- rownames( pt_matrix )
# Representative peak pseudotime for each hierarchical block.
cluster_peak_pseudotime <- tapply( gene_peak_pseudotime[names(row_clusters)], row_clusters, median, na.rm = TRUE )
if (
  length(cluster_peak_pseudotime) != 3L ||
    any(!is.finite(cluster_peak_pseudotime))
) {
  stop( "No se pudo calcular el pseudotime de pico de los tres bloques." )
}
# Early -> intermediate -> late.
desired_cluster_order <- as.integer( names( sort( cluster_peak_pseudotime, decreasing = FALSE ) ) )
current_gene_order <- labels( row_dend )
# Preserve the internal dendrogram order of genes within
# each block.
desired_gene_order <- unlist(
  lapply(
    desired_cluster_order,
    function(cluster_id) {
      current_gene_order[ row_clusters[current_gene_order] == cluster_id ]
    }
  ),
use.names = FALSE
)
if (
  length(desired_gene_order) != nrow(pt_matrix) ||
    anyDuplicated(desired_gene_order) ||
    !setequal(desired_gene_order, rownames(pt_matrix))
) {
  stop( "El nuevo orden no contiene exactamente todos los genes." )
}
# Safe branch rotation: dendrogram attributes are preserved.
row_dend <- dendextend::rotate( row_dend, order = desired_gene_order )
# Confirm that the three major blocks were ordered correctly.
major_block_gene_order <- labels( row_dend )
if (!identical(
  as.character(major_block_gene_order),
  as.character(desired_gene_order)
)) {
  stop( "No se pudo obtener el orden general ", "early-intermediate-late solicitado." )
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
  dend_labels <- labels( dend )
  missing_genes <- setdiff( c( block_to_move, reference_block ), dend_labels )
  if (length(missing_genes) > 0) {
    stop( "No se encuentran estos genes para la rotación '", description, "': ", paste(missing_genes, collapse = ", ") )
  }
  find_sibling_split <- function(node) {
    if (is.leaf(node)) {
      return(NULL)
    }
    left_genes <- labels( node[[1]] )
    right_genes <- labels( node[[2]] )
    move_in_left <- all( block_to_move %in% left_genes )
    move_in_right <- all( block_to_move %in% right_genes )
    reference_in_left <- all( reference_block %in% left_genes )
    reference_in_right <- all( reference_block %in% right_genes )
    # The two requested gene groups are located in
    # opposite immediate branches of this node.
    if (
      move_in_left &&
        reference_in_right
    ) {
      return( list( node_genes = labels(node), block_to_move_genes = left_genes, reference_branch_genes = right_genes ) )
    }
    if (
      move_in_right &&
        reference_in_left
    ) {
      return( list( node_genes = labels(node), block_to_move_genes = right_genes, reference_branch_genes = left_genes ) )
    }
    target_genes <- unique( c( block_to_move, reference_block ) )
    # Continue down the only branch containing both groups.
    if (all(target_genes %in% left_genes)) {
      return( find_sibling_split( node[[1]] ) )
    }
    if (all(target_genes %in% right_genes)) {
      return( find_sibling_split( node[[2]] ) )
    }
    NULL
  }
  sibling_split <- find_sibling_split( dend )
  if (is.null(sibling_split)) {
    stop( "No se encontraron ramas hermanas compatibles para: ", description )
  }
  current_order <- labels( dend )
  node_positions <- which( current_order %in% sibling_split$node_genes )
  if (
    length(node_positions) == 0L ||
      !identical( current_order[node_positions], sibling_split$node_genes )
  ) {
    stop( "No se pudo localizar de forma continua la rama para: ", description )
  }
  # Desired local order:
  # reference branch first, moved branch afterwards.
  desired_local_order <- c( sibling_split$reference_branch_genes, sibling_split$block_to_move_genes )
  desired_full_order <- current_order
  desired_full_order[node_positions] <-
    desired_local_order
  # Nothing to do if it is already correctly oriented.
  if (identical(
    current_order,
    desired_full_order
  )) {
    cat( "La rotación ya estaba aplicada: ", description, "\n", sep = "" )
    return(dend)
  }
  rotated_dend <- dendextend::rotate( dend, order = desired_full_order )
  if (!identical(
    as.character(labels(rotated_dend)),
    as.character(desired_full_order)
  )) {
    stop( "No se pudo aplicar la rotación local: ", description )
  }
  cat( "Rotación local aplicada: ", description, "\n", sep = "" )
  rotated_dend
}
# ------------------------------------------------------------
# 1. Move the Inmt/Fn1/Col15a1 branch below the
#    neighbouring early-pseudotime branch ending in Fbln2.
# ------------------------------------------------------------
row_dend <- place_sibling_block_after(
dend = row_dend,
block_to_move = c( "Inmt", "Fn1", "Col15a1" ),
reference_block = c( "Fbln2" ),
description = paste0( "rama Inmt/Fn1/Col15a1 después ", "de la rama que termina en Fbln2" )
)
# ------------------------------------------------------------
# 2. Move the Nampt/Tnc branch below the
#    Fmo2/Ddit4/Per1 branch.
# ------------------------------------------------------------
row_dend <- place_sibling_block_after(
dend = row_dend,
block_to_move = c( "Nampt", "Tnc" ),
reference_block = c( "Fmo2", "Ddit4", "Per1" ),
description = paste0( "rama Nampt/Tnc después ", "de la rama Fmo2/Ddit4/Per1" )
)
# ------------------------------------------------------------
# 3. Move the Notch3/Ptn branch below the
#    Pdgfrb/Macf1 branch.
# ------------------------------------------------------------
row_dend <- place_sibling_block_after(
dend = row_dend,
block_to_move = c( "Notch3", "Ptn" ),
reference_block = c( "Pdgfrb", "Macf1" ),
description = paste0( "rama Notch3/Ptn después ", "de la rama Pdgfrb/Macf1" )
)
# Final visual order after the two local rotations.
heatmap_gene_order <- labels( row_dend )
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
print( cluster_peak_pseudotime[ as.character(desired_cluster_order) ] )
cat(
  "\nOrientación manual del dendrograma aplicada:\n",
  "  - Intercambiados los dos clusters principales.\n",
  "  - Intercambiados los dos subclusters del nuevo cluster superior.\n",
sep = ""
)
cat( "\nGenes incluidos en el clustering jerárquico: ", length(final_gene_order), "\n", sep = "" )
Figure2F <- Heatmap(
  pt_matrix,
name = "Scaled Expression",
col = expression_color,
cluster_rows = row_dend,
cluster_columns = FALSE,
show_row_dend = TRUE,
show_row_names = TRUE,
show_column_names = FALSE,
row_names_gp = gpar( fontsize = 12, fontface = "italic" ),
top_annotation = column_annotation,
column_title = "Pseudotime",
column_title_gp = gpar( fontsize = 16, fontface = "bold" ),
use_raster = TRUE,
raster_quality = 3
)
pdf( file = file.path( figure_dir, "Figure2F_pseudotime_heatmap.pdf" ), width = 18, height = 12, useDingbats = FALSE )
draw( Figure2F, heatmap_legend_side = "right", annotation_legend_side = "right" )
dev.off()
png( filename = file.path( figure_dir, "Figure2F_pseudotime_heatmap.png" ), width = 4000, height = 2400, res = 250 )
draw( Figure2F, heatmap_legend_side = "right", annotation_legend_side = "right" )
dev.off()
selected_gene_table <- top_genes[ match(heatmap_gene_order, top_genes$Gene), , drop = FALSE ]
selected_gene_table$Heatmap_order <- match( selected_gene_table$Gene, final_gene_order )
selected_gene_table <- selected_gene_table[ order(selected_gene_table$Heatmap_order), , drop = FALSE ]
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
# Cytokines, ECM / collagen and metabolism
#
# Three vertical heatmaps arranged in three columns.
# The same number of genes is displayed in all panels.
#
# Cytokines / chemokines / growth factors:
#   all 22 explicitly named cytokines / growth factors from the manuscript
#   text are forced into the heatmap. Two additional TNF/TNFSF-family genes
#   are selected from Fibroblasts Activated 1 to represent the phrase
#   "Tnfa-related family genes". The cytokine panel therefore contains
#   24 genes in total.
#
# ECM / collagen:
#   5 distinctive genes are selected for each fibroblast state,
#   and the remaining 9 positions are filled by global differential ranking.
#
# Metabolism:
#   genes are ranked specifically by differential expression
#   of Fibroblasts Activated 2 versus the other two fibroblast
#   states, because this population shows the strongest
#   metabolic phenotype.
# ============================================================
figure2g_states <- fibro_states
figure2g_cells <- rownames(meta)[as.character(meta$Clustering) %in% figure2g_states]
fibro_only <- subset(fibroblasts, cells = figure2g_cells)
fibro_only$Clustering <- factor(as.character(fibro_only$Clustering), levels = figure2g_states)
DefaultAssay(fibro_only) <- "RNA"
fibro_only <- JoinLayers(fibro_only, assay = "RNA")
Idents(fibro_only) <- "Clustering"
cat("\nCélulas para los heatmaps de Figure 2G:\n")
print(table(fibro_only$Clustering, useNA = "ifany"))
# ============================================================
# Same number of genes in every heatmap
# ============================================================
figure2g_n_genes <- 24L
# ============================================================
# 1. Cytokines / chemokines
# ============================================================
cytokine_genes <- c(
  "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl7", "Ccl8", "Ccl11",
  "Ccl12", "Ccl17", "Ccl19", "Ccl20", "Ccl21", "Ccl22",
  "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl5", "Cxcl9", "Cxcl10",
  "Cxcl11", "Cxcl12", "Cxcl13", "Cxcl14", "Cxcl16",
  "Il1a", "Il1b", "Il6", "Il7", "Il11", "Il15", "Il18",
  "Il33", "Il34", "Tnf", "Tnfsf10", "Tnfsf11", "Tnfsf12",
  "Tnfsf13b", "Tnfsf14", "Lta", "Tgfb1", "Tgfb2", "Tgfb3",
  "Lif", "Osm", "Csf1", "Csf2", "Csf3", "Kitl", "Flt3l",
  "Areg", "Ereg", "Gdf15", "Vegfa", "Vegfb", "Pdgfa", "Pdgfb",
  "Fgf2", "Ptgs2"
)

# Genes explicitly named in the manuscript text.
#
# The first block corresponds to the first sentence of the supplied text;
# the second block to Activated fibroblast 1; the third block to
# Activated fibroblast 2.
manuscript_cytokine_genes <- c(
  "Cxcl9", "Cxcl10", "Cxcl5", "Cxcl12",
  "Csf1", "Il34", "Vegfa", "Tgfb1", "Tgfb2",
  "Osm", "Ccl5",
  "Ccl19", "Il15",
  "Cxcl16", "Cxcl13", "Ccl2", "Cxcl1", "Cxcl2", "Ccl3",
  "Il7", "Il11", "Il33"
)

manuscript_cytokine_groups <- c(
  "Cxcl9" = "Manuscript first block",
  "Cxcl10" = "Manuscript first block",
  "Cxcl5" = "Manuscript first block",
  "Cxcl12" = "Manuscript first block",
  "Csf1" = "Manuscript first block",
  "Il34" = "Manuscript first block",
  "Vegfa" = "Manuscript first block",
  "Tgfb1" = "Manuscript first block",
  "Tgfb2" = "Manuscript first block",
  "Osm" = "Activated fibroblast 1",
  "Ccl5" = "Activated fibroblast 1",
  "Ccl19" = "Activated fibroblast 2",
  "Il15" = "Activated fibroblast 2",
  "Cxcl16" = "Activated fibroblast 2",
  "Cxcl13" = "Activated fibroblast 2",
  "Ccl2" = "Activated fibroblast 2",
  "Cxcl1" = "Activated fibroblast 2",
  "Cxcl2" = "Activated fibroblast 2",
  "Ccl3" = "Activated fibroblast 2",
  "Il7" = "Activated fibroblast 2",
  "Il11" = "Activated fibroblast 2",
  "Il33" = "Activated fibroblast 2"
)

activated1_tnf_family_candidates <- c(
  "Tnf", "Tnfsf10", "Tnfsf11", "Tnfsf12",
  "Tnfsf13b", "Tnfsf14", "Lta"
)
# ============================================================
# 2. ECM / matrix-associated genes
# ============================================================
ecm_genes <- c(
  "Fn1", "Fbn1", "Fbn2", "Postn", "Tnc", "Thbs1", "Thbs2",
  "Sparc", "Sparcl1", "Vtn", "Nid1", "Nid2", "Emilin1", "Emilin2",
  "Mfap2", "Mfap4", "Mfap5", "Eln", "Efemp1", "Efemp2",
  "Fbln1", "Fbln2", "Fbln5", "Fbln7", "Dcn", "Lum", "Bgn",
  "Fmod", "Aspn", "Ogn", "Prelp", "Hspg2", "Sdc1", "Sdc2",
  "Sdc3", "Sdc4", "Gpc1", "Gpc3", "Gpc4", "Gpc6",
  "Ccn1", "Ccn2", "Ccn3", "Ccn4", "Ccn5",
  "Lox", "Loxl1", "Loxl2", "Loxl3", "Loxl4",
  "Serpine1", "Serpinh1", "Pcolce", "Pcolce2",
  "Plod1", "Plod2", "Plod3",
  "P4ha1", "P4ha2", "P4ha3",
  "Timp1", "Timp2", "Timp3", "Timp4",
  "Acta2", "Tagln", "Myl9", "Myh11", "Cnn1"
)
collagen_genes <- c(
  "Col1a1", "Col1a2", "Col2a1", "Col3a1",
  "Col4a1", "Col4a2", "Col4a3", "Col4a4", "Col4a5", "Col4a6",
  "Col5a1", "Col5a2", "Col5a3",
  "Col6a1", "Col6a2", "Col6a3", "Col6a4", "Col6a5", "Col6a6",
  "Col7a1", "Col8a1", "Col8a2",
  "Col9a1", "Col9a2", "Col9a3",
  "Col10a1", "Col11a1", "Col11a2", "Col12a1", "Col13a1",
  "Col14a1", "Col15a1", "Col16a1", "Col17a1", "Col18a1",
  "Col19a1", "Col20a1", "Col21a1", "Col22a1", "Col23a1",
  "Col24a1", "Col25a1", "Col26a1", "Col27a1", "Col28a1"
)
ecm_collagen_genes <- unique(c( ecm_genes, collagen_genes ))
# ============================================================
# 3. Metabolism
#
# Broad candidate set:
# glycolysis, pyruvate metabolism, TCA/OXPHOS, fatty-acid
# metabolism, lipid synthesis, pentose-phosphate pathway,
# amino-acid metabolism, NAD/redox and metabolic-stress
# regulators.
#
# Selection is data-driven and specifically prioritizes genes
# altered in Fibroblasts Activated 2 versus the other states.
# ============================================================
metabolic_genes <- unique(c(
  "Ddit4", "Nampt", "Txnip", "Klf15", "Per2",
  "Hif1a", "Mtor", "Rptor", "Rictor",
  "Prkaa1", "Prkaa2", "Sirt1", "Sirt3",
  "Ppara", "Ppard", "Ppargc1a",
  "Slc2a1", "Slc2a3", "Hk1", "Hk2",
  "Gpi1", "Pfkl", "Pfkm", "Pfkp",
  "Aldoa", "Aldoc", "Tpi1", "Gapdh",
  "Pgk1", "Pgam1", "Eno1", "Eno2",
  "Pkm", "Ldha", "Ldhb",
  "Pdk1", "Pdk2", "Pdk3", "Pdk4",
  "Mpc1", "Mpc2", "Pdha1", "Pdhb", "Dlat",
  "Cs", "Aco2", "Idh1", "Idh2", "Idh3a",
  "Ogdh", "Suclg1", "Suclg2",
  "Sdha", "Sdhb", "Sdhc", "Sdhd",
  "Fh1", "Mdh1", "Mdh2",
  "Ndufa9", "Ndufs1", "Ndufs2", "Ndufs3", "Ndufv1",
  "Uqcrc1", "Uqcrc2", "Cyc1",
  "Cox4i1", "Cox5a", "Cox5b", "Cox6a1",
  "Atp5f1a", "Atp5f1b", "Atp5mc1",
  "Cd36", "Lpl", "Fabp3", "Fabp4", "Fabp5",
  "Cpt1a", "Cpt2", "Acadl", "Acadm", "Acads",
  "Hadha", "Hadhb", "Echs1", "Acox1", "Acaa2",
  "Acly", "Acaca", "Acacb", "Fasn", "Scd1", "Elovl6",
  "Hmgcs1", "Hmgcr", "Sqle", "Fdps",
  "G6pd", "Pgd", "Tkt", "Taldo1",
  "Phgdh", "Psat1", "Psph",
  "Shmt1", "Shmt2", "Mthfd1", "Mthfd2",
  "Gls", "Glul", "Got1", "Got2", "Gpt2",
  "Slc1a5", "Slc7a5", "Bcat1", "Bcat2",
  "Asns", "Ass1",
  "Me1", "Nqo1", "Gclc", "Gclm",
  "Gpx1", "Gpx4", "Sod1", "Sod2",
  "Prdx1", "Prdx2", "Txn1", "Txnrd1"
))
cytokine_present <- intersect(cytokine_genes, rownames(fibro_only))
ecm_present <- intersect(ecm_collagen_genes, rownames(fibro_only))
metabolic_present <- intersect(metabolic_genes, rownames(fibro_only))

missing_manuscript_cytokines <- setdiff(
  manuscript_cytokine_genes,
  cytokine_present
)

if (length(missing_manuscript_cytokines) > 0) {
  stop(
    "Figure 2G: estos genes citados explícitamente en el manuscrito ",
    "no están disponibles en el objeto y no pueden omitirse silenciosamente: ",
    paste(missing_manuscript_cytokines, collapse = ", ")
  )
}
if (length(cytokine_present) < figure2g_n_genes) {
  stop("Figure 2G: no hay suficientes genes de cytokines presentes.")
}
if (length(ecm_present) < figure2g_n_genes) {
  stop("Figure 2G: no hay suficientes genes de ECM / collagen presentes.")
}
if (length(metabolic_present) < figure2g_n_genes) {
  stop("Figure 2G: no hay suficientes genes metabólicos presentes.")
}
# ============================================================
# Cache DE and average expression for Figure 2G
# ============================================================
figure2g_cache_file <- file.path( table_dir, "Figure2G_three_panels_v2_cache.rds" )
input_info <- file.info(input_object_file)
figure2g_cache_signature <- list(
input_size = unname(input_info$size),
input_mtime = as.numeric(input_info$mtime),
cells = colnames(fibro_only),
states = figure2g_states,
cytokine_genes = sort(cytokine_present),
ecm_collagen_genes = sort(ecm_present),
metabolic_genes = sort(metabolic_present)
)
figure2g_cache <- NULL
if (file.exists(figure2g_cache_file)) {
  candidate_cache <- tryCatch( readRDS(figure2g_cache_file), error = function(e) NULL )
  if (
    !is.null(candidate_cache) &&
      identical(candidate_cache$signature, figure2g_cache_signature)
  ) {
    figure2g_cache <- candidate_cache
    cat( "\nReutilizando cache de Figure 2G:\n", figure2g_cache_file, "\n", sep = "" )
  }
}
if (is.null(figure2g_cache)) {
  all_average_genes <- unique(c( cytokine_present, ecm_present, metabolic_present ))
  cat( "\nFigure 2G: calculando expresión media para ", length(all_average_genes), " genes candidatos...\n", sep = "" )
  avg_all <- AverageExpression(
object = fibro_only,
assays = "RNA",
features = all_average_genes,
group.by = "Clustering",
layer = "data",
verbose = FALSE
  )[["RNA"]]
  avg_all <- as.matrix(avg_all)
  storage.mode(avg_all) <- "double"
  state_keys <- make.names(figure2g_states)
  column_keys <- make.names(colnames(avg_all))
  column_index <- match(state_keys, column_keys)
  if (anyNA(column_index)) {
    stop(
      "Figure 2G: faltan estados en AverageExpression(): ",
      paste(figure2g_states[is.na(column_index)], collapse = ", ")
    )
  }
  avg_all <- avg_all[, column_index, drop = FALSE]
  colnames(avg_all) <- figure2g_states
  allstate_genes <- unique(c(cytokine_present, ecm_present))
  cat("Figure 2G: DE global para cytokines + ECM / collagen...\n")
  markers_allstates <- FindAllMarkers(
object = fibro_only,
assay = "RNA",
features = allstate_genes,
test.use = "wilcox",
only.pos = FALSE,
min.pct = 0,
logfc.threshold = 0,
return.thresh = 1.01,
verbose = FALSE
  )
  if (
    !"avg_log2FC" %in% colnames(markers_allstates) &&
      "avg_logFC" %in% colnames(markers_allstates)
  ) {
    markers_allstates <- markers_allstates %>%
      rename(avg_log2FC = avg_logFC)
  }
  if (nrow(markers_allstates) == 0) {
    stop( "Figure 2G: FindAllMarkers no devolvió resultados ", "para cytokines + ECM / collagen." )
  }
  cat( "Figure 2G: DE metabólico de Fibroblasts Activated 2 ", "versus los otros estados...\n", sep = "" )
  markers_metabolic <- FindMarkers(
object = fibro_only,
ident.1 = "Fibroblasts Activated 2",
ident.2 = c( "Fibroblasts Universal", "Fibroblasts Activated 1" ),
assay = "RNA",
features = metabolic_present,
test.use = "wilcox",
only.pos = FALSE,
min.pct = 0,
logfc.threshold = 0,
verbose = FALSE
  )
  markers_metabolic <- as.data.frame(markers_metabolic)
  markers_metabolic$gene <- rownames(markers_metabolic)
  if (
    !"avg_log2FC" %in% colnames(markers_metabolic) &&
      "avg_logFC" %in% colnames(markers_metabolic)
  ) {
    markers_metabolic <- markers_metabolic %>%
      rename(avg_log2FC = avg_logFC)
  }
  if (nrow(markers_metabolic) == 0) {
    stop( "Figure 2G: FindMarkers no devolvió resultados ", "para los genes metabólicos." )
  }
  figure2g_cache <- list(
signature = figure2g_cache_signature,
avg_all = avg_all,
markers_allstates = markers_allstates,
markers_metabolic = markers_metabolic
  )
  saveRDS( figure2g_cache, figure2g_cache_file )
  cat( "Cache de Figure 2G guardada en:\n", figure2g_cache_file, "\n", sep = "" )
}
figure2g_avg_all <- figure2g_cache$avg_all
markers_allstates <- figure2g_cache$markers_allstates
markers_metabolic <- figure2g_cache$markers_metabolic
# ============================================================
# Variation table
# ============================================================
get_expression_variation <- function(avg_mat, genes) {
  x <- avg_mat[ intersect(genes, rownames(avg_mat)), , drop = FALSE ]
  data.frame(
gene = rownames(x),
mean_expression = rowMeans(x, na.rm = TRUE),
expression_sd = apply(x, 1, sd, na.rm = TRUE),
expression_range = apply( x, 1, function(v) max(v, na.rm = TRUE) - min(v, na.rm = TRUE) ),
stringsAsFactors = FALSE
  ) %>%
    filter(
      mean_expression > 0,
      is.finite(expression_sd),
      expression_sd > 1e-8,
      is.finite(expression_range),
      expression_range > 1e-8
    )
}
# ============================================================
# Select cytokines required by the manuscript text
# ============================================================
select_manuscript_cytokine_genes <- function(
  genes,
  markers,
  avg_mat,
  mandatory_genes,
  tnf_family_candidates,
  n_genes,
  filename
) {
  variation <- get_expression_variation(avg_mat, genes)

  missing_from_variation <- setdiff(
    mandatory_genes,
    variation$gene
  )

  if (length(missing_from_variation) > 0) {
    stop(
      "Figure 2G cytokines: estos genes del manuscrito no tienen ",
      "expresión/variación suficiente para calcular el heatmap: ",
      paste(missing_from_variation, collapse = ", ")
    )
  }

  marker_table <- markers %>%
    filter(
      gene %in% genes,
      cluster %in% figure2g_states,
      !is.na(gene),
      gene != "",
      !is.na(avg_log2FC)
    ) %>%
    inner_join(
      variation,
      by = "gene"
    )

  n_tnf_needed <- n_genes - length(mandatory_genes)

  if (n_tnf_needed != 2L) {
    stop(
      "Figure 2G cytokines: se esperaban exactamente 2 genes TNF/TNFSF ",
      "adicionales, pero n_genes - mandatory_genes = ",
      n_tnf_needed,
      "."
    )
  }

  activated1_tnf <- marker_table %>%
    filter(
      cluster == "Fibroblasts Activated 1",
      gene %in% tnf_family_candidates,
      avg_log2FC > 0
    ) %>%
    distinct(
      gene,
      .keep_all = TRUE
    ) %>%
    arrange(
      desc(avg_log2FC),
      desc(expression_range),
      p_val_adj,
      gene
    )

  if (nrow(activated1_tnf) < n_tnf_needed) {
    stop(
      "Figure 2G cytokines: solo hay ",
      nrow(activated1_tnf),
      " genes TNF/TNFSF positivos en Fibroblasts Activated 1; ",
      "se necesitan ",
      n_tnf_needed,
      ". Candidatos evaluados: ",
      paste(
        tnf_family_candidates,
        collapse = ", "
      )
    )
  }

  activated1_tnf_selected <- activated1_tnf$gene[
    seq_len(n_tnf_needed)
  ]

  selected <- c(
    mandatory_genes,
    activated1_tnf_selected
  )

  if (
    length(selected) != n_genes ||
      anyDuplicated(selected)
  ) {
    stop(
      "Figure 2G cytokines: la selección final contiene ",
      length(selected),
      " genes únicos en lugar de ",
      n_genes,
      "."
    )
  }

  selected_table <- data.frame(
    Gene = selected,
    Manuscript_block = c(
      unname(
        manuscript_cytokine_groups[mandatory_genes]
      ),
      rep(
        "Activated fibroblast 1 - TNF/TNFSF family",
        length(activated1_tnf_selected)
      )
    ),
    Selection_reason = c(
      rep(
        "Explicitly named in manuscript text",
        length(mandatory_genes)
      ),
      rep(
        "Top positive Activated 1 TNF/TNFSF marker",
        length(activated1_tnf_selected)
      )
    ),
    stringsAsFactors = FALSE
  )

  selected_marker_summary <- marker_table %>%
    filter(gene %in% selected) %>%
    group_by(gene) %>%
    slice_max(
      order_by = abs(avg_log2FC),
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup() %>%
    select(
      gene,
      cluster,
      avg_log2FC,
      p_val_adj,
      mean_expression,
      expression_sd,
      expression_range
    )

  selected_table <- selected_table %>%
    left_join(
      selected_marker_summary,
      by = c("Gene" = "gene")
    )

  write.table(
    selected_table,
    file.path(
      table_dir,
      paste0(
        filename,
        "_manuscript_genes.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat(
    "\nCytokines / chemokines / growth factors seleccionadas para Figure 2G:\n",
    "  Genes explícitos del texto: ",
    length(mandatory_genes),
    "\n",
    "  Genes TNF/TNFSF Activated 1 añadidos: ",
    paste(
      activated1_tnf_selected,
      collapse = ", "
    ),
    "\n",
    "  Total: ",
    length(selected),
    "\n",
    sep = ""
  )

  selected
}

# ============================================================
# Select cytokines with representation from all fibroblast states
#
# Five positive one-vs-rest markers are reserved for each fibroblast state.
# The remaining positions are filled with the strongest
# remaining genes from the corresponding candidate set.
# ============================================================
select_balanced_state_genes <- function(
  genes,
  markers,
  avg_mat,
  n_genes,
  n_per_state = 5L,
  filename,
  panel_label
) {
  variation <- get_expression_variation(avg_mat, genes)
  marker_table <- markers %>%
    filter( gene %in% genes, cluster %in% figure2g_states, !is.na(gene), gene != "", !is.na(avg_log2FC) ) %>%
    inner_join(variation, by = "gene")
  positive_markers <- marker_table %>%
    filter(avg_log2FC > 0) %>%
    arrange( cluster, desc(avg_log2FC), desc(expression_range), p_val_adj, gene )
  selected <- character(0)
  selected_info <- list()
  for (state in figure2g_states) {
    state_candidates <- positive_markers %>%
      filter( cluster == state, !gene %in% selected ) %>%
      distinct(gene, .keep_all = TRUE) %>%
      slice_head(n = n_per_state)
    if (nrow(state_candidates) < n_per_state) {
      stop(
        "Figure 2G ", panel_label, ": solo hay ",
        nrow(state_candidates),
        " genes positivos distintivos disponibles para ",
        state,
        "; se necesitan ",
        n_per_state,
        "."
      )
    }
    selected <- c(selected, state_candidates$gene)
    selected_info[[state]] <- state_candidates %>%
      transmute(
        gene,
        Representative_state = as.character(cluster),
        Selection_reason = "State-specific",
        avg_log2FC,
        p_val_adj,
        mean_expression,
        expression_sd,
        expression_range
      )
  }
  global_ranking <- marker_table %>%
    group_by(gene) %>%
    summarise(
      max_abs_log2FC = max(abs(avg_log2FC), na.rm = TRUE),
      min_p_val_adj = min(p_val_adj, na.rm = TRUE),
      mean_expression = first(mean_expression),
      expression_sd = first(expression_sd),
      expression_range = first(expression_range),
      .groups = "drop"
    ) %>%
    arrange( desc(max_abs_log2FC), desc(expression_range), min_p_val_adj, gene )
  n_remaining <- n_genes - length(selected)
  if (n_remaining < 0) {
    stop( "Figure 2G ", panel_label, ": n_per_state genera más genes ", "de los permitidos por n_genes." )
  }
  if (n_remaining > 0) {
    fill_genes <- global_ranking %>%
      filter(!gene %in% selected) %>%
      slice_head(n = n_remaining) %>%
      pull(gene)
    if (length(fill_genes) < n_remaining) {
      stop( "Figure 2G ", panel_label, ": no hay suficientes genes para completar ", n_genes, " genes únicos." )
    }
    selected <- c(selected, fill_genes)
  }
  if (length(selected) != n_genes || anyDuplicated(selected)) {
    stop( "Figure 2G ", panel_label, ": la selección final no contiene exactamente ", n_genes, " genes únicos." )
  }
  selected_state_table <- bind_rows(selected_info)
  ranking_output <- global_ranking %>%
    mutate(
      Selected = gene %in% selected,
      Selection_reason = case_when(
        gene %in% selected_state_table$gene ~ "State-specific",
        gene %in% selected ~ "Global fill",
        TRUE ~ "Not selected"
      ),
      Representative_state = selected_state_table$Representative_state[ match(gene, selected_state_table$gene) ]
    )
  write.table(
    ranking_output,
    file.path( table_dir, paste0(filename, "_balanced_ranking.tsv") ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  cat(
    "\n", panel_label, " seleccionados:\n",
    "  ", n_per_state, " distintivas de cada estado (",
    n_per_state * length(figure2g_states), " genes)\n",
    "  ", n_remaining, " genes adicionales por DE global\n",
    "  Total: ", length(selected), "\n",
    sep = ""
  )
  selected
}
# ============================================================
# Select top genes from all-state DE
# ============================================================
select_allstate_genes <- function(
genes,
markers,
avg_mat,
n_genes,
filename
) {
  variation <- get_expression_variation(avg_mat, genes)
  ranking <- markers %>%
    filter( gene %in% genes, !is.na(gene), gene != "", !is.na(avg_log2FC) ) %>%
    group_by(gene) %>%
    summarise(
max_abs_log2FC = max(abs(avg_log2FC), na.rm = TRUE),
min_p_val_adj = min(p_val_adj, na.rm = TRUE),
.groups = "drop"
    ) %>%
    inner_join(variation, by = "gene") %>%
    arrange( desc(max_abs_log2FC), desc(expression_range), min_p_val_adj, gene )
  if (nrow(ranking) < n_genes) {
    stop( filename, ": solo quedan ", nrow(ranking), " genes utilizables; se necesitan ", n_genes, "." )
  }
  ranking$Selected <- FALSE
  ranking$Selected[seq_len(n_genes)] <- TRUE
  write.table(
    ranking,
    file.path( table_dir, paste0(filename, "_differential_ranking.tsv") ),
sep = "\t",
quote = FALSE,
row.names = FALSE
  )
  ranking$gene[seq_len(n_genes)]
}
# ============================================================
# Select metabolic genes specifically altered in Fib2
# ============================================================
select_fib2_metabolic_genes <- function(
genes,
markers,
avg_mat,
n_genes,
filename
) {
  variation <- get_expression_variation(avg_mat, genes)
  ranking <- markers %>%
    filter( gene %in% genes, !is.na(gene), gene != "", !is.na(avg_log2FC) ) %>%
    mutate( abs_log2FC_Fib2 = abs(avg_log2FC), Direction_in_Fib2 = ifelse( avg_log2FC >= 0, "Up", "Down" ) ) %>%
    inner_join(variation, by = "gene") %>%
    arrange( desc(abs_log2FC_Fib2), desc(expression_range), p_val_adj, gene )
  if (nrow(ranking) < n_genes) {
    stop( filename, ": solo quedan ", nrow(ranking), " genes metabólicos utilizables; se necesitan ", n_genes, "." )
  }
  ranking$Selected <- FALSE
  ranking$Selected[seq_len(n_genes)] <- TRUE
  write.table(
    ranking,
    file.path( table_dir, paste0(filename, "_Fib2_vs_rest_ranking.tsv") ),
sep = "\t",
quote = FALSE,
row.names = FALSE
  )
  selected <- ranking$gene[seq_len(n_genes)]
  cat( "\nMetabolism - top Fib2 genes:\n", paste(selected, collapse = ", "), "\n", sep = "" )
  selected
}
# ============================================================
# Select exactly the same number of genes for all three panels
# ============================================================
cytokine_genes_selected <- select_manuscript_cytokine_genes(
  genes = cytokine_present,
  markers = markers_allstates,
  avg_mat = figure2g_avg_all,
  mandatory_genes = manuscript_cytokine_genes,
  tnf_family_candidates = activated1_tnf_family_candidates,
  n_genes = figure2g_n_genes,
  filename = "Figure2G_cytokines_chemokines"
)
ecm_genes_selected <- select_balanced_state_genes(
  genes = ecm_present,
  markers = markers_allstates,
  avg_mat = figure2g_avg_all,
  n_genes = figure2g_n_genes,
  n_per_state = 5L,
  filename = "Figure2G_ECM_collagen",
  panel_label = "ECM / collagen"
)
metabolic_genes_selected <- select_fib2_metabolic_genes(
genes = metabolic_present,
markers = markers_metabolic,
avg_mat = figure2g_avg_all,
n_genes = figure2g_n_genes,
filename = "Figure2G_metabolism"
)
# ============================================================
# Common Z-score scale and legend
# ============================================================
figure2g_zscore_colors <- circlize::colorRamp2( c(-2.5, 0, 2.5), c("navy", "white", "firebrick3") )
figure2g_legend_breaks <- c(-2.5, -1.25, 0, 1.25, 2.5)
figure2g_legend <- ComplexHeatmap::Legend(
title = "Z-score",
title_position = "topcenter",
col_fun = figure2g_zscore_colors,
at = figure2g_legend_breaks
)
figure2g_padding <- grid::unit(c(7, 8, 10, 18), "mm")
# ============================================================
# Build one vertical heatmap
# ============================================================
make_figure2g_heatmap <- function(
  avg_mat,
  genes,
  title,
  heatmap_name,
  filename
) {
  # Average expression matrix: genes x fibroblast states
  panel_avg <- avg_mat[
    genes,
    figure2g_states,
    drop = FALSE
  ]

  # Gene-wise Z-score across the three fibroblast states.
  # This is calculated before transposing, so scaling remains
  # identical to the previous vertical version.
  gene_sd <- apply(
    panel_avg,
    1,
    sd,
    na.rm = TRUE
  )

  bad_genes <- rownames(panel_avg)[
    !is.finite(gene_sd) |
      gene_sd <= 1e-8
  ]

  if (length(bad_genes) > 0) {
    stop(
      "Figure 2G: genes sin variación entre estados: ",
      paste(
        bad_genes,
        collapse = ", "
      )
    )
  }

  gene_mean <- rowMeans(
    panel_avg,
    na.rm = TRUE
  )

  z_gene_state <- sweep(
    panel_avg,
    1,
    gene_mean,
    "-"
  )

  z_gene_state <- sweep(
    z_gene_state,
    1,
    gene_sd,
    "/"
  )

  if (any(!is.finite(z_gene_state))) {
    stop(
      "Figure 2G: se generaron Z-scores no finitos en ",
      title,
      "."
    )
  }

  z_gene_state[z_gene_state > 2.5] <- 2.5
  z_gene_state[z_gene_state < -2.5] <- -2.5

  # Cluster genes exactly as before, but display them as columns.
  gene_hc <- hclust(
    dist(z_gene_state),
    method = "complete"
  )

  # Horizontal orientation:
  # rows    = fibroblast states
  # columns = genes
  z_horizontal <- t(z_gene_state)

  gene_name_gp <- grid::gpar(
    fontsize = 10,
    fontface = "italic"
  )

  Figure2G <- ComplexHeatmap::Heatmap(
    z_horizontal,
    name = heatmap_name,
    col = figure2g_zscore_colors,
    show_heatmap_legend = FALSE,

    # Fibroblast states remain in the biological order already
    # defined in figure2g_states.
    cluster_rows = FALSE,
    show_row_dend = FALSE,

    # Preserve gene clustering, now along the horizontal axis.
    cluster_columns = gene_hc,
    show_column_dend = TRUE,

    show_row_names = TRUE,
    show_column_names = TRUE,

    row_names_side = "left",
    row_names_gp = grid::gpar(
      fontsize = 11
    ),

    column_names_side = "bottom",
    column_names_rot = 45,
    column_names_gp = gene_name_gp,

    column_title = title,
    column_title_side = "top",
    column_title_gp = grid::gpar(
      fontsize = 13,
      fontface = "bold"
    ),

    border = FALSE,
    use_raster = FALSE
  )

  # Individual horizontal heatmaps.
  pdf(
    file = file.path(
      figure_dir,
      paste0(
        filename,
        ".pdf"
      )
    ),
    width = 10.5,
    height = 3.6,
    useDingbats = FALSE
  )

  ComplexHeatmap::draw(
    Figure2G,
    heatmap_legend_list = list(
      figure2g_legend
    ),
    heatmap_legend_side = "right",
    padding = grid::unit(
      c(
        8,
        8,
        18,
        14
      ),
      "mm"
    )
  )

  dev.off()

  png(
    filename = file.path(
      figure_dir,
      paste0(
        filename,
        ".png"
      )
    ),
    width = 10.5,
    height = 3.6,
    units = "in",
    res = 300
  )

  ComplexHeatmap::draw(
    Figure2G,
    heatmap_legend_list = list(
      figure2g_legend
    ),
    heatmap_legend_side = "right",
    padding = grid::unit(
      c(
        8,
        8,
        18,
        14
      ),
      "mm"
    )
  )

  dev.off()

  # Keep the exported table in the intuitive format
  # genes x fibroblast states, even though the plot is transposed.
  write.table(
    data.frame(
      Gene = rownames(z_gene_state),
      z_gene_state,
      check.names = FALSE
    ),
    file.path(
      table_dir,
      paste0(
        filename,
        "_zscore.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat(
    title,
    ": ",
    nrow(z_gene_state),
    " genes\n",
    sep = ""
  )

  Figure2G
}

# ============================================================
# Generate the three HORIZONTAL Figure 2G heatmaps
# ============================================================

Figure2G_cytokines <- make_figure2g_heatmap(
  avg_mat = figure2g_avg_all,
  genes = cytokine_genes_selected,
  title = "Cytokines / chemokines / growth factors",
  heatmap_name = "Figure2G_cytokines",
  filename = "Figure2G_cytokines_chemokines"
)

Figure2G_ECM_collagen <- make_figure2g_heatmap(
  avg_mat = figure2g_avg_all,
  genes = ecm_genes_selected,
  title = "ECM / collagen",
  heatmap_name = "Figure2G_ECM_collagen",
  filename = "Figure2G_ECM_collagen"
)

Figure2G_metabolic <- make_figure2g_heatmap(
  avg_mat = figure2g_avg_all,
  genes = metabolic_genes_selected,
  title = "Metabolism",
  heatmap_name = "Figure2G_metabolism",
  filename = "Figure2G_metabolism"
)

# ============================================================
# Combined Figure 2G battery
#
# The individual heatmaps are horizontal, so the combined figure
# is arranged in THREE ROWS to keep the gene labels readable.
#
# Each heatmap is drawn independently inside its own viewport so
# that the three different gene sets retain their own clustering.
# A single common Z-score legend is drawn on the right.
# ============================================================

draw_figure2g_combined <- function(
  filename,
  device = c(
    "pdf",
    "png"
  )
) {
  device <- match.arg(device)

  if (device == "pdf") {
    pdf(
      file = file.path(
        figure_dir,
        filename
      ),
      width = 11.5,
      height = 10.5,
      useDingbats = FALSE
    )
  } else {
    png(
      filename = file.path(
        figure_dir,
        filename
      ),
      width = 11.5,
      height = 10.5,
      units = "in",
      res = 300
    )
  }

  grid::grid.newpage()

  outer_layout <- grid::grid.layout(
    nrow = 3,
    ncol = 2,
    heights = grid::unit(
      c(
        1,
        1,
        1
      ),
      "null"
    ),
    widths = grid::unit.c(
      grid::unit(
        1,
        "null"
      ),
      grid::unit(
        18,
        "mm"
      )
    )
  )

  grid::pushViewport(
    grid::viewport(
      layout = outer_layout
    )
  )

  heatmap_list <- list(
    Figure2G_cytokines,
    Figure2G_ECM_collagen,
    Figure2G_metabolic
  )

  for (i in seq_along(heatmap_list)) {
    grid::pushViewport(
      grid::viewport(
        layout.pos.row = i,
        layout.pos.col = 1
      )
    )

    ComplexHeatmap::draw(
      heatmap_list[[i]],
      newpage = FALSE,
      padding = grid::unit(
        c(
          5,
          5,
          14,
          9
        ),
        "mm"
      )
    )

    grid::upViewport()
  }

  # One common legend for all three heatmaps.
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 1:3,
      layout.pos.col = 2
    )
  )

  grid::grid.draw(
    ComplexHeatmap::packLegend(
      figure2g_legend,
      direction = "vertical"
    )
  )

  grid::upViewport()
  grid::upViewport()

  dev.off()
}

draw_figure2g_combined(
  filename = "Figure2G_combined_horizontal_3rows.pdf",
  device = "pdf"
)

draw_figure2g_combined(
  filename = "Figure2G_combined_horizontal_3rows.png",
  device = "png"
)

cat(
  "\nFigure 2G combinada horizontal guardada en:\n",
  file.path(
    figure_dir,
    "Figure2G_combined_horizontal_3rows.pdf"
  ),
  "\n",
  file.path(
    figure_dir,
    "Figure2G_combined_horizontal_3rows.png"
  ),
  "\n",
  sep = ""
)

capture.output( sessionInfo(), file = file.path(table_dir, "Figure2_sessionInfo.txt") )
cat("\nFigure 2 finalizada.\n")