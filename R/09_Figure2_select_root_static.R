# ============================================================
# 09_Figure2_select_root_static.R
# Static node map + manual root selection from the terminal
# Fibroblast-only Monocle trajectory
#
# First generate the node map:
#   Rscript R/09_Figure2_select_root_static.R
#
# Then save one or more root nodes:
#   Rscript R/09_Figure2_select_root_static.R Y_3
#   Rscript R/09_Figure2_select_root_static.R Y_3,Y_17
# ============================================================

suppressPackageStartupMessages({
  library(monocle3)
  library(ggplot2)
  library(igraph)
})

if (interactive()) {
  stop(
    "No ejecutes este helper con source().\n",
    "Sal de R y ejecútalo desde la terminal con:\n",
    "Rscript R/09_Figure2_select_root_static.R"
  )
}

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
table_dir <- file.path(rsdir, "tables", "Figure2")
figure_dir <- file.path(outdir, "Figure2")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

preorder_cds_file <- file.path(
  objects_dir,
  "Figure2_fibroblast_before_ordering.rds"
)

manual_root_nodes_file <- file.path(
  table_dir,
  "Figure2E_manual_root_pr_nodes_fibro_only.txt"
)

manual_selection_file <- file.path(
  table_dir,
  "Figure2E_manual_root_selection_fibro_only.tsv"
)

if (!file.exists(preorder_cds_file)) {
  stop(
    "No se encuentra:\n",
    preorder_cds_file,
    "\n\nEjecuta primero:\n",
    "Rscript R/09_Figure2.R"
  )
}

cat("\nCargando el CDS de fibroblastos...\n")
cds <- readRDS(preorder_cds_file)

if (is.null(principal_graph(cds)[["UMAP"]])) {
  stop("El CDS no contiene el principal graph de UMAP.")
}

expected_source <- "fibroblasts.Fibros.annotated.rds"
expected_analysis <- "fibroblast_only_seurat_umap"

if (!"Figure2_source_object" %in% colnames(colData(cds))) {
  stop(
    "El CDS no contiene información sobre el objeto de origen.\n",
    "Elimina el CDS y vuelve a ejecutar:\n",
    "Rscript R/09_Figure2.R"
  )
}

cds_source <- unique(
  as.character(
    colData(cds)$Figure2_source_object
  )
)

cds_source <- cds_source[
  !is.na(cds_source) & nzchar(cds_source)
]

if (
  length(cds_source) != 1L ||
    cds_source != expected_source
) {
  stop(
    "El CDS no procede del objeto de fibroblastos anotados.\n",
    "Objeto esperado: ",
    expected_source,
    "\nObjeto registrado: ",
    paste(cds_source, collapse = ", "),
    "\n\nBorra el CDS previo y vuelve a ejecutar 09_Figure2.R."
  )
}

if (!"Figure2_analysis" %in% colnames(colData(cds))) {
  stop(
    "El CDS no identifica el tipo de análisis.\n",
    "Borra el CDS previo y vuelve a ejecutar 09_Figure2.R."
  )
}

cds_analysis <- unique(
  as.character(
    colData(cds)$Figure2_analysis
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
    "El CDS no corresponde al análisis fibroblast_only_seurat_umap.\n",
    "Tipo registrado: ",
    paste(cds_analysis, collapse = ", "),
    "\n\nBorra el CDS previo y vuelve a ejecutar 09_Figure2.R."
  )
}

if (!"Clustering" %in% colnames(colData(cds))) {
  stop("El CDS no contiene la columna 'Clustering'.")
}

fibro_states <- c(
  "Fibroblasts Universal",
  "Fibroblasts Activated 1",
  "Fibroblasts Activated 2"
)

cds_states <- unique(
  as.character(
    colData(cds)$Clustering
  )
)

cds_states <- cds_states[
  !is.na(cds_states)
]

unexpected_states <- setdiff(
  cds_states,
  fibro_states
)

if (length(unexpected_states) > 0) {
  stop(
    "El CDS contiene poblaciones no permitidas: ",
    paste(unexpected_states, collapse = ", ")
  )
}

if (
  any(
    as.character(
      colData(cds)$Clustering
    ) == "Mesothelial"
  )
) {
  stop(
    "El mapa contiene células Mesothelial. ",
    "Borra el CDS y vuelve a ejecutar 09_Figure2.R."
  )
}

cat("\nPoblaciones presentes en el mapa de nodos:\n")
print(
  table(
    colData(cds)$Clustering,
    useNA = "ifany"
  )
)

trajectory_colors <- c(
  "Fibroblasts Universal" = "#FBCDC7",
  "Fibroblasts Activated 1" = "#F37B79",
  "Fibroblasts Activated 2" = "#6D1A1A"
)

graph_nodes <- igraph::V(
  principal_graph(cds)[["UMAP"]]
)$name

cat(
  "Nodos del principal graph: ",
  length(graph_nodes),
  "\n",
  sep = ""
)

# ============================================================
# Static node map containing fibroblasts only
# ============================================================

root_map <- plot_cells(
  cds,
  reduction_method = "UMAP",
  color_cells_by = "Clustering",
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "black",
  trajectory_graph_segment_size = 0.6,
  label_principal_points = TRUE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  label_roots = FALSE,
  graph_label_size = 2.5,
  cell_size = 0.7
) +
  scale_color_manual(
    values = trajectory_colors,
    breaks = fibro_states,
    drop = FALSE
  ) +
  labs(
    title = "Principal graph: select root node Y_*",
    x = "UMAP1",
    y = "UMAP2",
    color = NULL
  ) +
  theme_classic(base_size = 12)

root_map_pdf <- file.path(
  figure_dir,
  "Figure2E_root_nodes_fibroblast.pdf"
)

root_map_png <- file.path(
  figure_dir,
  "Figure2E_root_nodes_fibroblast.png"
)

ggsave(
  filename = root_map_pdf,
  plot = root_map,
  width = 8,
  height = 7,
  units = "in",
  useDingbats = FALSE
)

ggsave(
  filename = root_map_png,
  plot = root_map,
  width = 8,
  height = 7,
  units = "in",
  dpi = 300
)

cat(
  "\nMapa estático guardado en:\n",
  root_map_pdf,
  "\n",
  root_map_png,
  "\n",
  sep = ""
)

# ============================================================
# Manual selection supplied as command-line argument
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cat(
    "\nNo se ha seleccionado todavía ningún root.\n",
    "Abre el PNG, identifica el nodo o nodos Y_* y ejecuta, por ejemplo:\n\n",
    "  Rscript R/09_Figure2_select_root_static.R Y_17\n\n",
    "Para varios nodos:\n\n",
    "  Rscript R/09_Figure2_select_root_static.R Y_17,Y_42\n",
    sep = ""
  )

  quit(save = "no", status = 0)
}

root_input <- paste(args, collapse = ",")

selected_root_nodes <- trimws(
  unlist(
    strsplit(
      root_input,
      ",",
      fixed = TRUE
    )
  )
)

selected_root_nodes <- unique(
  selected_root_nodes[
    nzchar(selected_root_nodes)
  ]
)

if (length(selected_root_nodes) == 0) {
  stop("No se ha indicado ningún nodo raíz.")
}

invalid_nodes <- setdiff(
  selected_root_nodes,
  graph_nodes
)

if (length(invalid_nodes) > 0) {
  stop(
    "Estos nodos no existen en el principal graph: ",
    paste(invalid_nodes, collapse = ", ")
  )
}

cat(
  "\nNodo(s) seleccionado(s): ",
  paste(selected_root_nodes, collapse = ", "),
  "\nCalculando preview del pseudotime...\n",
  sep = ""
)

cds_preview <- order_cells(
  cds,
  reduction_method = "UMAP",
  root_pr_nodes = selected_root_nodes
)

pseudotime_values <- pseudotime(
  cds_preview,
  reduction_method = "UMAP"
)

writeLines(
  selected_root_nodes,
  manual_root_nodes_file
)

write.table(
  data.frame(
    Root_pr_node = selected_root_nodes,
    Selection_method = "manual_static_node_map_fibroblast_only",
    stringsAsFactors = FALSE
  ),
  manual_selection_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  cds_preview,
  file.path(
    objects_dir,
    "Figure2_fibroblast_manual_preview.rds"
  )
)

manual_preview <- plot_cells(
  cds_preview,
  reduction_method = "UMAP",
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "black",
  trajectory_graph_segment_size = 0.6,
  label_principal_points = FALSE,
  label_branch_points = FALSE,
  label_leaves = FALSE,
  label_roots = TRUE,
  graph_label_size = 3,
  cell_size = 0.7
) +
  scale_color_viridis_c(
    option = "plasma",
    name = "Pseudotime"
  ) +
  labs(
    title = paste0(
      "Manual root: ",
      paste(selected_root_nodes, collapse = ", ")
    ),
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme_classic(base_size = 12)

preview_pdf <- file.path(
  figure_dir,
  "Figure2E_root_preview_fibroblast.pdf"
)

preview_png <- file.path(
  figure_dir,
  "Figure2E_root_preview_fibroblast.png"
)

ggsave(
  filename = preview_pdf,
  plot = manual_preview,
  width = 7,
  height = 6,
  units = "in",
  useDingbats = FALSE
)

ggsave(
  filename = preview_png,
  plot = manual_preview,
  width = 7,
  height = 6,
  units = "in",
  dpi = 300
)

cat(
  "\nSelección guardada.\n",
  "Root node(s): ",
  paste(selected_root_nodes, collapse = ", "),
  "\nCélulas con pseudotime finito: ",
  sum(is.finite(pseudotime_values)),
  "\nCélulas con pseudotime infinito: ",
  sum(!is.finite(pseudotime_values)),
  "\nPreview:\n",
  preview_png,
  "\n\nPara generar la figura definitiva ejecuta:\n",
  "Rscript R/09_Figure2.R\n",
  sep = ""
)