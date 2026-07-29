# ============================================================
# 09_Figure2_select_root_static.R
# Static node map + manual root selection from the terminal
# No Shiny interface and no browser is launched
# Launch with default parameters "Rscript R/09_Figure2_select_root_static.R Y_3, Y_17" in terminal
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
  "Figure2_fibroblast_mesothelial_before_ordering.rds"
)

manual_root_nodes_file <- file.path(
  table_dir,
  "Figure2E_manual_root_pr_nodes.txt"
)

manual_selection_file <- file.path(
  table_dir,
  "Figure2E_manual_root_selection.tsv"
)

if (!file.exists(preorder_cds_file)) {
  stop(
    "No se encuentra:\n",
    preorder_cds_file,
    "\n\nEjecuta primero:\n",
    "Rscript R/09_Figure2.R"
  )
}

cat("\nCargando el CDS con fibroblastos y mesoteliales...\n")
cds <- readRDS(preorder_cds_file)

if (is.null(principal_graph(cds)[["UMAP"]])) {
  stop("El CDS no contiene el principal graph de UMAP.")
}

cat("\nPoblaciones presentes en el mapa de nodos:\n")
print(table(colData(cds)$Clustering, useNA = "ifany"))

if (sum(as.character(colData(cds)$Clustering) == "Mesothelial") == 0) {
  stop(
    "El mapa de selección no contiene células Mesothelial. ",
    "Borra el CDS previo y vuelve a ejecutar 09_Figure2.R."
  )
}

trajectory_colors <- c(
  "Fibroblasts Universal" = "#FBCDC7",
  "Fibroblasts Activated 1" = "#F37B79",
  "Fibroblasts Activated 2" = "#6D1A1A",
  "Mesothelial" = "#FED579"
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
# Static map containing fibroblasts + mesothelial cells
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
  scale_color_manual(values = trajectory_colors) +
  labs(
    title = "Principal graph: select root node Y_*",
    x = "UMAP1",
    y = "UMAP2",
    color = NULL
  ) +
  theme_classic(base_size = 12)

root_map_pdf <- file.path(
  figure_dir,
  "Figure2E_root_nodes_fibro_mesothelial.pdf"
)

root_map_png <- file.path(
  figure_dir,
  "Figure2E_root_nodes_fibro_mesothelial.png"
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
# Examples:
# Rscript R/09_Figure2_select_root_static.R Y_17
# Rscript R/09_Figure2_select_root_static.R Y_17,Y_42
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
  unlist(strsplit(root_input, ",", fixed = TRUE))
)
selected_root_nodes <- unique(
  selected_root_nodes[nzchar(selected_root_nodes)]
)

if (length(selected_root_nodes) == 0) {
  stop("No se ha indicado ningún nodo raíz.")
}

invalid_nodes <- setdiff(selected_root_nodes, graph_nodes)
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
    Selection_method = "manual_static_node_map",
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
    "Figure2_fibroblast_mesothelial_manual_preview.rds"
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
  "Figure2E_root_preview_fibro_mesothelial.pdf"
)

preview_png <- file.path(
  figure_dir,
  "Figure2E_root_preview_fibro_mesothelial.png"
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
