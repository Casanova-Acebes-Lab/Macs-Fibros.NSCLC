
# ============================================================
# 11_Figure4_CellChat.R
#
# Figure 4A-C: CellChat
#
# INPUT ÚNICO:
#   results/objects/fibroblasts.Not.Integrated.annotated.rds
#
# Este script:
#   1. crea las carpetas necesarias;
#   2. construye los cuatro objetos CellChat desde Seurat;
#   3. usa computeCommunProb(raw.use = TRUE), sin proyección PPI;
#   4. guarda los objetos CellChat terminados;
#   5. genera Figure 4A, 4B y 4C;
#   6. guarda las tablas correspondientes.
#
# No necesita barcodes externos, tablas previas, plots previos
# ni objetos CellChat preexistentes.
#
# IMPORTANTE:
# Ejecutar el archivo completo con Rscript. No pegar por bloques
# ni enviar selecciones parciales a radian.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(future)
  library(circlize)
})

# ============================================================
# Parámetros
# ============================================================

seed <- 1234L
min_cells <- 1L
bubble_pvalue_threshold <- 0.1

# Figure 4C:
# subsetCommunication() applies p < 0.05 by default.
# Set to 1 to recover every inferred ligand-receptor pair and let
# the chord ranking, rather than significance filtering, determine order.
chord_pvalue_threshold <- 1

# FALSE por defecto:
#   reconstruye todos los objetos exclusivamente desde Seurat.
#
# Después de completar una ejecución correcta puede cambiarse a TRUE
# para regenerar únicamente plots y tablas desde los RDS creados
# por este mismo script.
reuse_completed_cellchat <- FALSE

# Método rápido y estable usado en la versión que funcionó:
# no se aplica smoothData() ni projectData().
# CellChat calcula las probabilidades desde RNA/data:
# computeCommunProb(raw.use = TRUE).
use_ppi_projection <- FALSE

set.seed(seed)

# Evitar workers de background y bloqueos al ejecutar en servidor.
future::plan(future::sequential)
options(future.globals.maxSize = 8 * 1024^3)


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

values <- setNames(
  config$value,
  config$key
)

required_config_keys <- c(
  "rsdir",
  "outdir"
)

missing_config_keys <- setdiff(
  required_config_keys,
  names(values)
)

if (length(missing_config_keys) > 0L) {
  stop(
    "Faltan claves en config.tsv: ",
    paste(missing_config_keys, collapse = ", ")
  )
}

rsdir <- unname(values[["rsdir"]])
outdir <- unname(values[["outdir"]])

objects_dir <- file.path(
  rsdir,
  "objects"
)

figure_dir <- file.path(
  outdir,
  "Figure4",
  "CellChat"
)

table_dir <- file.path(
  rsdir,
  "tables",
  "Figure4",
  "CellChat"
)

cellchat_dir <- file.path(
  objects_dir,
  "Figure4_CellChat"
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

dir.create(
  cellchat_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

log_file <- file.path(
  table_dir,
  "Figure4_CellChat_run.log"
)

if (file.exists(log_file)) {
  file.remove(log_file)
}

log_message <- function(...) {

  current_text <- paste0(...)

  current_line <- paste0(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    current_text
  )

  message(current_line)

  cat(
    current_line,
    "\n",
    file = log_file,
    append = TRUE
  )

  flush.console()
}

save_tsv <- function(
  x,
  filename
) {

  write.table(
    x,
    file = file.path(
      table_dir,
      filename
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

# ============================================================
# Poblaciones utilizadas
# ============================================================

receiver_order <- c(
  "Fibroblasts Activated 1",
  "Fibroblasts Activated 2",
  "Fibroblasts Universal",
  "Mesothelial",
  "MoMacs",
  "Monocytes Ly6cHi",
  "Monocytes Ly6cLo",
  "mregDC",
  "pDC",
  "DC1",
  "DC2"
)

states_TRM <- c(
  "TRMs",
  receiver_order
)

states_moTRM <- c(
  "mo-TRMs",
  receiver_order
)

state_order_TRM <- c(
  "TRMs",
  "DC2",
  "DC1",
  "pDC",
  "mregDC",
  "Monocytes Ly6cLo",
  "Monocytes Ly6cHi",
  "MoMacs",
  "Mesothelial",
  "Fibroblasts Universal",
  "Fibroblasts Activated 2",
  "Fibroblasts Activated 1"
)

state_order_moTRM <- c(
  "mo-TRMs",
  setdiff(
    state_order_TRM,
    "TRMs"
  )
)

figure4_colors <- c(
  "TRMs" = "#2D6A4F",
  "mo-TRMs" = "#66CDAA",
  "Mesothelial" = "#FED579",
  "Fibroblasts Universal" = "#FBCDC7",
  "Fibroblasts Activated 1" = "#F37B79",
  "Fibroblasts Activated 2" = "#6D1A1A",
  "Monocytes Ly6cLo" = "#FF4D4D",
  "Monocytes Ly6cHi" = "#D21F1F",
  "MoMacs" = "#A8D5BA",
  "mregDC" = "#FFB84D",
  "DC1" = "#7F60B6",
  "DC2" = "#D0A0D5",
  "pDC" = "#3E8BC4"
)

# ============================================================
# Leer el único input
# ============================================================

seurat_file <- file.path(
  objects_dir,
  "fibroblasts.Not.Integrated.annotated.rds"
)

if (!file.exists(seurat_file)) {
  stop(
    "No se encuentra el objeto Seurat: ",
    seurat_file
  )
}

log_message(
  "Leyendo objeto Seurat: ",
  seurat_file
)

data <- readRDS(
  seurat_file
)

required_metadata <- c(
  "Clustering",
  "orig.ident"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(data[[]])
)

if (length(missing_metadata) > 0L) {
  stop(
    "Faltan columnas de metadatos en el objeto Seurat: ",
    paste(missing_metadata, collapse = ", ")
  )
}

DefaultAssay(data) <- "RNA"

data <- JoinLayers(
  data,
  assay = "RNA"
)

# Confirmar que existe la capa normalizada.
rna_data <- GetAssayData(
  data,
  assay = "RNA",
  layer = "data"
)

if (nrow(rna_data) == 0L || ncol(rna_data) == 0L) {
  stop(
    "La capa RNA/data está vacía. CellChat necesita expresión ",
    "normalizada no negativa."
  )
}

rm(rna_data)
gc()

# ============================================================
# Homogeneizar anotaciones
# ============================================================

cell_state <- trimws(
  as.character(
    data$Clustering
  )
)

annotation_aliases <- c(
  "emAMs" = "TRMs",
  "em-AMs" = "TRMs",
  "moAMs" = "mo-TRMs",
  "mo-AMs" = "mo-TRMs",
  "Fibroblasts_Universal" = "Fibroblasts Universal",
  "Fibroblasts_Activated 1" = "Fibroblasts Activated 1",
  "Fibroblasts_Activated 2" = "Fibroblasts Activated 2"
)

alias_index <- match(
  cell_state,
  names(annotation_aliases)
)

cells_to_rename <- !is.na(
  alias_index
)

cell_state[
  cells_to_rename
] <- unname(
  annotation_aliases[
    alias_index[
      cells_to_rename
    ]
  ]
)

data$Figure4_cell_state <- cell_state

allowed_conditions <- c(
  "WT",
  "CD169DTR",
  "CCR2KO"
)

all_figure4_states <- unique(
  c(
    states_TRM,
    states_moTRM
  )
)

cells_to_keep <- colnames(data)[
  as.character(data$orig.ident) %in%
    allowed_conditions &
    as.character(data$Figure4_cell_state) %in%
    all_figure4_states
]

if (length(cells_to_keep) == 0L) {
  stop(
    "No se encontraron células compatibles con las condiciones ",
    "y anotaciones definidas para Figure 4."
  )
}

data <- subset(
  data,
  cells = cells_to_keep
)

cell_count_table <- as.data.frame(
  table(
    Condition = data$orig.ident,
    Cell_state = data$Figure4_cell_state
  ),
  stringsAsFactors = FALSE
)

colnames(cell_count_table)[3] <- "Cells"

save_tsv(
  cell_count_table,
  "Figure4_CellChat_cell_counts.tsv"
)

log_message(
  "Células seleccionadas para CellChat: ",
  ncol(data)
)

print(
  table(
    data$orig.ident,
    data$Figure4_cell_state
  )
)

# ============================================================
# Funciones de preparación
# ============================================================

prepare_seurat_subset <- function(
  object,
  conditions,
  allowed_states,
  state_order,
  object_name
) {

  metadata <- object[[]]

  selected_cells <- rownames(metadata)[
    as.character(metadata$orig.ident) %in%
      conditions &
    as.character(metadata$Figure4_cell_state) %in%
      allowed_states
  ]

  if (length(selected_cells) == 0L) {
    stop(
      object_name,
      ": no se encontraron células."
    )
  }

  current_object <- subset(
    object,
    cells = selected_cells
  )

  present_states <- state_order[
    state_order %in%
      unique(
        as.character(
          current_object$Figure4_cell_state
        )
      )
  ]

  current_object$Figure4_cell_state <- factor(
    as.character(
      current_object$Figure4_cell_state
    ),
    levels = present_states
  )

  current_counts <- table(
    current_object$Figure4_cell_state
  )

  log_message(
    object_name,
    ": ",
    ncol(current_object),
    " células y ",
    length(current_counts),
    " poblaciones."
  )

  print(current_counts)

  low_count_states <- names(
    current_counts[
      current_counts < min_cells
    ]
  )

  if (length(low_count_states) > 0L) {
    warning(
      object_name,
      ": poblaciones con menos de ",
      min_cells,
      " células: ",
      paste(low_count_states, collapse = ", ")
    )
  }

  current_object
}

group_index <- function(
  cellchat,
  group_name
) {

  current_index <- match(
    group_name,
    levels(cellchat@idents)
  )

  if (is.na(current_index)) {
    stop(
      "La población '",
      group_name,
      "' no está presente en el objeto CellChat. Poblaciones: ",
      paste(
        levels(cellchat@idents),
        collapse = ", "
      )
    )
  }

  current_index
}

present_receiver_indices <- function(
  cellchat
) {

  present_receivers <- receiver_order[
    receiver_order %in%
      levels(cellchat@idents)
  ]

  missing_receivers <- setdiff(
    receiver_order,
    present_receivers
  )

  if (length(missing_receivers) > 0L) {
    log_message(
      "Poblaciones receptoras ausentes en este objeto: ",
      paste(missing_receivers, collapse = ", ")
    )
  }

  match(
    present_receivers,
    levels(cellchat@idents)
  )
}

# ============================================================
# Construcción autónoma de CellChat
# ============================================================

build_cellchat <- function(
  seurat_object,
  object_name,
  output_file
) {

  if (
    reuse_completed_cellchat &&
    file.exists(output_file)
  ) {

    log_message(
      object_name,
      ": cargando objeto CellChat terminado: ",
      output_file
    )

    return(
      readRDS(output_file)
    )
  }

  log_message(
    object_name,
    ": iniciando CellChat desde el objeto Seurat."
  )

  seurat_object <- JoinLayers(
    seurat_object,
    assay = "RNA"
  )

  expression_data <- GetAssayData(
    seurat_object,
    assay = "RNA",
    layer = "data"
  )

  metadata <- seurat_object[[]][
    colnames(expression_data),
    ,
    drop = FALSE
  ]

  metadata$Figure4_cell_state <- droplevels(
    factor(
      metadata$Figure4_cell_state
    )
  )

  # CellChat v2 espera una columna `samples`.
  # Aquí cada objeto se analiza como un dataset independiente.
  metadata$samples <- as.character(
    metadata$orig.ident
  )

  set.seed(seed)

  log_message(
    object_name,
    ": createCellChat()."
  )

  cellchat <- createCellChat(
    object = expression_data,
    meta = metadata,
    group.by = "Figure4_cell_state"
  )

  log_message(
    object_name,
    ": base de datos mouse Secreted Signaling."
  )

  cellchat@DB <- subsetDB(
    CellChatDB.mouse,
    search = "Secreted Signaling"
  )

  log_message(
    object_name,
    ": subsetData()."
  )

  cellchat <- subsetData(
    cellchat
  )

  log_message(
    object_name,
    ": identifyOverExpressedGenes()."
  )

  cellchat <- identifyOverExpressedGenes(
    cellchat
  )

  log_message(
    object_name,
    ": identifyOverExpressedInteractions()."
  )

  cellchat <- identifyOverExpressedInteractions(
    cellchat
  )

  if (use_ppi_projection) {

    log_message(
      object_name,
      ": smoothData(adj = PPI.mouse)."
    )

    # En CellChat v2, smoothData() sustituye a projectData().
    cellchat <- CellChat::smoothData(
      cellchat,
      adj = PPI.mouse
    )
  }

  set.seed(seed)

  log_message(
    object_name,
    ": computeCommunProb()."
  )

  cellchat <- computeCommunProb(
    cellchat,
    raw.use = !use_ppi_projection
  )

  log_message(
    object_name,
    ": filterCommunication(min.cells = ",
    min_cells,
    ")."
  )

  cellchat <- filterCommunication(
    cellchat,
    min.cells = min_cells
  )

  log_message(
    object_name,
    ": computeCommunProbPathway()."
  )

  cellchat <- computeCommunProbPathway(
    cellchat
  )

  log_message(
    object_name,
    ": aggregateNet()."
  )

  cellchat <- aggregateNet(
    cellchat
  )

  saveRDS(
    cellchat,
    output_file
  )

  log_message(
    object_name,
    ": objeto guardado en ",
    output_file
  )

  cellchat
}

# ============================================================
# Preparar los cuatro datasets
# ============================================================

data_WT <- prepare_seurat_subset(
  object = data,
  conditions = "WT",
  allowed_states = states_TRM,
  state_order = state_order_TRM,
  object_name = "WT"
)

data_CD169DTR <- prepare_seurat_subset(
  object = data,
  conditions = "CD169DTR",
  allowed_states = states_moTRM,
  state_order = state_order_moTRM,
  object_name = "CD169DTR"
)

data_CCR2KO <- prepare_seurat_subset(
  object = data,
  conditions = "CCR2KO",
  allowed_states = states_TRM,
  state_order = state_order_TRM,
  object_name = "CCR2KO"
)

# CD169+ TRMs:
# WT + CCR2KO, excluyendo CD169DTR por la depleción de TRMs.
data_WT_CCR2KO <- prepare_seurat_subset(
  object = data,
  conditions = c(
    "WT",
    "CCR2KO"
  ),
  allowed_states = states_TRM,
  state_order = state_order_TRM,
  object_name = "WT + CCR2KO"
)

# ============================================================
# Construir CellChat
# ============================================================

cellchat_WT <- build_cellchat(
  seurat_object = data_WT,
  object_name = "WT",
  output_file = file.path(
    cellchat_dir,
    "Figure4_CellChat_WT.rds"
  )
)

cellchat_CD169DTR <- build_cellchat(
  seurat_object = data_CD169DTR,
  object_name = "CD169DTR",
  output_file = file.path(
    cellchat_dir,
    "Figure4_CellChat_CD169DTR.rds"
  )
)

cellchat_CCR2KO <- build_cellchat(
  seurat_object = data_CCR2KO,
  object_name = "CCR2KO",
  output_file = file.path(
    cellchat_dir,
    "Figure4_CellChat_CCR2KO.rds"
  )
)

cellchat_WT_CCR2KO <- build_cellchat(
  seurat_object = data_WT_CCR2KO,
  object_name = "WT + CCR2KO",
  output_file = file.path(
    cellchat_dir,
    "Figure4_CellChat_WT_CCR2KO.rds"
  )
)

cellchat_objects <- list(
  WT = cellchat_WT,
  CD169DTR = cellchat_CD169DTR,
  CCR2KO = cellchat_CCR2KO,
  WT_CCR2KO = cellchat_WT_CCR2KO
)

for (current_name in names(cellchat_objects)) {

  current_table <- subsetCommunication(
    cellchat_objects[[current_name]]
  )

  save_tsv(
    current_table,
    paste0(
      "Figure4_CellChat_",
      current_name,
      "_all_communications.tsv"
    )
  )
}

# Liberar objetos Seurat secundarios antes de dibujar.
rm(
  data_WT,
  data_CD169DTR,
  data_CCR2KO,
  data_WT_CCR2KO
)

gc()

# ============================================================
# Figure 4A
# Señales salientes de CD169+ TRMs y fibroblastos Universal
# en el objeto WT + CCR2KO
# ============================================================

single_source_matrix <- function(
  cellchat,
  source_name
) {

  current_matrix <- cellchat@net$weight

  if (!source_name %in% rownames(current_matrix)) {
    stop(
      "Figure 4A: no se encuentra la población emisora ",
      source_name
    )
  }

  output_matrix <- matrix(
    0,
    nrow = nrow(current_matrix),
    ncol = ncol(current_matrix),
    dimnames = dimnames(current_matrix)
  )

  output_matrix[
    source_name,
  ] <- current_matrix[
    source_name,
  ]

  output_matrix
}

matrix_TRM <- single_source_matrix(
  cellchat_WT_CCR2KO,
  "TRMs"
)

matrix_Universal <- single_source_matrix(
  cellchat_WT_CCR2KO,
  "Fibroblasts Universal"
)

edge_weight_max <- max(
  c(
    matrix_TRM,
    matrix_Universal
  ),
  na.rm = TRUE
)

if (
  !is.finite(edge_weight_max) ||
  edge_weight_max <= 0
) {
  stop(
    "Figure 4A: no se encontraron pesos positivos."
  )
}

figure4A_groups <- rownames(
  cellchat_WT_CCR2KO@net$weight
)

figure4A_group_size <- as.numeric(
  table(
    factor(
      cellchat_WT_CCR2KO@idents,
      levels = figure4A_groups
    )
  )
)

figure4A_colors <- figure4_colors[
  figure4A_groups
]

if (any(is.na(figure4A_colors))) {
  stop(
    "Figure 4A: faltan colores para: ",
    paste(
      figure4A_groups[
        is.na(figure4A_colors)
      ],
      collapse = ", "
    )
  )
}

draw_figure4A <- function(
  filename,
  device = c(
    "pdf",
    "png"
  )
) {

  device <- match.arg(
    device
  )

  if (device == "pdf") {

    pdf(
      filename,
      width = 12,
      height = 6.5,
      useDingbats = FALSE
    )

  } else {

    png(
      filename,
      width = 3600,
      height = 1950,
      res = 300
    )
  }

  old_par <- par(
    mfrow = c(1, 2),
    xpd = TRUE,
    mar = c(1, 1, 3, 1)
  )

  netVisual_circle(
    matrix_TRM,
    vertex.weight = figure4A_group_size,
    weight.scale = TRUE,
    edge.weight.max = edge_weight_max,
    label.edge = FALSE,
    color.use = figure4A_colors,
    title.name = "TRMs to cells"
  )

  netVisual_circle(
    matrix_Universal,
    vertex.weight = figure4A_group_size,
    weight.scale = TRUE,
    edge.weight.max = edge_weight_max,
    label.edge = FALSE,
    color.use = figure4A_colors,
    title.name = "Universal fibroblasts to cells"
  )

  par(old_par)
  dev.off()
}

draw_figure4A(
  filename = file.path(
    figure_dir,
    "Figure4A_outgoing_TRMs_Universal.pdf"
  ),
  device = "pdf"
)

draw_figure4A(
  filename = file.path(
    figure_dir,
    "Figure4A_outgoing_TRMs_Universal.png"
  ),
  device = "png"
)

save_tsv(
  subsetCommunication(
    cellchat_WT_CCR2KO,
    sources.use = "TRMs"
  ),
  "Figure4A_TRMs_outgoing.tsv"
)

save_tsv(
  subsetCommunication(
    cellchat_WT_CCR2KO,
    sources.use = "Fibroblasts Universal"
  ),
  "Figure4A_Universal_outgoing.tsv"
)

# ============================================================
# Figure 4B
# Bubble plots de TRMs/mo-TRMs hacia todos los tipos celulares
# presentes, incluidos Fibroblasts Activated 1/2 y DCs.
# ============================================================

make_bubble_plot <- function(
  cellchat,
  source_name,
  title_name
) {

  source_index <- group_index(
    cellchat,
    source_name
  )

  target_indices <- present_receiver_indices(
    cellchat
  )

  if (length(target_indices) == 0L) {
    stop(
      title_name,
      ": no hay poblaciones receptoras para el bubble plot."
    )
  }

  netVisual_bubble(
    cellchat,
    sources.use = source_index,
    targets.use = target_indices,
    remove.isolate = FALSE,
    thresh = bubble_pvalue_threshold,
    title.name = title_name
  ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      )
    )
}

p_WT <- make_bubble_plot(
  cellchat = cellchat_WT,
  source_name = "TRMs",
  title_name = "WT"
)

p_CD169DTR <- make_bubble_plot(
  cellchat = cellchat_CD169DTR,
  source_name = "mo-TRMs",
  title_name = "CD169DTR"
)

p_CCR2KO <- make_bubble_plot(
  cellchat = cellchat_CCR2KO,
  source_name = "TRMs",
  title_name = "CCR2KO"
)

Figure4B <- plot_grid(
  p_WT,
  p_CD169DTR,
  p_CCR2KO,
  nrow = 1,
  align = "hv",
  axis = "tblr"
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure4B_TRM_outgoing_bubble_by_condition.pdf"
  ),
  plot = Figure4B,
  width = 24,
  height = 9,
  units = "in",
  useDingbats = FALSE,
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure4B_TRM_outgoing_bubble_by_condition.png"
  ),
  plot = Figure4B,
  width = 24,
  height = 9,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

get_bubble_table <- function(
  cellchat,
  source_name
) {

  source_index <- group_index(
    cellchat,
    source_name
  )

  target_indices <- present_receiver_indices(
    cellchat
  )

  output_table <- subsetCommunication(
    cellchat,
    sources.use = source_index,
    targets.use = target_indices
  )

  if ("pval" %in% colnames(output_table)) {

    output_table <- output_table[
      !is.na(output_table$pval) &
      output_table$pval <
        bubble_pvalue_threshold,
      ,
      drop = FALSE
    ]
  }

  output_table
}

bubble_table_WT <- get_bubble_table(
  cellchat_WT,
  "TRMs"
)

bubble_table_CD169DTR <- get_bubble_table(
  cellchat_CD169DTR,
  "mo-TRMs"
)

bubble_table_CCR2KO <- get_bubble_table(
  cellchat_CCR2KO,
  "TRMs"
)

save_tsv(
  bubble_table_WT,
  "Figure4B_WT_TRMs_outgoing_pval005.tsv"
)

save_tsv(
  bubble_table_CD169DTR,
  "Figure4B_CD169DTR_moTRMs_outgoing_pval005.tsv"
)

save_tsv(
  bubble_table_CCR2KO,
  "Figure4B_CCR2KO_TRMs_outgoing_pval005.tsv"
)

bubble_summary <- data.frame(
  Condition = c(
    "WT",
    "CD169DTR",
    "CCR2KO"
  ),
  Source = c(
    "TRMs",
    "mo-TRMs",
    "TRMs"
  ),
  LR_pairs_pval_lt_0.05 = c(
    nrow(bubble_table_WT),
    nrow(bubble_table_CD169DTR),
    nrow(bubble_table_CCR2KO)
  ),
  stringsAsFactors = FALSE
)

save_tsv(
  bubble_summary,
  "Figure4B_LR_pair_counts_pval005.tsv"
)

# Diagnóstico específico solicitado:
# CCR2KO TRMs -> Fibroblasts Activated 1.

CCR2KO_TRM_index <- group_index(
  cellchat_CCR2KO,
  "TRMs"
)

CCR2KO_Fib1_index <- group_index(
  cellchat_CCR2KO,
  "Fibroblasts Activated 1"
)

CCR2KO_TRM_to_Fib1_all <- subsetCommunication(
  cellchat_CCR2KO,
  sources.use = CCR2KO_TRM_index,
  targets.use = CCR2KO_Fib1_index
)

save_tsv(
  CCR2KO_TRM_to_Fib1_all,
  "Figure4B_CCR2KO_TRMs_to_Fibroblasts_Activated1_all.tsv"
)

CCR2KO_TRM_to_Fib1_sig <- CCR2KO_TRM_to_Fib1_all

if ("pval" %in% colnames(CCR2KO_TRM_to_Fib1_sig)) {

  CCR2KO_TRM_to_Fib1_sig <- CCR2KO_TRM_to_Fib1_sig[
    !is.na(CCR2KO_TRM_to_Fib1_sig$pval) &
    CCR2KO_TRM_to_Fib1_sig$pval <
      bubble_pvalue_threshold,
    ,
    drop = FALSE
  ]
}

save_tsv(
  CCR2KO_TRM_to_Fib1_sig,
  "Figure4B_CCR2KO_TRMs_to_Fibroblasts_Activated1_pval005.tsv"
)

log_message(
  "CCR2KO TRMs -> Fibroblasts Activated 1: ",
  nrow(CCR2KO_TRM_to_Fib1_all),
  " pares totales; ",
  nrow(CCR2KO_TRM_to_Fib1_sig),
  " con p < ",
  bubble_pvalue_threshold,
  "."
)

# ============================================================
# Figure 4C
# Chord direccional:
#   TRMs de WT + CCR2KO
#   -> Fibroblasts Universal
#
# Ajustes:
#   - flechas desde ligando hacia receptor;
#   - todas las interacciones con probabilidad positiva;
#   - ligandos ordenados por la suma de probabilidad CellChat;
#   - receptores ordenados por probabilidad acumulada;
#   - Pros1 y Axl destacados.
# ============================================================

Figure4C_source_index <- group_index(
  cellchat_WT_CCR2KO,
  "TRMs"
)

Figure4C_target_index <- group_index(
  cellchat_WT_CCR2KO,
  "Fibroblasts Universal"
)

figure4C_table <- subsetCommunication(
  cellchat_WT_CCR2KO,
  sources.use = Figure4C_source_index,
  targets.use = Figure4C_target_index,
  thresh = chord_pvalue_threshold
)

required_chord_columns <- c(
  "ligand",
  "receptor",
  "prob"
)

missing_chord_columns <- setdiff(
  required_chord_columns,
  colnames(figure4C_table)
)

if (length(missing_chord_columns) > 0L) {
  stop(
    "Figure 4C: faltan columnas en subsetCommunication(): ",
    paste(
      missing_chord_columns,
      collapse = ", "
    )
  )
}

figure4C_table <- figure4C_table[
  !is.na(figure4C_table$ligand) &
    !is.na(figure4C_table$receptor) &
    is.finite(figure4C_table$prob) &
    figure4C_table$prob > 0,
  ,
  drop = FALSE
]

# Para el chord se conservan todas las interacciones recuperadas con
# thresh = chord_pvalue_threshold y probabilidad positiva.
# Con chord_pvalue_threshold = 1 entran también las interacciones con
# p >= 0.05, que subsetCommunication() excluiría por defecto.

if (nrow(figure4C_table) == 0L) {
  stop(
    "Figure 4C: no hay interacciones inferidas entre ",
    "TRMs y Fibroblasts Universal."
  )
}

ligand_ranking <- figure4C_table |>
  group_by(ligand) |>
  summarise(
    Total_probability = sum(
      prob,
      na.rm = TRUE
    ),
    Best_pvalue = if (
      "pval" %in% names(cur_data())
    ) {
      min(
        pval,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },
    Number_of_receptors = n_distinct(receptor),
    .groups = "drop"
  ) |>
  arrange(
    desc(Total_probability),
    ligand
  ) |>
  mutate(
    Ligand_rank = row_number()
  ) |>
  select(
    Ligand_rank,
    ligand,
    Total_probability,
    Best_pvalue,
    Number_of_receptors
  )

receptor_ranking <- figure4C_table |>
  group_by(receptor) |>
  summarise(
    Total_probability = sum(
      prob,
      na.rm = TRUE
    ),
    Number_of_ligands = n_distinct(ligand),
    .groups = "drop"
  ) |>
  arrange(
    desc(Total_probability),
    receptor
  ) |>
  mutate(
    Receptor_rank = row_number()
  ) |>
  select(
    Receptor_rank,
    receptor,
    Total_probability,
    Number_of_ligands
  )

save_tsv(
  ligand_ranking,
  "Figure4C_ligand_ranking.tsv"
)

tgfb1_table <- figure4C_table[
  toupper(
    as.character(figure4C_table$ligand)
  ) == "TGFB1",
  ,
  drop = FALSE
]

save_tsv(
  tgfb1_table,
  "Figure4C_TRMs_TGFB1_interactions.tsv"
)

save_tsv(
  receptor_ranking,
  "Figure4C_receptor_ranking.tsv"
)

figure4C_table <- figure4C_table |>
  left_join(
    ligand_ranking |>
      select(
        ligand,
        Ligand_rank,
        Ligand_total_probability = Total_probability
      ),
    by = "ligand"
  ) |>
  left_join(
    receptor_ranking |>
      select(
        receptor,
        Receptor_rank,
        Receptor_total_probability = Total_probability
      ),
    by = "receptor"
  ) |>
  arrange(
    Ligand_rank,
    Receptor_rank,
    desc(prob)
  )

save_tsv(
  figure4C_table,
  "Figure4C_TRMs_to_Universal_communications_ranked.tsv"
)

pros1_axl_table <- figure4C_table[
  toupper(
    as.character(figure4C_table$ligand)
  ) == "PROS1" &
    toupper(
      as.character(figure4C_table$receptor)
    ) == "AXL",
  ,
  drop = FALSE
]

save_tsv(
  pros1_axl_table,
  "Figure4C_TRMs_PROS1_AXL.tsv"
)

ligand_order <- ligand_ranking$ligand
receptor_order <- receptor_ranking$receptor

chord_links <- figure4C_table |>
  transmute(
    from = paste0(
      "TRM_LIGAND::",
      ligand
    ),
    to = paste0(
      "UNIVERSAL_RECEPTOR::",
      receptor
    ),
    value = prob
  ) |>
  group_by(
    from,
    to
  ) |>
  summarise(
    value = sum(
      value,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

ligand_sectors <- paste0(
  "TRM_LIGAND::",
  ligand_order
)

receptor_sectors <- paste0(
  "UNIVERSAL_RECEPTOR::",
  receptor_order
)

sector_order <- c(
  ligand_sectors,
  receptor_sectors
)

sector_colors <- c(
  setNames(
    rep(
      figure4_colors[["TRMs"]],
      length(ligand_sectors)
    ),
    ligand_sectors
  ),
  setNames(
    rep(
      figure4_colors[["Fibroblasts Universal"]],
      length(receptor_sectors)
    ),
    receptor_sectors
  )
)

strip_sector_prefix <- function(
  sector_name
) {
  sub(
    "^[^:]+::",
    "",
    sector_name
  )
}

draw_directional_chord <- function(
  filename,
  device = c(
    "pdf",
    "png"
  )
) {

  device <- match.arg(device)

  if (device == "pdf") {
    pdf(
      filename,
      width = 10,
      height = 10,
      useDingbats = FALSE
    )
  } else {
    png(
      filename,
      width = 3000,
      height = 3000,
      res = 300
    )
  }

  circos.clear()

  gap_vector <- rep(
    1.5,
    length(sector_order)
  )

  gap_vector[
    length(ligand_sectors)
  ] <- 10

  gap_vector[
    length(sector_order)
  ] <- 10

  circos.par(
    start.degree = 90,
    gap.degree = gap_vector,
    track.margin = c(
      0.01,
      0.01
    ),
    cell.padding = c(
      0,
      0,
      0,
      0
    )
  )

  chordDiagram(
    x = chord_links,
    order = sector_order,
    grid.col = sector_colors,
    col = figure4_colors[["TRMs"]],
    transparency = 0.55,
    directional = 1,
    direction.type = c(
      "diffHeight",
      "arrows"
    ),
    diffHeight = mm_h(2),
    link.arr.type = "big.arrow",
    link.sort = TRUE,
    link.decreasing = TRUE,
    annotationTrack = "grid",
    preAllocateTracks = list(
      track.height = 0.14
    )
  )

  circos.trackPlotRegion(
    track.index = 1,
    bg.border = NA,
    panel.fun = function(
      x,
      y
    ) {

      sector_name <- get.cell.meta.data(
        "sector.index"
      )

      displayed_gene <- strip_sector_prefix(
        sector_name
      )

      x_center <- get.cell.meta.data(
        "xcenter"
      )

      y_limits <- get.cell.meta.data(
        "ylim"
      )

      label_color <- "black"
      label_font <- 1

      if (
        toupper(displayed_gene) == "PROS1"
      ) {
        label_color <- figure4_colors[["TRMs"]]
        label_font <- 2
      }

      if (
        toupper(displayed_gene) == "AXL"
      ) {
        label_color <- figure4_colors[[
          "Fibroblasts Activated 2"
        ]]
        label_font <- 2
      }

      circos.text(
        x = x_center,
        y = mean(y_limits),
        labels = displayed_gene,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(
          0,
          0.5
        ),
        cex = 0.72,
        col = label_color,
        font = label_font
      )
    }
  )

  legend(
    "topright",
    legend = c(
      "TRM ligands",
      "Universal fibroblast receptors"
    ),
    fill = c(
      figure4_colors[["TRMs"]],
      figure4_colors[["Fibroblasts Universal"]]
    ),
    border = NA,
    bty = "n",
    cex = 0.9
  )

  title(
    main = paste0(
      "TRMs to Universal fibroblasts\n",
      "Ligands ordered by CellChat probability"
    ),
    cex.main = 1
  )

  circos.clear()
  dev.off()
}

draw_directional_chord(
  filename = file.path(
    figure_dir,
    "Figure4C_TRMs_to_Universal_chord_directional.pdf"
  ),
  device = "pdf"
)

draw_directional_chord(
  filename = file.path(
    figure_dir,
    "Figure4C_TRMs_to_Universal_chord_directional.png"
  ),
  device = "png"
)

# ============================================================
# Información de sesión
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    table_dir,
    "Figure4_CellChat_sessionInfo.txt"
  )
)

log_message(
  "Análisis terminado correctamente."
)

log_message(
  "Figuras: ",
  figure_dir
)

log_message(
  "Tablas: ",
  table_dir
)

log_message(
  "Objetos CellChat: ",
  cellchat_dir
)