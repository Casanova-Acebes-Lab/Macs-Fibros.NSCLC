# ============================================================

# 07_Join.Clustering.R

# Join final de anotaciones por barcode

#

# Cell_Type:

#   Macros$macro_annotation

#   CD45NP$Clusters

#   Tcells$Cluster

#

# Cell_Type2:

#   Macros$macro_annotation

#   CD45NP$fibroblast_annotation

#   Tcells$Cluster

#

# Las B cells y los neutrofilos se recuperan al final desde

# los clusters globales, pero nunca sobrescriben una anotacion

# transferida desde los objetos fuente.

# ============================================================

suppressPackageStartupMessages({

  library(Seurat)

  library(ggplot2)

})

# ============================================================

# Paletas

# ============================================================





Nora_colors2 <- c(

  # ----------------

  # Mieloides

  # ----------------

  "Monocytes Ly6cLo" = "#FF4D4D",

  "Monocytes Ly6cHi" = "#B22222",

  "MoMacs" = "#2A9D8F",

  "em-AMs" = "#1B4332",

  "mo-AMs" = "#52B788",

  "mregDC" = "#F4A261",

  "DC1" = "#6A4C93",

  "DC2" = "#B565A7",

  "pDC" = "#75a2b6",

  "Neutrophils" = "#E9C46A",

  "B cells" = "#4CC9F0",

  # ----------------

  # CD4

  # ----------------

  "CD4 Early Activated" = "#4b7591",

  "CD4 Naive" = "#2ECC71",

  "Cd4 IL7R Memory Like" = "#008F5A",

  "Treg" = "#FFFF00",

  # ----------------

  # CD8 T cells (nueva gama separada)

  # ----------------

  "Immature T Cells"   = "#FFC300",

  "Stem like T Cells"  = "#00CED1",

  "Tpex Transitional"               = "#FF5A36",

  "CD8 Memory" = "#C95F82",

  "Cd8 Exhausted"      = "#6A0DAD",

  # ----------------

  # Innate-like

  # ----------------

  "NK" = "#1F77B4",

  "NKT" = "#003566",

  "ILC" = "#C77D00",

  "Tgd" = "#80ED99",

  # ----------------

  # Epiteliales

  # ----------------

  "Endothelial Vascular" = "#4A90E2",

  "Aerocytes" = "#1B4F72",

  "Epithelial Krt8+" = "#43AA8B",

  "Epithelial AT1" = "#2D6A4F",

  "Epithelial AT2" = "#95D5B2",

  "Ciliated" = "#9B5DE5",

  "Mesothelial" = "#F9C74F",

  "Prolif. Cancer Cells" = "#a8bfce",

  "EMT Cancer Cells" = "#bb6666",

  # ----------------

  # Fibroblastos

  # ----------------

  "Fibroblasts Universal"   = "#FBCDC7",

  "Fibroblasts Activated 1" = "#F94144",

  "Fibroblasts Activated 2" = "#6D1A1A"

)





Nora_colors <- c(

  # ----------------

  # Mieloides

  # ----------------

  "Monocytes Ly6cLo" = "#FF4D4D",

  "Monocytes Ly6cHi" = "#B22222",

  "MoMacs" = "#2A9D8F",

  "em-AMs" = "#1B4332",

  "mo-AMs" = "#52B788",

  "mregDC" = "#F4A261",

  "DC1" = "#6A4C93",

  "DC2" = "#B565A7",

  "pDC" = "#75a2b6",

  "Neutrophils" = "#E9C46A",

  "B cells" = "#4CC9F0",

  # ----------------

  # CD4

  # ----------------

  "CD4 Early Activated" = "#4b7591",

  "CD4 Naive" = "#2ECC71",

  "Cd4 IL7R Memory Like" = "#008F5A",

  "Treg" = "#FFFF00",

  # ----------------

  # CD8 T cells (nueva gama separada)

  # ----------------

  "Immature T Cells"   = "#FFC300",  # amarillo fuerte

  "Stem like T Cells"  = "#00CED1",  # turquesa brillante

  "Tpex Transitional"               = "#FF5A36",  # naranja-rojizo

  "CD8 Memory" = "#C95F82",  # rosa-morado

  "Cd8 Exhausted"      = "#6A0DAD",  # púrpura oscuro

  # ----------------

  # Innate-like

  # ----------------

  "NK" = "#1F77B4",

  "NKT" = "#003566",

  "ILC" = "#C77D00",

  "Tgd" = "#80ED99",

  # ----------------

  # Epiteliales

  # ----------------

  "Endothelial Vascular" = "#4A90E2",

  "Aerocytes" = "#1B4F72",

  "Epithelial Krt8+" = "#43AA8B",

  "Epithelial AT1" = "#2D6A4F",

  "Epithelial AT2" = "#95D5B2",

  "Ciliated" = "#9B5DE5",

  "Prolif. Cancer Cells" = "#a8bfce",

  "EMT Cancer Cells" = "#bb6666",

  "Mesenchymal" = "#FFD6A5",





  # ----------------

  # Fibroblastos

  # ----------------

  "Fibroblasts Universal"   = "#FBCDC7",

  "Fibroblasts Activated 1" = "#F94144",

  "Fibroblasts Activated 2" = "#6D1A1A",

  # ----------------

  # Mesoteliales

  # ----------------

  "Mesothelial" = "#F9C74F"

)













# ============================================================

# Configuracion

# ============================================================

config <- read.delim(

  "config.tsv",

stringsAsFactors = FALSE,

check.names = FALSE

)

if (!all(c("key", "value") %in% colnames(config))) {

  stop(

    "config.tsv must contain the columns 'key' and 'value'."

  )

}

values <- setNames(config$value, config$key)

if (!all(c("rsdir", "outdir") %in% names(values))) {

  stop(

    "config.tsv must contain the keys 'rsdir' and 'outdir'."

  )

}

# ============================================================

# Directorios

# ============================================================

rsdir <- unname(values[["rsdir"]])

outdir <- unname(values[["outdir"]])

objects_dir <- file.path(

  rsdir,

  "objects"

)

join_dir <- file.path(

  outdir,

  "Join"

)

dir.create(

  objects_dir,

recursive = TRUE,

showWarnings = FALSE

)

dir.create(

  join_dir,

recursive = TRUE,

showWarnings = FALSE

)

cat("\nDirectorios del Join:\n")

cat("Objetos: ", objects_dir, "\n")

cat("Plots:   ", join_dir, "\n")

# ============================================================

# Archivos de entrada

# ============================================================

input_files <- c(

General = file.path(

    objects_dir,

    "fibroblasts.QC.rds"

  ),

Macros = file.path(

    objects_dir,

    "fibroblasts.Macrophage_DC.annotated.rds"

  ),

CD45NP = file.path(

    objects_dir,

    "fibroblasts.CD45NP.annotated.rds"

  ),

Tcells = file.path(

    objects_dir,

    "fibroblasts.CD45.Tcells.annotated.rds"

  )

)

missing_files <- input_files[

  !file.exists(input_files)

]

if (length(missing_files) > 0) {

  stop(

    "Missing input files:\n",

    paste(

      names(missing_files),

      missing_files,

sep = ": ",

collapse = "\n"

    )

  )

}

# ============================================================

# Cargar objetos

# ============================================================

fibroblasts.Not.Integrated <- readRDS(

  input_files[["General"]]

)

Macros <- readRDS(

  input_files[["Macros"]]

)

Cd45neg <- readRDS(

  input_files[["CD45NP"]]

)

Tcells <- readRDS(

  input_files[["Tcells"]]

)

cat(

  "\nLoaded objects:\n",

  "General: ", ncol(fibroblasts.Not.Integrated), "\n",

  "Macros/DC: ", ncol(Macros), "\n",

  "CD45 negative: ", ncol(Cd45neg), "\n",

  "T/NK/ILC: ", ncol(Tcells), "\n",

sep = ""

)

# ============================================================

# Preparar columnas de anotación

# ============================================================

if (!"subcluster" %in% colnames(Macros@meta.data)) {

  stop(

    "El objeto Macros no contiene la columna 'subcluster'."

  )

}

if (!"cancer_annotation" %in% colnames(Cd45neg@meta.data)) {

  stop(

    paste0(

      "El objeto Cd45neg no contiene la columna ",

      "'cancer_annotation'."

    )

  )

}

if (!"Cluster" %in% colnames(Tcells@meta.data)) {

  stop(

    "El objeto Tcells no contiene la columna 'Cluster'."

  )

}

Macros$Clustering <- as.character(

  Macros$subcluster

)

Cd45neg$Clustering <- as.character(

  Cd45neg$cancer_annotation

)

Tcells$Clustering <- as.character(

  Tcells$Cluster

)

cat("\nAnotación de Macros/DC utilizada en el Join:\n")

print(

  table(

    Macros$Clustering,

useNA = "ifany"

  )

)

cat("\nAnotación CD45 negativa utilizada en el Join:\n")

print(

  table(

    Cd45neg$Clustering,

useNA = "ifany"

  )

)

cat("\nAnotación de T cells utilizada en el Join:\n")

print(

  table(

    Tcells$Clustering,

useNA = "ifany"

  )

)





# ============================================================

# Añadir anotación de Macros/DC

# ============================================================

fibroblasts.Not.Integrated$Clustering <- NA_character_

macro_cells <- intersect(

  colnames(Macros),

  colnames(fibroblasts.Not.Integrated)

)

fibroblasts.Not.Integrated@meta.data[

  macro_cells,

  "Clustering"

] <- as.character(

  Macros@meta.data[

    macro_cells,

    "Clustering"

  ]

)

cat("\nCélulas en Macros:", ncol(Macros), "\n")

cat("Células transferidas:", length(macro_cells), "\n")

cat(

  "Células de Macros no encontradas:",

  ncol(Macros) - length(macro_cells),

  "\n"

)

cat("\nAnotaciones transferidas desde Macros:\n")

print(

  table(

    fibroblasts.Not.Integrated@meta.data[

      macro_cells,

      "Clustering"

    ],

useNA = "ifany"

  )

)













# ============================================================

# Añadir anotación CD45 negativa

# ============================================================

cd45neg_cells <- intersect(

  colnames(Cd45neg),

  colnames(fibroblasts.Not.Integrated)

)

fibroblasts.Not.Integrated@meta.data[

  cd45neg_cells,

  "Clustering"

] <- as.character(

  Cd45neg@meta.data[

    cd45neg_cells,

    "Clustering"

  ]

)

cat("\nCélulas en Cd45neg:", ncol(Cd45neg), "\n")

cat("Células transferidas:", length(cd45neg_cells), "\n")

cat(

  "Células de Cd45neg no encontradas:",

  ncol(Cd45neg) - length(cd45neg_cells),

  "\n"

)

cat("\nAnotaciones transferidas desde Cd45neg:\n")

print(

  table(

    fibroblasts.Not.Integrated@meta.data[

      cd45neg_cells,

      "Clustering"

    ],

useNA = "ifany"

  )

)

# Adding Tcells 





# Células presentes en ambos objetos

tcell_cells <- intersect(

  colnames(Tcells),

  colnames(fibroblasts.Not.Integrated)

)

# Transferir la anotación a la misma columna

fibroblasts.Not.Integrated@meta.data[

  tcell_cells,

  "Clustering"

] <- as.character(

  Tcells@meta.data[

    tcell_cells,

    "Cluster"

  ]

)

# Comprobaciones

cat("Células en Tcells:", ncol(Tcells), "\n")

cat("Células transferidas:", length(tcell_cells), "\n")

cat(

  "Células de Tcells no encontradas:",

  ncol(Tcells) - length(tcell_cells),

  "\n"

)

table(

  fibroblasts.Not.Integrated$Clustering,

useNA = "ifany"

)













# ============================================================

# AÑADIR B CELLS Y NEUTRÓFILOS A CLUSTERING

#

# Se recuperan desde los clusters globales únicamente para

# células que siguen sin anotación. Nunca se sobrescriben las

# anotaciones transferidas desde Macros/DC, CD45NP o T cells.

# ============================================================

unannotated_cells <- is.na(

  fibroblasts.Not.Integrated$Clustering

)

b_cell_mask <- (

  unannotated_cells &

    fibroblasts.Not.Integrated$seurat_clusters %in% c(5, 10, 22)

)

fibroblasts.Not.Integrated$Clustering[

  b_cell_mask

] <- "B cells"

# Recalcular NA después de recuperar B cells.

unannotated_cells <- is.na(

  fibroblasts.Not.Integrated$Clustering

)

neutrophil_mask <- (

  unannotated_cells &

    fibroblasts.Not.Integrated$seurat_clusters == 4

)

fibroblasts.Not.Integrated$Clustering[

  neutrophil_mask

] <- "Neutrophils"

cat(

  "\nB cells recuperadas desde clusters globales: ",

  sum(b_cell_mask),

  "\n",

sep = ""

)

cat(

  "Neutrófilos recuperados desde cluster global 4: ",

  sum(neutrophil_mask),

  "\n",

sep = ""

)

table(

  fibroblasts.Not.Integrated$Clustering,

useNA = "ifany"

)

# ============================================================

# Crear anotación con compartimento mesenquimal agrupado

# ============================================================

fibroblasts.Not.Integrated$Clustering.Mesenchymal <-

  as.character(

    fibroblasts.Not.Integrated$Clustering

  )

mesenchymal_labels <- c(

  "Fibroblasts Universal",

  "Fibroblasts Activated 1",

  "Fibroblasts Activated 2",

  "Mesothelial"

)

fibroblasts.Not.Integrated$Clustering.Mesenchymal[

  fibroblasts.Not.Integrated$Clustering %in%

    mesenchymal_labels

] <- "Mesenchymal"

cat("\nClustering original:\n")

print(

  table(

    fibroblasts.Not.Integrated$Clustering,

useNA = "ifany"

  )

)

cat("\nClustering con compartimento mesenquimal agrupado:\n")

print(

  table(

    fibroblasts.Not.Integrated$Clustering.Mesenchymal,

useNA = "ifany"

  )

)













# ============================================================

# ELIMINAR CÉLULAS SIN ANOTACIÓN FINAL

# ============================================================

cells_with_na <- colnames(

  fibroblasts.Not.Integrated

)[

  is.na(fibroblasts.Not.Integrated$Clustering) |

    is.na(

      fibroblasts.Not.Integrated$Clustering.Mesenchymal

    )

]

cat(

  "\nCélulas eliminadas por contener NA en alguna anotación:",

  length(cells_with_na),

  "\n"

)

# Guardar sus metadatos antes de eliminarlas

if (length(cells_with_na) > 0) {

  excluded_na_metadata <-

    fibroblasts.Not.Integrated@meta.data[

      cells_with_na,

      ,

      drop = FALSE

    ]

  excluded_na_metadata$barcode <-

    rownames(excluded_na_metadata)

  excluded_na_metadata <- excluded_na_metadata[

    ,

    c(

      "barcode",

      setdiff(

        colnames(excluded_na_metadata),

        "barcode"

      )

    ),

    drop = FALSE

  ]

  write.table(

    excluded_na_metadata,

file = file.path(

      join_dir,

      "excluded_cells_with_NA_annotation.tsv"

    ),

sep = "\t",

quote = FALSE,

row.names = FALSE

  )

  # Eliminar del objeto principal

  fibroblasts.Not.Integrated <- subset(

    fibroblasts.Not.Integrated,

cells = setdiff(

      colnames(fibroblasts.Not.Integrated),

      cells_with_na

    )

  )

}

# ============================================================

# COMPROBACIONES

# ============================================================

if (

  anyNA(

    fibroblasts.Not.Integrated$Clustering

  )

) {

  stop(

    "Todavía quedan valores NA en Clustering."

  )

}

if (

  anyNA(

    fibroblasts.Not.Integrated$

      Clustering.Mesenchymal

  )

) {

  stop(

    paste0(

      "Todavía quedan valores NA en ",

      "Clustering.Mesenchymal."

    )

  )

}

cat(

  "Células finales retenidas:",

  ncol(fibroblasts.Not.Integrated),

  "\n"

)

cat(

  "\nClustering final:\n"

)

print(

  table(

    fibroblasts.Not.Integrated$Clustering,

useNA = "ifany"

  )

)

cat(

  "\nClustering.Mesenchymal final:\n"

)

print(

  table(

    fibroblasts.Not.Integrated$

      Clustering.Mesenchymal,

useNA = "ifany"

  )

)



# ============================================================

# TARGETED CLEANING USING GLOBAL SEURAT CLUSTERS

#

# Only explicitly identified contamination is removed.

#

# For each target population:

#   1. identify the global Seurat cluster(s) dominated by it

#   2. remove only the specified contaminating annotations

#      inside those clusters

#

# No nearest-neighbour analysis is performed.

# ============================================================

cat(

  "\n========================================\n",

  "Targeted cleaning by Seurat cluster\n",

  "========================================\n",

sep = ""

)

# ============================================================

# Explicit contamination rules

# ============================================================

cleaning_rules <- list(

  # ----------------------------------------------------------

  # Immature T cells

  # ----------------------------------------------------------

  "Immature T Cells" = c(

    "Stem like T Cells"

  ),

  # ----------------------------------------------------------

  # Stem-like T cells

  # ----------------------------------------------------------

  "Stem like T Cells" = c(

    "CD4 Naive",

    "CD4 Early Activated",

    "Immature T Cells",

    "CD8 Memory",

    "Tpex Transitional"

  ),

  # ----------------------------------------------------------

  # Tpex

  # ----------------------------------------------------------

  "Tpex Transitional" = c(

    "Stem like T Cells"

  ),

  # ----------------------------------------------------------

  # B cells

  # ----------------------------------------------------------

  "B cells" = c(

    "CD8 Memory"

  ),

  # ----------------------------------------------------------

  # Neutrophils

  # ----------------------------------------------------------

  "Neutrophils" = c(

    "CD4 Naive",

    "CD8 Memory",

    "Tgd",

    "Endothelial Vascular"

  ),

  # ----------------------------------------------------------

  # Tgd

  # ----------------------------------------------------------

  "Tgd" = c(

    "Stem like T Cells",

    "CD8 Memory",

    "Mesothelial"

  ),

  # ----------------------------------------------------------

  # Aerocytes

  # ----------------------------------------------------------

  "Aerocytes" = c(

    "Endothelial Vascular"

  ),

  # ----------------------------------------------------------

  # Endothelial vascular

  # ----------------------------------------------------------

  "Endothelial Vascular" = c(

    "EMT Cancer Cells"

  )

)

# ============================================================

# Cluster x annotation table

# ============================================================

cluster_annotation_table <- table(

seurat_cluster = as.character(

    fibroblasts.Not.Integrated$seurat_clusters

  ),

annotation = as.character(

    fibroblasts.Not.Integrated$Clustering

  )

)

cat(

  "\nSeurat cluster x final annotation:\n"

)

print(

  cluster_annotation_table

)

write.table(

  as.data.frame(

    cluster_annotation_table

  ),

file = file.path(

    join_dir,

    "seurat_cluster_vs_annotation_before_cleaning.tsv"

  ),

sep = "\t",

quote = FALSE,

row.names = FALSE

)

# ============================================================

# Determine dominant annotation of each Seurat cluster

# ============================================================

cluster_size <- rowSums(

  cluster_annotation_table

)

dominant_index <- max.col(

  cluster_annotation_table,

ties.method = "first"

)

dominant_annotation <- colnames(

  cluster_annotation_table

)[

  dominant_index

]

dominant_n <- cluster_annotation_table[

  cbind(

    seq_len(

      nrow(cluster_annotation_table)

    ),

    dominant_index

  )

]

dominant_fraction <- dominant_n /

  cluster_size

cluster_info <- data.frame(

seurat_cluster = rownames(

    cluster_annotation_table

  ),

dominant_annotation =

    dominant_annotation,

dominant_n =

    dominant_n,

cluster_size =

    cluster_size,

dominant_fraction =

    dominant_fraction,

stringsAsFactors = FALSE

)

cat(

  "\nDominant annotation of global clusters:\n"

)

print(

  cluster_info,

row.names = FALSE

)

write.table(

  cluster_info,

file = file.path(

    join_dir,

    "seurat_cluster_dominant_annotation.tsv"

  ),

sep = "\t",

quote = FALSE,

row.names = FALSE

)

# ============================================================

# Select contaminating cells

#

# We only clean clusters where the target annotation is:

#   - the dominant annotation

#   - at least 60% of the cluster

#

# This prevents cleaning genuinely mixed clusters.

# ============================================================

min_target_cluster_fraction <- 0.60

cells_to_remove <- character(0)

cleaning_summary_list <- list()

for (

  target_annotation

  in names(cleaning_rules)

) {

  contaminating_annotations <-

    cleaning_rules[[

      target_annotation

    ]]

  target_clusters <- cluster_info$seurat_cluster[

    cluster_info$dominant_annotation ==

      target_annotation &

      cluster_info$dominant_fraction >=

        min_target_cluster_fraction

  ]

  cat(

    "\n",

    target_annotation,

    " dominant Seurat cluster(s): ",

    paste(

      target_clusters,

collapse = ", "

    ),

    "\n",

sep = ""

  )

  if (length(target_clusters) == 0) {

    cat(

      "No suitable dominant cluster found; skipping.\n"

    )

    next

  }

  current_cells <- colnames(

    fibroblasts.Not.Integrated

  )[

    as.character(

      fibroblasts.Not.Integrated$seurat_clusters

    ) %in% target_clusters &

      as.character(

        fibroblasts.Not.Integrated$Clustering

      ) %in% contaminating_annotations

  ]

  if (length(current_cells) == 0) {

    next

  }

  current_metadata <- data.frame(

barcode = current_cells,

seurat_cluster = as.character(

      fibroblasts.Not.Integrated$seurat_clusters[

        match(

          current_cells,

          colnames(

            fibroblasts.Not.Integrated

          )

        )

      ]

    ),

annotation = as.character(

      fibroblasts.Not.Integrated$Clustering[

        match(

          current_cells,

          colnames(

            fibroblasts.Not.Integrated

          )

        )

      ]

    ),

target_cluster =

      target_annotation,

stringsAsFactors = FALSE

  )

  cleaning_summary_list[[

    length(cleaning_summary_list) + 1L

  ]] <- current_metadata

  cells_to_remove <- c(

    cells_to_remove,

    current_cells

  )

}

cells_to_remove <- unique(

  cells_to_remove

)

# ============================================================

# Cleaning summary

# ============================================================

cat(

  "\nCells selected for targeted cleaning: ",

  length(cells_to_remove),

  "\n",

sep = ""

)

if (length(cleaning_summary_list) > 0) {

  cleaning_cells <- do.call(

    rbind,

    cleaning_summary_list

  )

  cleaning_summary <- aggregate(

    barcode ~

      target_cluster +

      seurat_cluster +

      annotation,

data = cleaning_cells,

FUN = length

  )

  colnames(

    cleaning_summary

  )[

    colnames(cleaning_summary) ==

      "barcode"

  ] <- "n_cells"

  cat(

    "\nCells removed by rule:\n"

  )

  print(

    cleaning_summary,

row.names = FALSE

  )

  write.table(

    cleaning_summary,

file = file.path(

      join_dir,

      "targeted_cluster_cleaning_summary.tsv"

    ),

sep = "\t",

quote = FALSE,

row.names = FALSE

  )

  write.table(

    cleaning_cells,

file = file.path(

      join_dir,

      "targeted_cluster_cleaning_cells.tsv"

    ),

sep = "\t",

quote = FALSE,

row.names = FALSE

  )

}

# ============================================================

# Palette for diagnostic plots

# ============================================================

cleaning_colors <- Nora_colors2

present_cleaning_labels <- sort(
  unique(
    as.character(
      fibroblasts.Not.Integrated$Clustering
    )
  )
)

missing_cleaning_colors <- setdiff(
  present_cleaning_labels,
  names(cleaning_colors)
)

if (length(missing_cleaning_colors) > 0) {
  stop(
    "Faltan colores en Nora_colors2 para estas anotaciones:\n",
    paste(
      missing_cleaning_colors,
      collapse = "\n"
    )
  )
}



# ============================================================

# UMAP BEFORE CLEANING

# Cells to remove highlighted

# ============================================================

umap_coordinates <- Embeddings(

  fibroblasts.Not.Integrated,

reduction = "umap"

)

diagnostic_umap <- data.frame(

Cell = rownames(

    umap_coordinates

  ),

UMAP_1 = umap_coordinates[, 1],

UMAP_2 = umap_coordinates[, 2],

Clustering = as.character(

    fibroblasts.Not.Integrated$Clustering[

      rownames(umap_coordinates)

    ]

  ),

orig.ident = as.character(

    fibroblasts.Not.Integrated$orig.ident[

      rownames(umap_coordinates)

    ]

  ),

stringsAsFactors = FALSE

)

diagnostic_umap$remove <-

  diagnostic_umap$Cell %in%

    cells_to_remove

p_cleaning_candidates <- ggplot(

  diagnostic_umap,

  aes(

    UMAP_1,

    UMAP_2

  )

) +

  geom_point(

data = diagnostic_umap[

      !diagnostic_umap$remove,

      ,

      drop = FALSE

    ],

color = "grey85",

size = 0.2,

alpha = 0.4

  ) +

  geom_point(

data = diagnostic_umap[

      diagnostic_umap$remove,

      ,

      drop = FALSE

    ],

    aes(

color = Clustering

    ),

size = 1.1

  ) +

  scale_color_manual(

values = cleaning_colors

  ) +

  coord_equal() +

  theme_classic() +

  ggtitle(

    "Cells removed by targeted cleaning"

  )

ggsave(

  file.path(

    join_dir,

    "UMAP_targeted_cleaning_candidates.png"

  ),

  p_cleaning_candidates,

width = 11,

height = 8,

dpi = 300

)

ggsave(

  file.path(

    join_dir,

    "UMAP_targeted_cleaning_candidates.pdf"

  ),

  p_cleaning_candidates,

width = 11,

height = 8,

useDingbats = FALSE

)

# ============================================================

# REMOVE CELLS

# ============================================================

if (length(cells_to_remove) > 0) {

  fibroblasts.Not.Integrated <- subset(

    fibroblasts.Not.Integrated,

cells = setdiff(

      colnames(

        fibroblasts.Not.Integrated

      ),

      cells_to_remove

    )

  )

}

cat(

  "\nCells retained after cleaning: ",

  ncol(fibroblasts.Not.Integrated),

  "\n",

sep = ""

)



# ============================================================

# RESIDUAL CLEANING OF SPECIFIC CONTAMINANTS

#

# Segunda pasada únicamente para las contaminaciones que pueden

# quedar después de la limpieza principal por seurat_clusters.

#

# Se mantienen todas las reglas anteriores.

#

# Residuales:

#   Tgd                  -> Neutrophils

#   CD8 Memory           -> Neutrophils

#   Endothelial Vascular -> Neutrophils

#   Endothelial Vascular -> Aerocytes

# ============================================================

cat(

  "\n========================================\n",

  "Residual targeted cleaning\n",

  "========================================\n",

sep = ""

)

residual_rules <- data.frame(

annotation = c(

    "Tgd",

    "CD8 Memory",

    "Endothelial Vascular",

    "CD4 Naive",

    "Aerocytes",

    "Endothelial Vascular"

  ),

target_annotation = c(

    "Neutrophils",

    "Neutrophils",

    "Neutrophils",

    "Neutrophils",

    "Neutrophils",

    "Aerocytes"

  ),

stringsAsFactors = FALSE

)

residual_k <- 20L

residual_min_target_fraction <- 0.65

if (!requireNamespace("FNN", quietly = TRUE)) {

  stop(

    "El paquete R 'FNN' es necesario para la limpieza residual."

  )

}

annotation_residual <- as.character(

  fibroblasts.Not.Integrated$Clustering

)

names(annotation_residual) <- colnames(

  fibroblasts.Not.Integrated

)

umap_residual <- Embeddings(

  fibroblasts.Not.Integrated,

reduction = "umap"

)

cells_residual_remove <- character(0)

residual_results <- list()

for (rule_i in seq_len(nrow(residual_rules))) {

  source_annotation <- residual_rules$annotation[

    rule_i

  ]

  target_annotation <- residual_rules$target_annotation[

    rule_i

  ]

  # ----------------------------------------------------------

  # Compare only source and target populations

  # ----------------------------------------------------------

  pair_cells <- names(annotation_residual)[

    annotation_residual %in%

      c(

        source_annotation,

        target_annotation

      )

  ]

  source_cells <- names(annotation_residual)[

    annotation_residual ==

      source_annotation

  ]

  if (

    length(source_cells) == 0 ||

      length(pair_cells) <= residual_k

  ) {

    next

  }

  pair_umap <- umap_residual[

    pair_cells,

    1:2,

    drop = FALSE

  ]

  pair_annotation <- annotation_residual[

    pair_cells

  ]

  knn_pair <- FNN::get.knn(

    pair_umap,

k = residual_k

  )

  # ----------------------------------------------------------

  # Evaluate source cells

  # ----------------------------------------------------------

  for (cell_name in source_cells) {

    cell_index <- match(

      cell_name,

      pair_cells

    )

    if (is.na(cell_index)) {

      next

    }

    neighbour_cells <- pair_cells[

      knn_pair$nn.index[

        cell_index,

      ]

    ]

    neighbour_annotations <- pair_annotation[

      neighbour_cells

    ]

    target_fraction <- mean(

      neighbour_annotations ==

        target_annotation

    )

    remove_cell <- (

      target_fraction >=

        residual_min_target_fraction

    )

    residual_results[[

      length(residual_results) + 1L

    ]] <- data.frame(

barcode = cell_name,

annotation = source_annotation,

target_annotation = target_annotation,

target_fraction = target_fraction,

remove = remove_cell,

stringsAsFactors = FALSE

    )

    if (remove_cell) {

      cells_residual_remove <- c(

        cells_residual_remove,

        cell_name

      )

    }

  }

}

cells_residual_remove <- unique(

  cells_residual_remove

)

cat(

  "\nResidual contaminating cells detected: ",

  length(cells_residual_remove),

  "\n",

sep = ""

)

# ============================================================

# Residual cleaning summary

# ============================================================

if (length(residual_results) > 0) {

  residual_results <- do.call(

    rbind,

    residual_results

  )

  residual_selected <- residual_results[

    residual_results$remove,

    ,

    drop = FALSE

  ]

  if (nrow(residual_selected) > 0) {

    residual_summary <- aggregate(

      barcode ~

        annotation +

        target_annotation,

data = residual_selected,

FUN = length

    )

    colnames(residual_summary)[

      colnames(residual_summary) ==

        "barcode"

    ] <- "n_cells"

    cat(

      "\nResidual cells removed:\n"

    )

    print(

      residual_summary,

row.names = FALSE

    )

    write.table(

      residual_summary,

file = file.path(

        join_dir,

        "targeted_residual_cleaning_summary.tsv"

      ),

sep = "\t",

quote = FALSE,

row.names = FALSE

    )

    write.table(

      residual_selected,

file = file.path(

        join_dir,

        "targeted_residual_cleaning_cells.tsv"

      ),

sep = "\t",

quote = FALSE,

row.names = FALSE

    )

  }

}

# ============================================================

# Remove residual contaminants

# ============================================================

if (length(cells_residual_remove) > 0) {

  fibroblasts.Not.Integrated <- subset(

    fibroblasts.Not.Integrated,

cells = setdiff(

      colnames(

        fibroblasts.Not.Integrated

      ),

      cells_residual_remove

    )

  )

}

cat(

  "\nCells retained after residual cleaning: ",

  ncol(fibroblasts.Not.Integrated),

  "\n",

sep = ""

)

cat(

  "\nCells removed by main targeted cleaning: ",

  length(cells_to_remove),

  "\n",

sep = ""

)

cat(

  "Cells removed by residual cleaning: ",

  length(cells_residual_remove),

  "\n",

sep = ""

)



# ============================================================

# CLEAN UMAP

# ============================================================

p_umap_cleaned <- DimPlot(

  fibroblasts.Not.Integrated,

reduction = "umap",

group.by = "Clustering",

cols = cleaning_colors,

label = TRUE,

repel = TRUE,

label.size = 3,

pt.size = 0.3,

shuffle = TRUE,

seed = 42

) +

  ggtitle(

    "Final cell annotation after targeted cleaning"

  ) +

  theme(

plot.title = element_text(

hjust = 0.5,

face = "bold"

    ),

legend.title =

      element_blank()

  )

ggsave(

  file.path(

    join_dir,

    "UMAP_after_targeted_cleaning_annotation.png"

  ),

  p_umap_cleaned,

width = 15,

height = 11,

dpi = 300,

limitsize = FALSE

)

ggsave(

  file.path(

    join_dir,

    "UMAP_after_targeted_cleaning_annotation.pdf"

  ),

  p_umap_cleaned,

width = 15,

height = 11,

useDingbats = FALSE,

limitsize = FALSE

)

# ============================================================

# CLEAN UMAP SPLIT BY CONDITION

# ============================================================

p_umap_cleaned_split <- DimPlot(

  fibroblasts.Not.Integrated,

reduction = "umap",

group.by = "Clustering",

split.by = "orig.ident",

cols = cleaning_colors,

label = FALSE,

pt.size = 0.25,

shuffle = TRUE,

seed = 42,

ncol = 3

) +

  ggtitle(

    "Final cell annotation after targeted cleaning by condition"

  )

ggsave(

  file.path(

    join_dir,

    "UMAP_after_targeted_cleaning_annotation_split_orig.ident.png"

  ),

  p_umap_cleaned_split,

width = 22,

height = 8,

dpi = 300,

limitsize = FALSE

)

ggsave(

  file.path(

    join_dir,

    "UMAP_after_targeted_cleaning_annotation_split_orig.ident.pdf"

  ),

  p_umap_cleaned_split,

width = 22,

height = 8,

useDingbats = FALSE,

limitsize = FALSE

)

# ============================================================

# Guardar objeto final

# ============================================================

Idents(

  fibroblasts.Not.Integrated

) <- "Clustering"

output_file <- file.path(

  objects_dir,

  "fibroblasts.Not.Integrated.annotated.rds"

)

saveRDS(

  fibroblasts.Not.Integrated,

  output_file

)

cat(

  "\nObjeto final guardado en:\n",

  output_file,

  "\n"

)

cat(

  "\nPlots guardados en:\n",

  join_dir,

  "\n"

)