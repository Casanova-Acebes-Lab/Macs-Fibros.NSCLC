# ============================================================
# 06_JOIN.R
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
  "Cd4 Early Activated" = "#4b7591",
  "Cd4 Naive" = "#2ECC71",
  "Cd4 IL7R Memory Like" = "#008F5A",
  "Treg" = "#FFFF00",

  # ----------------
  # CD8 T cells (nueva gama separada)
  # ----------------
  "Inmature T cells"   = "#FFC300",
  "Stem like T cells"  = "#00CED1",
  "Tpex Transitional"               = "#FF5A36",
  "Cd8 Memory" = "#C95F82",
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
  "Cd4 Early Activated" = "#4b7591",
  "Cd4 Naive" = "#2ECC71",
  "Cd4 IL7R Memory Like" = "#008F5A",
  "Treg" = "#FFFF00",

  # ----------------
  # CD8 T cells (nueva gama separada)
  # ----------------
  "Inmature T cells"   = "#FFC300",  # amarillo fuerte
  "Stem like T cells"  = "#00CED1",  # turquesa brillante
  "Tpex Transitional"               = "#FF5A36",  # naranja-rojizo
  "Cd8 Memory" = "#C95F82",  # rosa-morado
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

objects_dir <- file.path(values[["rsdir"]], "objects")
join_dir <- file.path(values[["outdir"]], "Join")

dir.create(join_dir, recursive = TRUE, showWarnings = FALSE)



# ============================================================
# Cargar objetos
# ============================================================

input_files <- c(
  General = file.path(
    objects_dir,
    "fibroblasts.Not.Integrated.rds"
  ),
  Macros = file.path(
    objects_dir,
    "fibroblasts.Macrophage_DC.annotated.rds"
  ),
  CD45NP = file.path(
    objects_dir,
    "fibroblasts.CD45NP.final.rds"
  ),
  Tcells = file.path(
    objects_dir,
    "fibroblasts.CD45.Tcells.annotated.rds"
  )
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing input files:\n",
    paste(missing_files, collapse = "\n")
  )
}

fibroblasts.Not.Integrated <- readRDS(
  input_files[["General"]]
)

Macros <- readRDS(
  input_files[["Macros"]]
)

CD45NP <- readRDS(
  input_files[["CD45NP"]]
)

Tcells <- readRDS(
  input_files[["Tcells"]]
)

cat(
  "\nLoaded objects:\n",
  "General: ", ncol(fibroblasts.Not.Integrated), "\n",
  "Macros/DC: ", ncol(Macros), "\n",
  "CD45 negative: ", ncol(CD45NP), "\n",
  "T/NK/ILC: ", ncol(Tcells), "\n",
  sep = ""
)




 ### Cell Type1
# Paso 1: Leer cada archivo, con barcodes como rownames
Macros <- readRDS(paste0(objects_dir, "/fibroblasts.Macrophage_DC.annotated.rds"))
Cd45neg <- readRDS(paste0(objects_dir, "/fibroblasts.CD45NP.final.rds"))
Tcells <- readRDS(paste0(objects_dir, "/fibroblasts.CD45.Tcells.annotated.rds"))


Macros$Clustering <- Macros$macro_annotation
Cd45neg$Clustering <-  Cd45neg$fibroblast_annotation
Cd45neg$Clustering.broad <-  Cd45neg$Clusters
Tcells$Clustering <-  Tcells$Cluster




# Adding Macros

# Crear la nueva columna una sola vez
fibroblasts.Not.Integrated$Clustering <- NA_character_

# Células presentes en ambos objetos
macro_cells <- intersect(
  colnames(Macros),
  colnames(fibroblasts.Not.Integrated)
)

# Transferir la anotación
fibroblasts.Not.Integrated$Clustering[macro_cells] <-
  as.character(Macros$macro_annotation[macro_cells])

# Comprobaciones
cat("Células en Macros:", ncol(Macros), "\n")
cat("Células transferidas:", length(macro_cells), "\n")
cat(
  "Células de Macros no encontradas:",
  ncol(Macros) - length(macro_cells),
  "\n"
)

table(
  fibroblasts.Not.Integrated$Clustering,
  useNA = "ifany"
)



# Adding Cd45neg

# Células presentes en ambos objetos
cd45neg_cells <- intersect(
  colnames(Cd45neg),
  colnames(fibroblasts.Not.Integrated)
)

# Transferir la anotación a la misma columna
fibroblasts.Not.Integrated$Clustering[cd45neg_cells] <-
  as.character(Cd45neg$fibroblast_annotation[cd45neg_cells])

# Comprobaciones
cat("Células en Cd45neg:", ncol(Cd45neg), "\n")
cat("Células transferidas:", length(cd45neg_cells), "\n")
cat(
  "Células de Cd45neg no encontradas:",
  ncol(Cd45neg) - length(cd45neg_cells),
  "\n"
)

table(
  fibroblasts.Not.Integrated$Clustering,
  useNA = "ifany"
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
# ============================================================

fibroblasts.Not.Integrated@meta.data$Clustering[
  fibroblasts.Not.Integrated@meta.data$seurat_clusters %in% c(5, 10, 22)
] <- "B cells"

fibroblasts.Not.Integrated@meta.data$Clustering[
  fibroblasts.Not.Integrated@meta.data$seurat_clusters == 4
] <- "Neutrophils"

# Comprobar
table(
  fibroblasts.Not.Integrated$Clustering,
  useNA = "ifany"
)


# ============================================================
# CREAR ANOTACIÓN AGRUPADA MESENCHYMAL
# ============================================================

# Copiar la anotación completa
fibroblasts.Not.Integrated$Clustering.Mesenchymal <-
  as.character(fibroblasts.Not.Integrated$Clustering)

# Agrupar fibroblastos y células mesoteliales
mesenchymal_labels <- c(
  "Fibroblasts Universal",
  "Fibroblasts Activated 1",
  "Fibroblasts Activated 2",
  "Mesothelial"
)

fibroblasts.Not.Integrated$Clustering.Mesenchymal[
  fibroblasts.Not.Integrated$Clustering %in% mesenchymal_labels
] <- "Mesenchymal"

# Comprobar
table(
  fibroblasts.Not.Integrated$Clustering.Mesenchymal,
  useNA = "ifany"
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
# GUARDAR OBJETO FINAL
# ============================================================

Idents(fibroblasts.Not.Integrated) <- "Clustering"

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


# ============================================================
# CORREGIR NOMBRES DE LAS PALETAS
# Deben coincidir exactamente con las anotaciones
# ============================================================

names(Nora_colors2)[
  names(Nora_colors2) == "Cd4 Early Activated"
] <- "CD4 Early Activated"

names(Nora_colors2)[
  names(Nora_colors2) == "Cd4 Naive"
] <- "CD4 Naive"

names(Nora_colors2)[
  names(Nora_colors2) == "Cd8 Memory"
] <- "CD8 Memory"

names(Nora_colors2)[
  names(Nora_colors2) == "Inmature T cells"
] <- "Immature T Cells"

names(Nora_colors2)[
  names(Nora_colors2) == "Stem like T cells"
] <- "Stem like T Cells"


names(Nora_colors)[
  names(Nora_colors) == "Cd4 Early Activated"
] <- "CD4 Early Activated"

names(Nora_colors)[
  names(Nora_colors) == "Cd4 Naive"
] <- "CD4 Naive"

names(Nora_colors)[
  names(Nora_colors) == "Cd8 Memory"
] <- "CD8 Memory"

names(Nora_colors)[
  names(Nora_colors) == "Inmature T cells"
] <- "Immature T Cells"

names(Nora_colors)[
  names(Nora_colors) == "Stem like T cells"
] <- "Stem like T Cells"

# ============================================================
# PLOTS FINALES
# ============================================================

if (!"umap" %in% Reductions(fibroblasts.Not.Integrated)) {
  stop(
    "El objeto principal no contiene una reducción llamada 'umap'."
  )
}

# Comprobar si quedan células sin anotación
cat(
  "\nCélulas sin anotación en Clustering:",
  sum(is.na(fibroblasts.Not.Integrated$Clustering)),
  "\n"
)

cat(
  "Células sin anotación en Clustering.Mesenchymal:",
  sum(is.na(fibroblasts.Not.Integrated$Clustering.Mesenchymal)),
  "\n"
)

# ============================================================
# Función para preparar la paleta
#
# Detiene el script indicando exactamente qué anotaciones
# no tienen un color definido.
# ============================================================

prepare_palette <- function(annotation, palette, palette_name) {

  present_labels <- sort(
    unique(
      na.omit(
        as.character(annotation)
      )
    )
  )

  missing_colors <- setdiff(
    present_labels,
    names(palette)
  )

  if (length(missing_colors) > 0) {
    stop(
      "Faltan colores en ",
      palette_name,
      " para las siguientes anotaciones:\n",
      paste(
        missing_colors,
        collapse = "\n"
      )
    )
  }

  palette[
    names(palette) %in% present_labels
  ]
}

# Paleta para la anotación fina
clustering_colors <- prepare_palette(
  annotation = fibroblasts.Not.Integrated$Clustering,
  palette = Nora_colors2,
  palette_name = "Nora_colors2"
)

# Paleta para la anotación mesenquimal agrupada
mesenchymal_colors <- prepare_palette(
  annotation =
    fibroblasts.Not.Integrated$Clustering.Mesenchymal,
  palette = Nora_colors,
  palette_name = "Nora_colors"
)

# ============================================================
# Convertir anotaciones en factores con el orden de las paletas
# ============================================================

fibroblasts.Not.Integrated$Clustering <- factor(
  as.character(
    fibroblasts.Not.Integrated$Clustering
  ),
  levels = names(clustering_colors)
)

fibroblasts.Not.Integrated$Clustering.Mesenchymal <- factor(
  as.character(
    fibroblasts.Not.Integrated$Clustering.Mesenchymal
  ),
  levels = names(mesenchymal_colors)
)

# ============================================================
# UMAP con anotación fina
# ============================================================

p_clustering <- DimPlot(
  object = fibroblasts.Not.Integrated,
  reduction = "umap",
  group.by = "Clustering",
  cols = clustering_colors,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  pt.size = 0.3,
  shuffle = TRUE,
  seed = 42
) +
  ggtitle(
    "Final cell annotation"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  )

# Necesario para que aparezca cuando se ejecuta con source()
print(p_clustering)

ggsave(
  filename = file.path(
    join_dir,
    "UMAP_Clustering.pdf"
  ),
  plot = p_clustering,
  width = 15,
  height = 11,
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    join_dir,
    "UMAP_Clustering.png"
  ),
  plot = p_clustering,
  width = 15,
  height = 11,
  dpi = 300,
  limitsize = FALSE
)

# ============================================================
# UMAP con compartimento mesenquimal agrupado
# ============================================================

p_clustering_mesenchymal <- DimPlot(
  object = fibroblasts.Not.Integrated,
  reduction = "umap",
  group.by = "Clustering.Mesenchymal",
  cols = mesenchymal_colors,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  pt.size = 0.3,
  shuffle = TRUE,
  seed = 42
) +
  ggtitle(
    "Final cell annotation — Mesenchymal grouped"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  )

print(p_clustering_mesenchymal)

ggsave(
  filename = file.path(
    join_dir,
    "UMAP_Clustering_Mesenchymal.pdf"
  ),
  plot = p_clustering_mesenchymal,
  width = 15,
  height = 11,
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    join_dir,
    "UMAP_Clustering_Mesenchymal.png"
  ),
  plot = p_clustering_mesenchymal,
  width = 15,
  height = 11,
  dpi = 300,
  limitsize = FALSE
)

# ============================================================
# UMAP divididos por condición
# ============================================================

if (
  "orig.ident" %in%
    colnames(fibroblasts.Not.Integrated@meta.data)
) {

  p_clustering_condition <- DimPlot(
    object = fibroblasts.Not.Integrated,
    reduction = "umap",
    group.by = "Clustering",
    split.by = "orig.ident",
    cols = clustering_colors,
    label = FALSE,
    pt.size = 0.25,
    shuffle = TRUE,
    seed = 42,
    ncol = 3
  ) +
    ggtitle(
      "Final cell annotation by condition"
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      legend.title = element_blank()
    )

  print(p_clustering_condition)

  ggsave(
    filename = file.path(
      join_dir,
      "UMAP_Clustering_by_condition.pdf"
    ),
    plot = p_clustering_condition,
    width = 22,
    height = 8,
    limitsize = FALSE
  )

  p_mesenchymal_condition <- DimPlot(
    object = fibroblasts.Not.Integrated,
    reduction = "umap",
    group.by = "Clustering.Mesenchymal",
    split.by = "orig.ident",
    cols = mesenchymal_colors,
    label = FALSE,
    pt.size = 0.25,
    shuffle = TRUE,
    seed = 42,
    ncol = 3
  ) +
    ggtitle(
      "Mesenchymal-grouped annotation by condition"
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      legend.title = element_blank()
    )

  print(p_mesenchymal_condition)

  ggsave(
    filename = file.path(
      join_dir,
      "UMAP_Clustering_Mesenchymal_by_condition.pdf"
    ),
    plot = p_mesenchymal_condition,
    width = 22,
    height = 8,
    limitsize = FALSE
  )
}

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