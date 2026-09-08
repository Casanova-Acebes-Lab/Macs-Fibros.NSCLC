# ============================================================
# 15_SuppFigure4.R
# Supplementary Figure 4
#
# A. Global CD45+ UMAP
# B. CD45+ composition
# C. CD45+ marker bubbleplot
# D. Immune Response Potential Score (IRPS)
#
# Input:
#   results/objects/fibroblasts.Not.Integrated.annotated.rds
#
# Final cell annotation:
#   Clustering
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})

set.seed(42)

# ============================================================
# Paths
# ============================================================

input_file <- file.path(
  "results",
  "objects",
  "fibroblasts.Not.Integrated.annotated.rds"
)

plot_dir <- file.path(
  "results",
  "plots",
  "SuppFigure4"
)

table_dir <- file.path(
  "results",
  "tables",
  "SuppFigure4"
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop(
    "Input object not found: ",
    input_file
  )
}

# ============================================================
# Load final object
# ============================================================

main_obj <- readRDS(
  input_file
)

stopifnot(
  inherits(
    main_obj,
    "Seurat"
  )
)

required_metadata <- c(
  "orig.ident",
  "group",
  "Clustering"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(
    main_obj@meta.data
  )
)

if (length(missing_metadata) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(
      missing_metadata,
      collapse = ", "
    )
  )
}

if (!"umap" %in% Reductions(main_obj)) {
  stop(
    "Global UMAP reduction 'umap' not found."
  )
}

cat(
  "\nLoaded final object:\n",
  input_file,
  "\nCells: ",
  ncol(main_obj),
  "\nGenes: ",
  nrow(main_obj),
  "\n",
  sep = ""
)

cat(
  "\nFinal Clustering:\n"
)

print(
  table(
    main_obj$Clustering
  )
)

# ============================================================
# Nora palette
#
# Colors preserved.
# Names adapted to the exact labels in Clustering.
#
# ILC is absent because it is not present in the final object.
# ============================================================

Nora_colors2 <- c(

  # ----------------
  # Myeloid
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
  # CD8 / T cells
  # ----------------

  "Immature T Cells" = "#FFC300",
  "Stem like T Cells" = "#00CED1",
  "Tpex Transitional" = "#FF5A36",
  "CD8 Memory" = "#C95F82",
  "Cd8 Exhausted" = "#6A0DAD",

  # ----------------
  # Innate-like
  # ----------------

  "NK" = "#1F77B4",
  "NKT" = "#003566",
  "Tgd" = "#80ED99",

  # ----------------
  # CD45-negative
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
  # Fibroblasts
  # ----------------

  "Fibroblasts Universal" = "#FBCDC7",
  "Fibroblasts Activated 1" = "#F94144",
  "Fibroblasts Activated 2" = "#6D1A1A"
)

# ============================================================
# Immune cluster order
#
# Historical Nora_colors2 order, excluding ILC.
# Used in the UMAP legend and stacked composition.
# ============================================================

immune_cluster_order <- c(

  "Monocytes Ly6cLo",
  "Monocytes Ly6cHi",
  "MoMacs",
  "em-AMs",
  "mo-AMs",
  "mregDC",
  "DC1",
  "DC2",
  "pDC",
  "Neutrophils",
  "B cells",

  "CD4 Early Activated",
  "CD4 Naive",
  "Cd4 IL7R Memory Like",
  "Treg",

  "Immature T Cells",
  "Stem like T Cells",
  "Tpex Transitional",
  "CD8 Memory",
  "Cd8 Exhausted",

  "NK",
  "NKT",
  "Tgd"
)

# ============================================================
# Validate cluster names
# ============================================================

clustering_named <- as.character(
  main_obj$Clustering
)

names(
  clustering_named
) <- colnames(
  main_obj
)

cluster_names_present <- unique(
  clustering_named
)

missing_immune_clusters <- setdiff(
  immune_cluster_order,
  cluster_names_present
)

if (length(missing_immune_clusters) > 0) {
  stop(
    "Expected immune populations missing from Clustering: ",
    paste(
      missing_immune_clusters,
      collapse = ", "
    )
  )
}

missing_colors <- setdiff(
  immune_cluster_order,
  names(Nora_colors2)
)

if (length(missing_colors) > 0) {
  stop(
    "Missing colors for: ",
    paste(
      missing_colors,
      collapse = ", "
    )
  )
}

cat(
  "\nAll ",
  length(immune_cluster_order),
  " immune populations are present.\n",
  sep = ""
)

# ============================================================
# Select ALL CD45+ cells
#
# Selection is based on exact final Clustering names.
# ============================================================

cd45pos_cells <- names(
  clustering_named
)[
  clustering_named %in%
    immune_cluster_order
]

if (length(cd45pos_cells) == 0) {
  stop(
    "No CD45+ cells found."
  )
}

cd45pos <- subset(
  main_obj,
  cells = cd45pos_cells
)

cd45pos$immune_annotation <- factor(
  clustering_named[
    colnames(cd45pos)
  ],
  levels = immune_cluster_order
)

cat(
  "\nCD45+ cells selected: ",
  ncol(cd45pos),
  "\n",
  sep = ""
)

cat(
  "\nCD45+ cluster counts:\n"
)

print(
  table(
    cd45pos$immune_annotation
  )
)

# ============================================================
# FIGURE S4A
# CD45+ global UMAP
# ============================================================

umap_all <- Embeddings(
  main_obj,
  reduction = "umap"
)

if (
  !all(
    colnames(cd45pos) %in%
      rownames(umap_all)
  )
) {
  stop(
    "Some CD45+ cells are missing from the global UMAP."
  )
}

umap_cd45 <- umap_all[
  colnames(cd45pos),
  ,
  drop = FALSE
]

umap_table <- data.frame(

  cell =
    rownames(
      umap_cd45
    ),

  UMAP_1 =
    umap_cd45[, 1],

  UMAP_2 =
    umap_cd45[, 2],

  annotation =
    factor(
      as.character(
        cd45pos$immune_annotation
      ),
      levels =
        immune_cluster_order
    ),

  orig.ident =
    as.character(
      cd45pos$orig.ident
    ),

  group =
    as.character(
      cd45pos$group
    ),

  stringsAsFactors =
    FALSE
)

# ============================================================
# Visual cleaning of panel A
#
# Remove ONLY cells satisfying BOTH:
#
#   UMAP_2 > 0
#   UMAP_1 < 5
#
# This filtering is purely graphical and affects ONLY S4A.
#
# Complete CD45+ populations are retained for S4B and S4C.
# S4D uses the complete main object.
# ============================================================

umap_remove <- (
  umap_table$UMAP_2 > 0 &
    umap_table$UMAP_1 < 5
)

umap_table$plot_in_final <- !umap_remove

umap_plot_table <- umap_table[
  !umap_remove,
  ,
  drop = FALSE
]

excluded_umap_cells <- umap_table[
  umap_remove,
  ,
  drop = FALSE
]

cat(
  "\nUMAP visual filtering:\n",
  "Rule: UMAP_2 > 0 AND UMAP_1 < 5\n",
  "Cells before filtering: ",
  nrow(umap_table),
  "\n",
  "Cells excluded ONLY from panel A: ",
  sum(umap_remove),
  "\n",
  "Cells displayed in panel A: ",
  nrow(umap_plot_table),
  "\n",
  sep = ""
)

# ============================================================
# Export excluded cells
# ============================================================

excluded_umap_export <- excluded_umap_cells

excluded_umap_export$annotation <- as.character(
  excluded_umap_export$annotation
)

write.table(
  excluded_umap_export,
  file = file.path(
    table_dir,
    "SuppFigure4A_excluded_UMAP_cells.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Final S4A UMAP
#
# Small cell points.
# Enlarged legend symbols.
# Coordinates visible.
# ============================================================

p_s4a <- ggplot(
  umap_plot_table,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = annotation
  )
) +
  geom_point(
    size = 0.35,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = Nora_colors2[
      immune_cluster_order
    ],
    breaks = immune_cluster_order,
    limits = immune_cluster_order,
    drop = FALSE,
    name = NULL
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  ) +
  coord_equal() +
  labs(
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    legend.position = "right"
  )

ggsave(
  filename = file.path(
    plot_dir,
    "SuppFigure4A_CD45pos_UMAP.pdf"
  ),
  plot = p_s4a,
  width = 11,
  height = 9
)

# ============================================================
# S4A source data
#
# Contains all CD45+ cells.
# plot_in_final indicates which cells are displayed.
# ============================================================

umap_export <- umap_table

umap_export$annotation <- as.character(
  umap_export$annotation
)

write.table(
  umap_export,
  file = file.path(
    table_dir,
    "SuppFigure4A_source_data.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nFigure S4A completed.\n"
)

# ============================================================
# FIGURE S4B
# CD45+ stacked composition
#
# Uses COMPLETE cd45pos.
# UMAP visual filtering has no effect here.
#
# Cluster order follows Nora_colors2.
# ============================================================

composition_metadata <- data.frame(

  genotype =
    factor(
      as.character(
        cd45pos$orig.ident
      ),
      levels = c(
        "WT",
        "CD169DTR",
        "CCR2KO"
      )
    ),

  condition =
    factor(
      as.character(
        cd45pos$group
      ),
      levels = c(
        "Naive",
        "Tumor"
      )
    ),

  annotation =
    factor(
      as.character(
        cd45pos$immune_annotation
      ),
      levels =
        immune_cluster_order
    )
)

composition_table <- as.data.frame(
  table(
    composition_metadata$genotype,
    composition_metadata$condition,
    composition_metadata$annotation
  ),
  stringsAsFactors = FALSE
)

colnames(
  composition_table
) <- c(
  "genotype",
  "condition",
  "annotation",
  "n"
)

composition_table$genotype <- factor(
  composition_table$genotype,
  levels = c(
    "WT",
    "CD169DTR",
    "CCR2KO"
  )
)

composition_table$condition <- factor(
  composition_table$condition,
  levels = c(
    "Naive",
    "Tumor"
  )
)

composition_table$annotation <- factor(
  composition_table$annotation,
  levels = immune_cluster_order
)

composition_table <- composition_table %>%
  group_by(
    genotype,
    condition
  ) %>%
  mutate(
    total = sum(n),
    Proportion = n / total
  ) %>%
  ungroup()

composition_check <- composition_table %>%
  group_by(
    genotype,
    condition
  ) %>%
  summarise(
    sum_proportion =
      sum(
        Proportion
      ),
    .groups = "drop"
  )

cat(
  "\nComposition checks:\n"
)

print(
  composition_check
)

if (
  any(
    abs(
      composition_check$sum_proportion -
        1
    ) > 1e-8
  )
) {
  stop(
    "At least one composition bar does not sum to 1."
  )
}

p_s4b <- ggplot(
  composition_table,
  aes(
    x = condition,
    y = Proportion,
    fill = annotation
  )
) +
  geom_col(
    width = 0.8
  ) +
  facet_grid(
    . ~ genotype
  ) +
  scale_fill_manual(
    values = Nora_colors2[
      immune_cluster_order
    ],
    breaks = immune_cluster_order,
    limits = immune_cluster_order,
    drop = FALSE,
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = c(
      0,
      0.25,
      0.50,
      0.75,
      1
    ),
    expand = c(
      0,
      0
    )
  ) +
  labs(
    x = NULL,
    y = "Proportion"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    strip.background =
      element_blank(),

    strip.text =
      element_text(
        face = "bold"
      ),

    axis.text.x =
      element_text(
        angle = 0,
        hjust = 0.5
      ),

    legend.position =
      "right"
  )

ggsave(
  filename = file.path(
    plot_dir,
    "SuppFigure4B_CD45pos_Composition.pdf"
  ),
  plot = p_s4b,
  width = 10,
  height = 7
)

composition_export <- composition_table %>%
  mutate(

    genotype =
      as.character(
        genotype
      ),

    condition =
      as.character(
        condition
      ),

    annotation =
      as.character(
        annotation
      )
  )

write.table(
  composition_export,
  file = file.path(
    table_dir,
    "SuppFigure4B_source_data.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nFigure S4B completed.\n"
)

# ============================================================
# FIGURE S4C
# CD45+ marker bubbleplot
#
# Uses COMPLETE cd45pos.
#
# Historical marker definition retained.
# ============================================================

genes_manual_cd45pos <- list(

  "B cells" =
    c(
      "Cd19",
      "Cd79a"
    ),

  "CD4 Early Activated" =
    c(
      "Itm2a",
      "Nab2"
    ),

  "CD4 Naive" =
    c(
      "Igfbp4",
      "Ccr7"
    ),

  "Cd4 IL7R Memory Like" =
    c(
      "Il7r",
      "Maf"
    ),

  "CD8 Memory" =
    c(
      "Ctsw",
      "Adgrg5"
    ),

  "Cd8 Exhausted" =
    c(
      "Pdcd1",
      "Lag3"
    ),

  "DC1" =
    c(
      "Itgae",
      "Naaa"
    ),

  "DC2" =
    c(
      "Plbd1",
      "Clec9a"
    ),

  "em-AMs" =
    c(
      "Mrc1",
      "Marco",
      "Cd68"
    ),

  "Immature T Cells" =
    c(
      "Rag1",
      "Cd4",
      "Cd8a"
    ),

  "mo-AMs" =
    c(
      "Plet1",
      "Mpeg1"
    ),

  "MoMacs" =
    c(
      "C1qb",
      "Lyz1"
    ),

  "Monocytes Ly6cHi" =
    c(
      "Ly6c2",
      "Ccr2"
    ),

  "Monocytes Ly6cLo" =
    c(
      "Ace",
      "Adgre4"
    ),

  "mregDC" =
    c(
      "Fscn1",
      "Mreg"
    ),

  "Neutrophils" =
    c(
      "Csf3r",
      "S100a8"
    ),

  "NK" =
    c(
      "Nkg7",
      "Gzma"
    ),

  "NKT" =
    c(
      "Zfp683",
      "Zbtb16"
    ),

  "pDC" =
    c(
      "Siglech",
      "Bst2"
    ),

  "Stem like T Cells" =
    c(
      "Top2a",
      "Mki67"
    ),

  "Tgd" =
    c(
      "Myl10",
      "Tcrg-C1"
    ),

  "Tpex Transitional" =
    c(
      "Id3",
      "Cd28"
    ),

  "Treg" =
    c(
      "Foxp3",
      "Ctla4"
    )
)

bubble_cluster_order <- names(
  genes_manual_cd45pos
)

if (
  !setequal(
    bubble_cluster_order,
    immune_cluster_order
  )
) {
  stop(
    "Bubbleplot and CD45+ annotations do not contain the same populations."
  )
}

genes_plot <- unique(
  unlist(
    genes_manual_cd45pos,
    use.names = FALSE
  )
)

obj_dot <- cd45pos

DefaultAssay(
  obj_dot
) <- "RNA"

obj_dot <- tryCatch(
  {
    JoinLayers(
      object = obj_dot,
      assay = "RNA"
    )
  },
  error = function(e) {

    message(
      "JoinLayers not required/possible: ",
      conditionMessage(e)
    )

    obj_dot
  }
)

if (
  !"data" %in%
    Layers(
      obj_dot[["RNA"]]
    )
) {
  stop(
    "Normalized RNA data layer is missing for DotPlot."
  )
}

genes_present <- genes_plot[
  genes_plot %in%
    rownames(
      obj_dot
    )
]

genes_missing <- setdiff(
  genes_plot,
  genes_present
)

cat(
  "\nBubbleplot genes present: ",
  length(genes_present),
  "/",
  length(genes_plot),
  "\n",
  sep = ""
)

cat(
  "Bubbleplot genes missing: ",
  ifelse(
    length(genes_missing) == 0,
    "none",
    paste(
      genes_missing,
      collapse = ", "
    )
  ),
  "\n",
  sep = ""
)

if (length(genes_present) == 0) {
  stop(
    "None of the bubbleplot genes are present."
  )
}

obj_dot$immune_annotation <- factor(
  as.character(
    obj_dot$immune_annotation
  ),
  levels =
    bubble_cluster_order
)

Idents(
  obj_dot
) <- "immune_annotation"

dot_raw <- DotPlot(
  object = obj_dot,
  features = genes_present,
  group.by = "immune_annotation",
  assay = "RNA",
  scale = FALSE
)

dot_table <- dot_raw$data %>%
  transmute(

    Gene =
      as.character(
        features.plot
      ),

    annotation =
      as.character(
        id
      ),

    AvgExpression =
      avg.exp,

    PctExpressed =
      pct.exp
  )

color_cap <- as.numeric(
  quantile(
    dot_table$AvgExpression,
    probs = 0.95,
    na.rm = TRUE
  )
)

dot_table$AvgExpression_plot <- pmin(
  dot_table$AvgExpression,
  color_cap
)

dot_table$Gene <- factor(
  dot_table$Gene,
  levels =
    genes_present
)

dot_table$annotation <- factor(
  dot_table$annotation,
  levels =
    rev(
      bubble_cluster_order
    )
)

bubble_cols <- colorRampPalette(
  rev(
    RColorBrewer::brewer.pal(
      n = 7,
      name = "RdYlBu"
    )
  )
)(100)

p_s4c <- ggplot(
  dot_table,
  aes(
    x = Gene,
    y = annotation
  )
) +
  geom_point(
    aes(
      size = PctExpressed,
      color = AvgExpression_plot
    ),
    alpha = 1
  ) +
  scale_size(
    range = c(
      0.2,
      7
    ),
    limits = c(
      0,
      100
    ),
    breaks = c(
      0,
      25,
      50,
      75,
      100
    ),
    name =
      "Percentage Expressed"
  ) +
  scale_color_gradientn(
    colors =
      bubble_cols,
    name =
      "Average Expression"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(
    base_size = 13
  ) +
  theme(

    panel.grid.major =
      element_line(
        color = "grey90",
        linewidth = 0.45
      ),

    panel.grid.minor =
      element_blank(),

    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        size = 7
      ),

    axis.text.y =
      element_text(
        size = 8
      ),

    legend.title =
      element_text(
        size = 10
      ),

    legend.text =
      element_text(
        size = 9
      )
  )

ggsave(
  filename = file.path(
    plot_dir,
    "SuppFigure4C_CD45pos_Bubbleplot.pdf"
  ),
  plot = p_s4c,
  width = 18,
  height = 9
)

dot_export <- dot_table %>%
  mutate(

    Gene =
      as.character(
        Gene
      ),

    annotation =
      as.character(
        annotation
      ),

    color_cap =
      color_cap
  )

write.table(
  dot_export,
  file = file.path(
    table_dir,
    "SuppFigure4C_source_data.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nFigure S4C completed.\n"
)

# ============================================================
# FIGURE S4D
# Immune Response Potential Score (IRPS)
#
# Uses FULL main_obj.
#
# Required identities are already present in Clustering:
#
#   Tpex Transitional
#   Cd8 Exhausted
#   Treg
#   Fibroblasts Universal
#   Fibroblasts Activated 1
#   Fibroblasts Activated 2
# ============================================================

data_irps <- main_obj

data_irps$irps_identity <- as.character(
  data_irps$Clustering
)

irps_required <- c(
  "Tpex Transitional",
  "Cd8 Exhausted",
  "Fibroblasts Universal",
  "Fibroblasts Activated 1",
  "Fibroblasts Activated 2",
  "Treg"
)

irps_counts <- table(
  factor(
    data_irps$irps_identity,
    levels =
      irps_required
  )
)

cat(
  "\nIRPS population counts:\n"
)

print(
  irps_counts
)

missing_irps_populations <- names(
  irps_counts
)[
  irps_counts == 0
]

if (
  length(
    missing_irps_populations
  ) > 0
) {
  stop(
    "Missing IRPS populations: ",
    paste(
      missing_irps_populations,
      collapse = ", "
    )
  )
}

# ============================================================
# Condition groups
# ============================================================

data_irps$Condition_Group <- paste(
  as.character(
    data_irps$orig.ident
  ),
  as.character(
    data_irps$group
  ),
  sep = "_"
)

condition_groups <- c(
  "WT_Naive",
  "WT_Tumor",
  "CD169DTR_Naive",
  "CD169DTR_Tumor",
  "CCR2KO_Naive",
  "CCR2KO_Tumor"
)

missing_conditions <- setdiff(
  condition_groups,
  unique(
    data_irps$Condition_Group
  )
)

if (
  length(
    missing_conditions
  ) > 0
) {
  stop(
    "Missing condition groups: ",
    paste(
      missing_conditions,
      collapse = ", "
    )
  )
}

# ============================================================
# RNA expression
# ============================================================

DefaultAssay(
  data_irps
) <- "RNA"

data_irps <- tryCatch(
  {
    JoinLayers(
      object = data_irps,
      assay = "RNA"
    )
  },
  error = function(e) {

    message(
      "JoinLayers not required/possible: ",
      conditionMessage(e)
    )

    data_irps
  }
)

if (
  !"data" %in%
    Layers(
      data_irps[["RNA"]]
    )
) {
  stop(
    "Normalized RNA data layer is missing."
  )
}

expr_data <- GetAssayData(
  object = data_irps,
  assay = "RNA",
  layer = "data"
)

# ============================================================
# Mean marker expression per population
# ============================================================

mean_expr_per_cell <- function(
    genes,
    cells
) {

  if (
    length(
      cells
    ) == 0
  ) {
    return(0)
  }

  genes <- intersect(
    genes,
    rownames(
      expr_data
    )
  )

  if (
    length(
      genes
    ) == 0
  ) {
    return(0)
  }

  mat <- expr_data[
    genes,
    cells,
    drop = FALSE
  ]

  cell_means <- Matrix::colMeans(
    mat
  )

  mean(
    cell_means
  )
}

# ============================================================
# IRPS gene sets
# ============================================================

Tpex_genes <- c(
  "Id3",
  "Cd28",
  "Lef1"
)

Exh_genes <- c(
  "Pdcd1",
  "Lag3",
  "Havcr2"
)

FibUniv_genes <- c(
  "Col1a2",
  "Scara5",
  "Pi16",
  "Clec3b",
  "Pcolce"
)

FibA1_genes <- c(
  "Hhip",
  "Grem2",
  "Bmp5",
  "Tagln"
)

FibA2_genes <- c(
  "Ddit4",
  "Klf15",
  "Fmo2",
  "Acta2"
)

Treg_genes <- c(
  "Foxp3",
  "Ctla4",
  "Maf"
)

eps <- 1e-6

all_irps_genes <- unique(
  c(
    Tpex_genes,
    Exh_genes,
    FibUniv_genes,
    FibA1_genes,
    FibA2_genes,
    Treg_genes
  )
)

missing_irps_genes <- setdiff(
  all_irps_genes,
  rownames(
    expr_data
  )
)

cat(
  "\nMissing IRPS genes: ",
  ifelse(
    length(
      missing_irps_genes
    ) == 0,
    "none",
    paste(
      missing_irps_genes,
      collapse = ", "
    )
  ),
  "\n",
  sep = ""
)

# ============================================================
# Named vectors
# ============================================================

irps_identity_named <- as.character(
  data_irps$irps_identity
)

names(
  irps_identity_named
) <- colnames(
  data_irps
)

condition_named <- as.character(
  data_irps$Condition_Group
)

names(
  condition_named
) <- colnames(
  data_irps
)

# ============================================================
# Calculate raw IRPS components
# ============================================================

irps_results <- lapply(
  condition_groups,
  function(
      condition_group
  ) {

    condition_cells <- names(
      condition_named
    )[
      condition_named ==
        condition_group
    ]

    n_total <- length(
      condition_cells
    )

    if (
      n_total == 0
    ) {
      stop(
        "No cells for ",
        condition_group
      )
    }

    get_population_cells <- function(
        population
    ) {

      condition_cells[
        irps_identity_named[
          condition_cells
        ] == population
      ]
    }

    cells_Tpex <-
      get_population_cells(
        "Tpex Transitional"
      )

    cells_Exh <-
      get_population_cells(
        "Cd8 Exhausted"
      )

    cells_FibUniv <-
      get_population_cells(
        "Fibroblasts Universal"
      )

    cells_FibA1 <-
      get_population_cells(
        "Fibroblasts Activated 1"
      )

    cells_FibA2 <-
      get_population_cells(
        "Fibroblasts Activated 2"
      )

    cells_Treg <-
      get_population_cells(
        "Treg"
      )

    prop_Tpex <-
      length(
        cells_Tpex
      ) /
      n_total

    prop_Exh <-
      length(
        cells_Exh
      ) /
      n_total

    prop_FibUniv <-
      length(
        cells_FibUniv
      ) /
      n_total

    prop_FibA1 <-
      length(
        cells_FibA1
      ) /
      n_total

    prop_FibA2 <-
      length(
        cells_FibA2
      ) /
      n_total

    prop_Treg <-
      length(
        cells_Treg
      ) /
      n_total

    expr_Tpex <-
      mean_expr_per_cell(
        Tpex_genes,
        cells_Tpex
      )

    expr_Exh <-
      mean_expr_per_cell(
        Exh_genes,
        cells_Exh
      )

    expr_FibUniv <-
      mean_expr_per_cell(
        FibUniv_genes,
        cells_FibUniv
      )

    expr_FibA1 <-
      mean_expr_per_cell(
        FibA1_genes,
        cells_FibA1
      )

    expr_FibA2 <-
      mean_expr_per_cell(
        FibA2_genes,
        cells_FibA2
      )

    expr_Treg <-
      mean_expr_per_cell(
        Treg_genes,
        cells_Treg
      )

    data.frame(

      Condition_Group =
        condition_group,

      n_total =
        n_total,

      n_Tpex =
        length(
          cells_Tpex
        ),

      n_Exh =
        length(
          cells_Exh
        ),

      n_FibUniv =
        length(
          cells_FibUniv
        ),

      n_FibA1 =
        length(
          cells_FibA1
        ),

      n_FibA2 =
        length(
          cells_FibA2
        ),

      n_Treg =
        length(
          cells_Treg
        ),

      prop_Tpex =
        prop_Tpex,

      prop_Exh =
        prop_Exh,

      prop_FibUniv =
        prop_FibUniv,

      prop_FibA1 =
        prop_FibA1,

      prop_FibA2 =
        prop_FibA2,

      prop_Treg =
        prop_Treg,

      expr_Tpex =
        expr_Tpex,

      expr_Exh =
        expr_Exh,

      expr_FibUniv =
        expr_FibUniv,

      expr_FibA1 =
        expr_FibA1,

      expr_FibA2 =
        expr_FibA2,

      expr_Treg =
        expr_Treg,

      stringsAsFactors =
        FALSE
    )
  }
)

IRPS_df <- bind_rows(
  irps_results
)

# ============================================================
# Historical abundance / expression definitions
# ============================================================

IRPS_df <- IRPS_df %>%
  mutate(

    Tpex_abundance =
      prop_Tpex /
      (
        prop_Exh +
          eps
      ),

    Fibro_abundance =
      prop_FibUniv /
      (
        prop_FibA1 +
          prop_FibA2 +
          eps
      ),

    Treg_abundance =
      -prop_Treg,

    Tpex_expression =
      log2(
        (
          expr_Tpex +
            eps
        ) /
          (
            expr_Exh +
              eps
          )
      ),

    Fibro_expression =
      log2(
        (
          expr_FibUniv +
            eps
        ) /
          (
            expr_FibA1 +
              expr_FibA2 +
              eps
          )
      ),

    Treg_expression =
      -log2(
        expr_Treg +
          eps
      )
  )

# ============================================================
# Historical normalization
#
# Deliberately reproduce:
#
#   (x - min(x)) / max(x)
#
# NOT standard min-max normalization.
# ============================================================

normalize <- function(x) {

  if (
    length(
      unique(
        x
      )
    ) == 1
  ) {

    return(
      rep(
        0.5,
        length(
          x
        )
      )
    )
  }

  (
    x -
      min(
        x
      )
  ) /
    max(
      x
    )
}

IRPS_df <- IRPS_df %>%
  mutate(

    Tpex_abundance_norm =
      normalize(
        Tpex_abundance
      ),

    Fibro_abundance_norm =
      normalize(
        Fibro_abundance
      ),

    Treg_abundance_norm =
      -normalize(
        -Treg_abundance
      ),

    Tpex_expr_norm =
      normalize(
        Tpex_expression
      ),

    Fibro_expr_norm =
      normalize(
        Fibro_expression
      ),

    Treg_expr_norm =
      -normalize(
        -Treg_expression
      )
  )

# ============================================================
# Historical weights
# ============================================================

weight_tpex <- 0.3
weight_fibro <- 0.4
weight_treg <- 0.25

IRPS_df <- IRPS_df %>%
  mutate(

    Abundance_score =
      Tpex_abundance_norm *
        weight_tpex +
      Fibro_abundance_norm *
        weight_fibro +
      Treg_abundance_norm *
        weight_treg,

    Expression_score =
      Tpex_expr_norm *
        weight_tpex +
      Fibro_expr_norm *
        weight_fibro +
      Treg_expr_norm *
        weight_treg,

    Combined_score =
      sign(
        Abundance_score +
          Expression_score
      ) *
      sqrt(
        abs(
          Abundance_score *
            Expression_score
        )
      )
  )

# ============================================================
# IRPS plot data
# ============================================================

condition_order_plot <- c(
  "WT Naive",
  "WT Tumor",
  "CD169DTR Naive",
  "CD169DTR Tumor",
  "CCR2KO Naive",
  "CCR2KO Tumor"
)

IRPS_df$Condition_Group_plot <- gsub(
  "_",
  " ",
  IRPS_df$Condition_Group
)

IRPS_df$Condition_Group_plot <- factor(
  IRPS_df$Condition_Group_plot,
  levels =
    condition_order_plot
)

IRPS_long <- IRPS_df %>%
  select(
    Condition_Group_plot,
    Abundance_score,
    Expression_score,
    Combined_score
  ) %>%
  pivot_longer(
    cols = c(
      Abundance_score,
      Expression_score,
      Combined_score
    ),
    names_to =
      "Score_type",
    values_to =
      "Score"
  ) %>%
  mutate(

    Score_type =
      recode(
        Score_type,

        "Abundance_score" =
          "Cell Ratio",

        "Expression_score" =
          "Markers Average Expression",

        "Combined_score" =
          "Combined"
      ),

    Score_type =
      factor(
        Score_type,
        levels = c(
          "Cell Ratio",
          "Markers Average Expression",
          "Combined"
        )
      ),

    Condition_Group =
      factor(
        as.character(
          Condition_Group_plot
        ),
        levels =
          condition_order_plot
      ),

    Genotype =
      case_when(

        grepl(
          "^WT ",
          as.character(
            Condition_Group
          )
        ) ~
          "WT",

        grepl(
          "^CD169DTR ",
          as.character(
            Condition_Group
          )
        ) ~
          "CD169DTR",

        grepl(
          "^CCR2KO ",
          as.character(
            Condition_Group
          )
        ) ~
          "CCR2KO",

        TRUE ~
          NA_character_
      )
  )

geno_colors <- c(
  "WT" = "#A6CEE3",
  "CD169DTR" = "#1F78B4",
  "CCR2KO" = "#08306B"
)

# ============================================================
# S4D IRPS lollipop
# ============================================================

p_s4d <- ggplot(
  IRPS_long,
  aes(
    x = Score,
    y = Condition_Group,
    color = Genotype
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = Score,
      y = Condition_Group,
      yend = Condition_Group
    ),
    linewidth = 0.7
  ) +
  geom_point(
    size = 4
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  facet_wrap(
    ~Score_type,
    ncol = 1,
    scales = "fixed"
  ) +
  scale_color_manual(
    values =
      geno_colors
  ) +
  coord_cartesian(
    xlim = c(
      -1,
      1
    )
  ) +
  scale_y_discrete(
    limits =
      rev(
        condition_order_plot
      )
  ) +
  theme_classic(
    base_size = 14
  ) +
  theme(

    legend.position =
      "none",

    strip.text =
      element_text(
        size = 14
      ),

    plot.title =
      element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      )
  ) +
  labs(
    title =
      "Immune Response Potential Score (IRPS)",
    x =
      "Score",
    y =
      NULL
  )

ggsave(
  filename = file.path(
    plot_dir,
    "SuppFigure4D_IRPS.pdf"
  ),
  plot = p_s4d,
  width = 8,
  height = 8
)

# ============================================================
# Export IRPS source data
# ============================================================

IRPS_export <- IRPS_df

IRPS_export$Condition_Group_plot <-
  as.character(
    IRPS_export$Condition_Group_plot
  )

write.table(
  IRPS_export,
  file = file.path(
    table_dir,
    "SuppFigure4D_IRPS_scores.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

IRPS_long_export <- IRPS_long

IRPS_long_export$Condition_Group_plot <-
  as.character(
    IRPS_long_export$Condition_Group_plot
  )

IRPS_long_export$Condition_Group <-
  as.character(
    IRPS_long_export$Condition_Group
  )

IRPS_long_export$Score_type <-
  as.character(
    IRPS_long_export$Score_type
  )

write.table(
  IRPS_long_export,
  file = file.path(
    table_dir,
    "SuppFigure4D_IRPS_plot_source.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nFigure S4D completed.\n"
)

# ============================================================
# Final summary
# ============================================================

cat(
  "\n",
  "============================================\n",
  "SUPPLEMENTARY FIGURE 4 COMPLETE\n",
  "============================================\n",

  "Input object: ",
  input_file,
  "\n",

  "Annotation: Clustering\n",

  "CD45+ clusters: ",
  length(
    immune_cluster_order
  ),
  "\n",

  "CD45+ cells: ",
  ncol(
    cd45pos
  ),
  "\n",

  "Stacked order: Nora_colors2 order\n",

  "ILC: absent from final Clustering and excluded\n",

  "UMAP reduction: umap\n",

  "UMAP visual filter: UMAP_2 > 0 AND UMAP_1 < 5\n",

  "UMAP cells excluded ONLY from S4A: ",
  sum(
    umap_remove
  ),
  "\n",

  "UMAP cells displayed: ",
  nrow(
    umap_plot_table
  ),
  "\n",

  "UMAP cell point size: 0.35\n",

  "UMAP legend point size: 3\n",

  "Missing bubble genes: ",
  ifelse(
    length(
      genes_missing
    ) == 0,
    "none",
    paste(
      genes_missing,
      collapse = ", "
    )
  ),
  "\n",

  "IRPS annotation: Clustering\n",

  "IRPS object: full main_obj\n",

  "============================================\n",
  sep = ""
)