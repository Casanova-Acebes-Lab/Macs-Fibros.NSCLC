# ============================================================
# 08_Figure1.R
# Figure 1: panels F–I
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)

source("R/Functions.R")

set.seed(1234)

# ============================================================
# Paths
# ============================================================

config <- read.delim(
  "config.tsv",
  stringsAsFactors = FALSE
)

config <- setNames(
  config$value,
  config$key
)

rsdir <- config[["rsdir"]]
outdir <- config[["outdir"]]

figdir <- file.path(
  outdir,
  "Figure1"
)

tabledir <- file.path(
  rsdir,
  "tables",
  "Figure1"
)

dir.create(
  figdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  tabledir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Load final annotated object
# ============================================================

fibroblasts <- readRDS(
  file.path(
    rsdir,
    "objects",
    "fibroblasts.Not.Integrated.annotated.rds"
  )
)

# ============================================================
# Figure 1F: general clustering
# ============================================================

Idents(fibroblasts) <- "Clustering.Mesenchymal"

clustering_levels <- unique(
  as.character(fibroblasts$Clustering.Mesenchymal)
)

clustering_colors <- Nora_colors[
  clustering_levels
]

Figure1F <- DimPlot(
  object = fibroblasts,
  reduction = "umap",
  group.by = "Clustering.Mesenchymal",
  cols = clustering_colors,
  label = FALSE,
  pt.size = 0.30,
  shuffle = TRUE,
  seed = 1234,
  raster = FALSE
) +
  labs(
    title = NULL,
    color = NULL
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 3)
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(
      size = 7,
      colour = "black"
    ),
    legend.key.height = grid::unit(0.35, "cm"),
    axis.title = element_text(
      size = 10,
      colour = "black"
    ),
    axis.text = element_text(
      size = 8,
      colour = "black"
    ),
    axis.line = element_line(
      linewidth = 0.4,
      colour = "black"
    )
  )

Figure1F

ggsave(
  filename = file.path(
    figdir,
    "Figure1F_general_clustering.pdf"
  ),
  plot = Figure1F,
  width = 11,
  height = 6.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

# ============================================================
# Figure 1F: FindAllMarkers table
# ============================================================

DefaultAssay(fibroblasts) <- "RNA"

fibroblasts <- JoinLayers(
  fibroblasts,
  assay = "RNA"
)

Idents(fibroblasts) <- "Clustering.Mesenchymal"

markers_Figure1F <- FindAllMarkers(
  object = fibroblasts,
  assay = "RNA",
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25
) %>%
  filter(
    p_val_adj < 0.05
  ) %>%
  arrange(
    cluster,
    desc(avg_log2FC)
  )

write.table(
  markers_Figure1F,
  file = file.path(
    tabledir,
    "Figure1F_FindAllMarkers_significant.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



# ============================================================
# Figure 1H: mesenchymal quantification
# ============================================================

fibro_colors <- c(
  "Mesothelial"             = "#F9C74F",
  "Fibroblasts Universal"   = "#FBCDC7",
  "Fibroblasts Activated 1" = "#F94144",
  "Fibroblasts Activated 2" = "#6D1A1A"
)

fibro_order <- c(
  "Mesothelial",
  "Fibroblasts Universal",
  "Fibroblasts Activated 1",
  "Fibroblasts Activated 2"
)

# Subset mesenquimal
mesenchymal <- subset(
  fibroblasts,
  subset = Clustering.Mesenchymal == "Mesenchymal"
)

# Tabla de frecuencias y proporciones
table_Figure1H <- mesenchymal[[]] %>%
  filter(Clustering %in% fibro_order) %>%
  count(
    orig.ident,
    group,
    Clustering,
    name = "Freq"
  ) %>%
  tidyr::complete(
    orig.ident,
    group,
    Clustering = fibro_order,
    fill = list(Freq = 0)
  ) %>%
  group_by(
    orig.ident,
    group
  ) %>%
  mutate(
    Proportion = Freq / sum(Freq)
  ) %>%
  ungroup() %>%
  rename(
    Group = orig.ident,
    Type = group,
    Cell_Types = Clustering
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = c("WT", "CD169DTR", "CCR2KO")
    ),
    Type = factor(
      Type,
      levels = c("Naive", "Tumor")
    ),
    Cell_Types = factor(
      Cell_Types,
      levels = fibro_order
    )
  ) %>%
  arrange(
    Group,
    Type,
    Cell_Types
  )

# ============================================================
# Plot
# ============================================================

Figure1H <- ggplot(
  table_Figure1H,
  aes(
    x = Type,
    y = Proportion,
    fill = Cell_Types
  )
) +
  geom_col(
    width = 0.92
  ) +
  facet_grid(
    . ~ Group
  ) +
  scale_fill_manual(
    values = fibro_colors,
    breaks = fibro_order,
    drop = FALSE
  ) +
  scale_x_discrete(
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = scales::number_format(accuracy = 0.01),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Proportion",
    fill = NULL,
    title = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_blank(),

    strip.background = element_rect(
      fill = "white",
      colour = "black",
      linewidth = 0.5
    ),
    strip.text = element_text(
      size = 9,
      colour = "black"
    ),

    panel.spacing.x = grid::unit(0.35, "cm"),

    axis.text.x = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 10,
      colour = "black"
    ),

    legend.title = element_blank(),
    legend.text = element_text(
      size = 8,
      colour = "black"
    ),
    legend.key.height = grid::unit(0.4, "cm")
  )

Figure1H

# ============================================================
# Save plot
# ============================================================

ggsave(
  filename = file.path(
    figdir,
    "Figure1H_fibroblast_quantification.pdf"
  ),
  plot = Figure1H,
  width = 9,
  height = 10,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

# ============================================================
# Save source table
# ============================================================

write.table(
  table_Figure1H,
  file = file.path(
    tabledir,
    "Figure1H_fibroblast_quantification.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)




# ============================================================
# Figure 1G: mesenchymal marker bubble plot
# ============================================================

gene_list_Figure1G <- c(
  "Dpt",
  "Pdgfra",
  "Pi16",
  "Postn",
  "Col1a1",
  "Hhip",
  "Grem2",
  "Tnc",
  "Myh11",
  "Acta2",
  "Bmp5",
  "Ddit4",
  "Cdo1",
  "Npr3",
  "Pdpn",
  "Upk3b"
)

# La anotación de los subtipos está en Clustering
Idents(mesenchymal) <- "Clustering"

# Extraer expresión normalizada
expression_Figure1G <- FetchData(
  object = mesenchymal,
  vars = gene_list_Figure1G,
  layer = "data"
)

expression_Figure1G$cluster <- as.character(
  Idents(mesenchymal)[rownames(expression_Figure1G)]
)

# Formato largo
long_Figure1G <- expression_Figure1G %>%
  tidyr::pivot_longer(
    cols = all_of(gene_list_Figure1G),
    names_to = "gene",
    values_to = "expression"
  )

# Expresión media y porcentaje de células positivas
table_Figure1G <- long_Figure1G %>%
  group_by(cluster, gene) %>%
  summarise(
    avg_expression = mean(expression),
    pct_expressed = mean(expression > 0) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    cluster = dplyr::recode(
      cluster,
      "Fibroblasts Universal"   = "F. Universal",
      "Fibroblasts Activated 1" = "F. Activated 1",
      "Fibroblasts Activated 2" = "F. Activated 2"
    )
  )

# Orden vertical:
# arriba Universal, luego Activated 1, Activated 2 y abajo Mesothelial
cluster_order_Figure1G <- c(
  "Mesothelial",
  "F. Activated 2",
  "F. Activated 1",
  "F. Universal"
)

table_Figure1G <- table_Figure1G %>%
  mutate(
    cluster = factor(
      cluster,
      levels = cluster_order_Figure1G
    ),
    gene = factor(
      gene,
      levels = gene_list_Figure1G
    ),

    # Mismo z-score global utilizado en la figura original
    zscore = (
      avg_expression -
        mean(avg_expression, na.rm = TRUE)
    ) /
      sd(avg_expression, na.rm = TRUE)
  ) %>%
  arrange(cluster, gene)

# ============================================================
# Plot
# ============================================================

Figure1G <- ggplot(
  table_Figure1G,
  aes(
    x = gene,
    y = cluster
  )
) +
  geom_point(
    aes(
      size = pct_expressed,
      color = zscore
    )
  ) +
  scale_size_continuous(
    name = "Percentage Expressed",
    breaks = c(0, 25, 50, 75),
    range = c(0.5, 8)
  ) +
  scale_color_gradient(
    name = "Z-Score",
    low = "blue",
    high = "red"
  ) +
  labs(
    x = "Gene",
    y = "Cluster",
    title = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_blank(),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 10,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 10,
      colour = "black"
    ),

    axis.title.x = element_text(
      size = 11,
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 11,
      colour = "black"
    ),

    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.5
    ),
    panel.grid.minor = element_blank(),

    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),

    legend.box.spacing = grid::unit(0.4, "cm")
  )

Figure1G

# ============================================================
# Save plot
# ============================================================

ggsave(
  filename = file.path(
    figdir,
    "Figure1G_mesenchymal_bubbleplot.pdf"
  ),
  plot = Figure1G,
  width = 10,
  height = 6,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

# ============================================================
# Save source table
# ============================================================

write.table(
  table_Figure1G,
  file = file.path(
    tabledir,
    "Figure1G_mesenchymal_bubbleplot.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)





# ============================================================
# Figure 1I
# Differential expression + GO Biological Process fGSEA
#
# Differential:
#   - Fibroblasts Universal
#   - Fibroblasts Activated 1
#   - Fibroblasts Activated 2
#   - Mesothelial
#
# fGSEA:
#   - Fibroblasts Universal
#   - Fibroblasts Activated 1
#   - Fibroblasts Activated 2
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(msigdbr)
library(fgsea)

set.seed(1234)

# ============================================================
# 0. Directorios y funciones auxiliares
# ============================================================

tabledir_Figure1I <- file.path(
  tabledir,
  "Figure1I"
)

dir.create(
  tabledir_Figure1I,
  recursive = TRUE,
  showWarnings = FALSE
)

write_tsv_base <- function(data, filename) {

  write.table(
    data,
    file = filename,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

format_fgsea_table <- function(data) {

  output <- as.data.frame(data)

  if ("leadingEdge" %in% colnames(output)) {

    output$leadingEdge <- vapply(
      output$leadingEdge,
      function(genes) {

        genes <- unlist(
          genes,
          use.names = FALSE
        )

        genes <- genes[
          !is.na(genes) &
          genes != ""
        ]

        if (length(genes) == 0) {
          return("")
        }

        paste(
          genes,
          collapse = ";"
        )
      },
      character(1)
    )
  }

  output
}

normalise_pathway_name <- function(x) {

  x <- toupper(x)

  x <- gsub(
    "^GOBP_",
    "",
    x
  )

  x <- gsub(
    "_",
    " ",
    x
  )

  x <- gsub(
    "[[:space:]]+",
    " ",
    x
  )

  trimws(x)
}

normalise_pathway_key <- function(x) {

  x <- normalise_pathway_name(x)

  gsub(
    "[^A-Z0-9]",
    "",
    x
  )
}

# ============================================================
# 1. Definir estados
# ============================================================

# Estados utilizados en el diferencial
cluster_order_Figure1I <- c(
  "Fibroblasts Universal",
  "Fibroblasts Activated 1",
  "Fibroblasts Activated 2",
  "Mesothelial"
)

# Estados utilizados en fGSEA
fgsea_cluster_order_Figure1I <- c(
  "Fibroblasts Universal",
  "Fibroblasts Activated 1",
  "Fibroblasts Activated 2"
)

cluster_labels_Figure1I <- c(
  "Fibroblasts Universal"   = "F. Universal",
  "Fibroblasts Activated 1" = "F. Activated 1",
  "Fibroblasts Activated 2" = "F. Activated 2",
  "Mesothelial"             = "Mesothelial"
)

cluster_files_Figure1I <- c(
  "Fibroblasts Universal"   = "Universal",
  "Fibroblasts Activated 1" = "Activated1",
  "Fibroblasts Activated 2" = "Activated2",
  "Mesothelial"             = "Mesothelial"
)

# ============================================================
# 2. Preparar el subset mesenquimal
# ============================================================

mesenchymal_gsea <- subset(
  fibroblasts,
  subset = Clustering.Mesenchymal == "Mesenchymal"
)

mesenchymal_gsea <- subset(
  mesenchymal_gsea,
  subset = Clustering %in% cluster_order_Figure1I
)

DefaultAssay(mesenchymal_gsea) <- "RNA"

mesenchymal_gsea <- JoinLayers(
  mesenchymal_gsea,
  assay = "RNA"
)

mesenchymal_gsea$Clustering <- factor(
  as.character(
    mesenchymal_gsea$Clustering
  ),
  levels = cluster_order_Figure1I
)

Idents(mesenchymal_gsea) <- "Clustering"

# ============================================================
# 3. Diferencial
# Cada estado frente al resto del subset mesenquimal
# ============================================================

markers_Figure1I <- FindAllMarkers(
  object = mesenchymal_gsea,
  assay = "RNA",
  test.use = "wilcox",
  only.pos = FALSE,
  min.pct = 0,
  logfc.threshold = 0,
  return.thresh = 1.01
)

# Compatibilidad entre versiones de Seurat
if (
  !"avg_log2FC" %in% colnames(markers_Figure1I) &&
  "avg_logFC" %in% colnames(markers_Figure1I)
) {

  markers_Figure1I <- markers_Figure1I %>%
    rename(
      avg_log2FC = avg_logFC
    )
}

markers_Figure1I <- markers_Figure1I %>%
  mutate(
    cluster = factor(
      as.character(cluster),
      levels = cluster_order_Figure1I
    )
  ) %>%
  arrange(
    cluster,
    p_val_adj,
    desc(abs(avg_log2FC))
  )

# ============================================================
# 4. Separar diferenciales
# ============================================================

DEG.FBU <- markers_Figure1I %>%
  filter(
    cluster == "Fibroblasts Universal"
  )

DEG.FB1 <- markers_Figure1I %>%
  filter(
    cluster == "Fibroblasts Activated 1"
  )

DEG.FB2 <- markers_Figure1I %>%
  filter(
    cluster == "Fibroblasts Activated 2"
  )

DEG.Mesothelial <- markers_Figure1I %>%
  filter(
    cluster == "Mesothelial"
  )

markers_top50_Figure1I <- markers_Figure1I %>%
  group_by(cluster) %>%
  arrange(
    p_val_adj,
    desc(abs(avg_log2FC)),
    .by_group = TRUE
  ) %>%
  slice_head(n = 50) %>%
  ungroup()

table.ISO.CD45.fibros2 <- markers_top50_Figure1I

# ============================================================
# 5. Guardar tablas del diferencial
# ============================================================

write_tsv_base(
  markers_Figure1I,
  file.path(
    tabledir_Figure1I,
    "Figure1I_FindAllMarkers_all.tsv"
  )
)

write_tsv_base(
  markers_top50_Figure1I,
  file.path(
    tabledir_Figure1I,
    "Figure1I_FindAllMarkers_top50_per_cluster.tsv"
  )
)

write_tsv_base(
  DEG.FBU,
  file.path(
    tabledir_Figure1I,
    "Figure1I_DEG_Universal.tsv"
  )
)

write_tsv_base(
  DEG.FB1,
  file.path(
    tabledir_Figure1I,
    "Figure1I_DEG_Activated1.tsv"
  )
)

write_tsv_base(
  DEG.FB2,
  file.path(
    tabledir_Figure1I,
    "Figure1I_DEG_Activated2.tsv"
  )
)

write_tsv_base(
  DEG.Mesothelial,
  file.path(
    tabledir_Figure1I,
    "Figure1I_DEG_Mesothelial.tsv"
  )
)

# ============================================================
# 6. GO Biological Process de MSigDB mouse
# ============================================================

msig_GO_BP_Figure1I <- msigdbr(
  db_species = "MM",
  species = "Mus musculus",
  collection = "M5",
  subcollection = "GO:BP"
)

msig_GO_BP_Figure1I <- msig_GO_BP_Figure1I %>%
  filter(
    !is.na(gs_name),
    gs_name != "",
    !is.na(gene_symbol),
    gene_symbol != ""
  )

pathways_GO_BP_Figure1I <- split(
  msig_GO_BP_Figure1I$gene_symbol,
  msig_GO_BP_Figure1I$gs_name
)

pathways_GO_BP_Figure1I <- lapply(
  pathways_GO_BP_Figure1I,
  unique
)

write_tsv_base(
  msig_GO_BP_Figure1I,
  file.path(
    tabledir_Figure1I,
    "Figure1I_MSigDB_mouse_GO_BP_gene_sets.tsv"
  )
)

# ============================================================
# 7. Pathways seleccionados
# ============================================================

selected_pathways_Figure1I <- c(
  "EXTRACELLULAR MATRIX DISASSEMBLY",
  "COLLAGEN METABOLIC PROCESS",
  "POSITIVE REGULATION OF EPITHELIAL CELL DIFFERENTIATION",
  "VASODILATION",
  "REGULATION OF TUBE SIZE",
  "ARTERY MORPHOGENESIS",
  "LUNG MORPHOGENESIS",
  "REGULATION OF GLUCOSE METABOLIC PROCESS",
  "P38MAPK CASCADE",
  "AMINO ACID CATABOLIC PROCESS",
  "HYDROGEN PEROXIDE METABOLIC PROCESS"
)

selected_pathway_table_Figure1I <- data.frame(
  Description = selected_pathways_Figure1I
) %>%
  mutate(
    Description_key = normalise_pathway_key(
      Description
    )
  )

# ============================================================
# 8. Construir rankings
# Solo los tres estados de fibroblastos
# ============================================================

make_rank_Figure1I <- function(
    markers,
    cluster_name
) {

  rank_table <- markers %>%
    filter(
      as.character(cluster) == cluster_name,
      !is.na(gene),
      gene != "",
      !is.na(avg_log2FC),
      !is.na(p_val)
    ) %>%
    mutate(
      p_val_rank = pmax(
        p_val,
        1e-300
      ),

      metric_base =
        sign(avg_log2FC) *
        (
          -log10(p_val_rank) +
          abs(avg_log2FC) * 1e-6
        )
    ) %>%
    arrange(
      desc(metric_base),
      desc(avg_log2FC),
      gene
    ) %>%
    distinct(
      gene,
      .keep_all = TRUE
    ) %>%
    mutate(
      # Desempate determinista
      tie_break =
        (
          (n() + 1) / 2 -
          row_number()
        ) * 1e-12,

      metric =
        metric_base +
        tie_break
    ) %>%
    arrange(
      desc(metric)
    )

  rank_vector <- rank_table$metric

  names(rank_vector) <- rank_table$gene

  rank_vector <- sort(
    rank_vector,
    decreasing = TRUE
  )

  list(
    table = rank_table,
    stats = rank_vector
  )
}

rank_results_Figure1I <- lapply(
  fgsea_cluster_order_Figure1I,
  function(cluster_name) {

    make_rank_Figure1I(
      markers = markers_Figure1I,
      cluster_name = cluster_name
    )
  }
)

names(rank_results_Figure1I) <- fgsea_cluster_order_Figure1I

# ============================================================
# 9. Guardar rankings
# ============================================================

rank_tables_all_Figure1I <- bind_rows(
  lapply(
    fgsea_cluster_order_Figure1I,
    function(cluster_name) {

      rank_results_Figure1I[[cluster_name]]$table %>%
        mutate(
          comparison_cluster = cluster_name,
          .before = 1
        )
    }
  )
)

write_tsv_base(
  rank_tables_all_Figure1I,
  file.path(
    tabledir_Figure1I,
    "Figure1I_fgsea_rankings_all_clusters.tsv"
  )
)

for (cluster_name in fgsea_cluster_order_Figure1I) {

  cluster_file <- cluster_files_Figure1I[[cluster_name]]

  write_tsv_base(
    rank_results_Figure1I[[cluster_name]]$table,
    file.path(
      tabledir_Figure1I,
      paste0(
        "Figure1I_fgsea_ranking_",
        cluster_file,
        ".tsv"
      )
    )
  )
}

# ============================================================
# 10. Diagnóstico de tamaño y solapamiento de pathways
# ============================================================

pathway_lookup_Figure1I <- data.frame(
  pathway = names(pathways_GO_BP_Figure1I)
) %>%
  mutate(
    Description = normalise_pathway_name(pathway),
    Description_key = normalise_pathway_key(pathway)
  ) %>%
  filter(
    Description_key %in%
      selected_pathway_table_Figure1I$Description_key
  ) %>%
  distinct(
    Description_key,
    .keep_all = TRUE
  )

minSize_Figure1I <- 15
maxSize_Figure1I <- 500

pathway_diagnostics_Figure1I <- bind_rows(
  lapply(
    seq_len(nrow(pathway_lookup_Figure1I)),
    function(i) {

      pathway_name <- pathway_lookup_Figure1I$pathway[[i]]
      pathway_key <- pathway_lookup_Figure1I$Description_key[[i]]

      gene_set <- pathways_GO_BP_Figure1I[[pathway_name]]

      bind_rows(
        lapply(
          fgsea_cluster_order_Figure1I,
          function(cluster_name) {

            ranking_genes <- names(
              rank_results_Figure1I[[cluster_name]]$stats
            )

            overlap_genes <- intersect(
              gene_set,
              ranking_genes
            )

            data.frame(
              Description_key = pathway_key,
              pathway = pathway_name,
              cluster = cluster_name,
              original_gene_set_size = length(gene_set),
              overlap_with_ranking = length(overlap_genes),
              eligible_for_fgsea =
                length(overlap_genes) >= minSize_Figure1I &&
                length(overlap_genes) <= maxSize_Figure1I,
              overlap_genes = paste(
                overlap_genes,
                collapse = ";"
              )
            )
          }
        )
      )
    }
  )
)

pathway_diagnostics_Figure1I <- selected_pathway_table_Figure1I %>%
  left_join(
    pathway_diagnostics_Figure1I,
    by = "Description_key"
  ) %>%
  arrange(
    match(
      Description,
      selected_pathways_Figure1I
    ),
    match(
      cluster,
      fgsea_cluster_order_Figure1I
    )
  )

write_tsv_base(
  pathway_diagnostics_Figure1I,
  file.path(
    tabledir_Figure1I,
    "Figure1I_selected_pathways_diagnostics.tsv"
  )
)

# ============================================================
# 11. Función fGSEA
#
# minSize = 15 permite incluir pathways más pequeños,
# como EXTRACELLULAR MATRIX DISASSEMBLY.
#
# fgseaSimple reproduce la lógica del análisis original
# con 10.000 permutaciones.
# ============================================================

run_fgsea_Figure1I <- function(
    rank_vector,
    cluster_name
) {

  result <- fgsea::fgseaSimple(
    pathways = pathways_GO_BP_Figure1I,
    stats = rank_vector,
    nperm = 10000,
    minSize = minSize_Figure1I,
    maxSize = maxSize_Figure1I,
    scoreType = "std",
    nproc = 1
  )

  result <- as.data.frame(result)

  result$cluster <- cluster_name

  result %>%
    relocate(cluster) %>%
    arrange(
      desc(NES)
    )
}

# ============================================================
# 12. Ejecutar fGSEA
# No se incluyen mesoteliales
# ============================================================

fgsea_results_Figure1I <- lapply(
  fgsea_cluster_order_Figure1I,
  function(cluster_name) {

    run_fgsea_Figure1I(
      rank_vector = rank_results_Figure1I[[cluster_name]]$stats,
      cluster_name = cluster_name
    )
  }
)

names(fgsea_results_Figure1I) <- fgsea_cluster_order_Figure1I

Universal_fgsea <- fgsea_results_Figure1I[["Fibroblasts Universal"]]
Activated1_fgsea <- fgsea_results_Figure1I[["Fibroblasts Activated 1"]]
Activated2_fgsea <- fgsea_results_Figure1I[["Fibroblasts Activated 2"]]

fgsea_all_Figure1I <- bind_rows(
  fgsea_results_Figure1I
)

# ============================================================
# 13. Guardar todos los resultados fGSEA
# ============================================================

write_tsv_base(
  format_fgsea_table(
    fgsea_all_Figure1I
  ),
  file.path(
    tabledir_Figure1I,
    "Figure1I_fgsea_all_results.tsv"
  )
)

write_tsv_base(
  format_fgsea_table(
    Universal_fgsea
  ),
  file.path(
    tabledir_Figure1I,
    "Figure1I_fgsea_all_Universal.tsv"
  )
)

write_tsv_base(
  format_fgsea_table(
    Activated1_fgsea
  ),
  file.path(
    tabledir_Figure1I,
    "Figure1I_fgsea_all_Activated1.tsv"
  )
)

write_tsv_base(
  format_fgsea_table(
    Activated2_fgsea
  ),
  file.path(
    tabledir_Figure1I,
    "Figure1I_fgsea_all_Activated2.tsv"
  )
)

# ============================================================
# 14. Limpiar nombres de resultados
# ============================================================

df_all_Figure1I <- fgsea_all_Figure1I %>%
  mutate(
    Description = normalise_pathway_name(
      pathway
    ),

    Description_key = normalise_pathway_key(
      pathway
    ),

    condition = unname(
      cluster_labels_Figure1I[
        as.character(cluster)
      ]
    )
  )

write_tsv_base(
  format_fgsea_table(
    df_all_Figure1I
  ),
  file.path(
    tabledir_Figure1I,
    "Figure1I_fgsea_all_results_clean_names.tsv"
  )
)

# ============================================================
# 15. Recuperar pathways seleccionados
# ============================================================

condition_levels_Figure1I <- unname(
  cluster_labels_Figure1I[
    fgsea_cluster_order_Figure1I
  ]
)

df_selected_found_Figure1I <- df_all_Figure1I %>%
  filter(
    Description_key %in%
      selected_pathway_table_Figure1I$Description_key
  ) %>%
  select(
    cluster,
    condition,
    pathway,
    Description_key,
    pval,
    padj,
    ES,
    NES,
    size,
    leadingEdge
  ) %>%
  distinct(
    condition,
    Description_key,
    .keep_all = TRUE
  )

# Todas las combinaciones pathway × fibroblasto
df_plot_Figure1I <- tidyr::expand_grid(
  Description = selected_pathways_Figure1I,
  condition = condition_levels_Figure1I
) %>%
  mutate(
    Description_key = normalise_pathway_key(
      Description
    )
  ) %>%
  left_join(
    df_selected_found_Figure1I,
    by = c(
      "Description_key",
      "condition"
    )
  ) %>%
  mutate(
    condition = factor(
      condition,
      levels = condition_levels_Figure1I
    ),

    Description = factor(
      Description,
      levels = rev(
        selected_pathways_Figure1I
      )
    ),

    padj_plot = case_when(
      is.na(padj) ~ NA_real_,
      padj <= 0 ~ 1e-300,
      TRUE ~ padj
    ),

    minus_log10_padj = case_when(
      is.na(padj_plot) ~ NA_real_,
      TRUE ~ -log10(padj_plot)
    )
  )

# ============================================================
# 16. Comprobar qué pathways aparecen en fGSEA
# ============================================================

pathway_status_Figure1I <- tidyr::expand_grid(
  Description = selected_pathways_Figure1I,
  condition = condition_levels_Figure1I
) %>%
  mutate(
    Description_key = normalise_pathway_key(
      Description
    )
  ) %>%
  left_join(
    df_selected_found_Figure1I %>%
      select(
        condition,
        Description_key,
        pathway,
        pval,
        padj,
        NES,
        size
      ),
    by = c(
      "condition",
      "Description_key"
    )
  ) %>%
  mutate(
    found_in_fgsea = !is.na(pathway),
    has_NES = !is.na(NES)
  )

write_tsv_base(
  pathway_status_Figure1I,
  file.path(
    tabledir_Figure1I,
    "Figure1I_selected_pathways_status.tsv"
  )
)

# ============================================================
# 17. Guardar tabla fuente del plot
# ============================================================

write_tsv_base(
  format_fgsea_table(
    df_plot_Figure1I
  ),
  file.path(
    tabledir_Figure1I,
    "Figure1I_fgsea_selected_pathways_plot.tsv"
  )
)

write_tsv_base(
  format_fgsea_table(
    df_selected_found_Figure1I
  ),
  file.path(
    tabledir_Figure1I,
    "Figure1I_fgsea_selected_pathways_found.tsv"
  )
)

# ============================================================
# 18. Escala simétrica de NES
# ============================================================

NES_values_Figure1I <- df_plot_Figure1I$NES[
  is.finite(df_plot_Figure1I$NES)
]

if (length(NES_values_Figure1I) == 0) {

  NES_limit_Figure1I <- 1

} else {

  NES_limit_Figure1I <- max(
    abs(NES_values_Figure1I)
  )
}

if (
  !is.finite(NES_limit_Figure1I) ||
  NES_limit_Figure1I == 0
) {

  NES_limit_Figure1I <- 1
}

# ============================================================
# 19. Bubble plot
# ============================================================

Figure1I <- ggplot(
  df_plot_Figure1I,
  aes(
    x = condition,
    y = Description
  )
) +
  geom_point(
    aes(
      size = minus_log10_padj,
      color = NES
    ),
    na.rm = TRUE
  ) +
  scale_size_continuous(
    name = "-log10(padj)",
    range = c(2, 8)
  ) +
  scale_color_gradientn(
    name = "NES",
    colors = c(
      "blue",
      "lightblue",
      "white",
      "pink",
      "red"
    ),
    values = scales::rescale(
      c(
        -NES_limit_Figure1I,
        -NES_limit_Figure1I / 2,
        0,
        NES_limit_Figure1I / 2,
        NES_limit_Figure1I
      )
    ),
    limits = c(
      -NES_limit_Figure1I,
      NES_limit_Figure1I
    ),
    oob = scales::squish
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_blank(),

    axis.text.y = element_text(
      size = 9,
      colour = "black"
    ),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9,
      colour = "black"
    ),

    axis.title = element_blank(),

    panel.grid.major = element_line(
      colour = "grey90",
      linewidth = 0.4
    ),

    panel.grid.minor = element_blank(),

    legend.title = element_text(
      size = 9
    ),

    legend.text = element_text(
      size = 8
    )
  )

Figure1I

# ============================================================
# 20. Guardar Figura 1I
# ============================================================

ggsave(
  filename = file.path(
    figdir,
    "Figure1I_fgsea_bubbleplot.pdf"
  ),
  plot = Figure1I,
  width = 10,
  height = 6,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)