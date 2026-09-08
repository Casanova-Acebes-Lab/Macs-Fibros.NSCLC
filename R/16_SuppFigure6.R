# ============================================================
# 16_SuppFigure6.R
# Supplementary Figure 6
#
# Panel A:
#   Kaplan-Meier plots generated externally with KMplot.
#   Not generated in this repository.
#
# Panel B:
#   EGFR/KRAS mutation-context heatmap
#
# Panel C:
#   TRMs score vs fibroblast signatures in early-stage tumors
#   (stages 1-2)
#
# Human score terminology:
#   Em_Ams -> TRMs
#   MoMacs -> MDMs
#
# Input:
#   data/Figure6/firmas_fibros_en_frescos_16abr26_v3.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(42)

# ============================================================
# Configuration
# ============================================================

if (!file.exists("config.tsv")) {
  stop(
    "No se encuentra config.tsv. ",
    "Ejecuta el script desde la raíz del repositorio."
  )
}

config <- read.delim(
  "config.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

config_values <- setNames(
  config$value,
  config$key
)

rsdir <- unname(
  config_values[["rsdir"]]
)

outdir <- unname(
  config_values[["outdir"]]
)

# ============================================================
# Output directories
# ============================================================

plot_dir <- file.path(
  outdir,
  "SuppFigure6"
)

table_dir <- file.path(
  rsdir,
  "tables",
  "SuppFigure6"
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

write_tsv <- function(data, filename) {

  write.table(
    data,
    file = filename,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

# ============================================================
# Human score input
#
# NOTE:
# Despite the .csv extension, this file is semicolon-separated.
# ============================================================

human_scores_file <- file.path(
  "data",
  "Figure6",
  "firmas_fibros_en_frescos_16abr26_v3.csv"
)

if (!file.exists(human_scores_file)) {
  stop(
    "Human score file not found: ",
    human_scores_file
  )
}

df_all <- read.delim(
  human_scores_file,
  sep = ";",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(df_all) <- trimws(
  colnames(df_all)
)

cat(
  "\nHuman score input:\n",
  human_scores_file,
  "\nSamples: ",
  nrow(df_all),
  "\n",
  sep = ""
)

# ============================================================
# Validate required columns
# ============================================================

required_cols <- c(
  "Sample",
  "Fib_Univ",
  "MoMacs",
  "fib_act1",
  "fib_act2",
  "Em_Ams",
  "subtype",
  "stage",
  "EGFR",
  "KRAS"
)

missing_cols <- setdiff(
  required_cols,
  colnames(df_all)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}

cat(
  "\nRequired columns found.\n"
)

# ============================================================
# Helpers
# ============================================================

zscore_manual <- function(x) {

  x <- as.numeric(x)

  s <- stats::sd(
    x,
    na.rm = TRUE
  )

  if (
    is.na(s) ||
      s == 0
  ) {
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }

  (
    x -
      mean(
        x,
        na.rm = TRUE
      )
  ) / s
}

clean_mutation_status <- function(x) {

  x <- trimws(
    as.character(x)
  )

  dplyr::case_when(

    x %in% c(
      "WT",
      "wt",
      "Wt"
    ) ~ "WT",

    x %in% c(
      "MUT",
      "mut",
      "Mut"
    ) ~ "MUT",

    TRUE ~ NA_character_
  )
}

# ============================================================
# Prepare human scores
# ============================================================

df_human <- df_all

score_vars <- c(
  "Fib_Univ",
  "MoMacs",
  "fib_act1",
  "fib_act2",
  "Em_Ams"
)

df_human[
  score_vars
] <- lapply(
  df_human[
    score_vars
  ],
  function(x) {
    as.numeric(
      as.character(x)
    )
  }
)

df_human$stage_num <- as.numeric(
  as.character(
    df_human$stage
  )
)

# ============================================================
# Cohort-wide z-scores
# ============================================================

df_human$Fib_Univ_z <- zscore_manual(
  df_human$Fib_Univ
)

df_human$MoMacs_z <- zscore_manual(
  df_human$MoMacs
)

df_human$fib_act1_z <- zscore_manual(
  df_human$fib_act1
)

df_human$fib_act2_z <- zscore_manual(
  df_human$fib_act2
)

df_human$Em_Ams_z <- zscore_manual(
  df_human$Em_Ams
)

# ============================================================
# Scores used in Supplementary Figure 6B
# ============================================================

df_human$Activated_CAF_score <- rowMeans(
  df_human[
    ,
    c(
      "fib_act1_z",
      "fib_act2_z"
    )
  ],
  na.rm = TRUE
)

df_human$Universal_Activated_balance <-
  df_human$Fib_Univ_z -
  df_human$Activated_CAF_score

df_human$TRM_MDM_balance <-
  df_human$Em_Ams_z -
  df_human$MoMacs_z

# ============================================================
# ============================================================
# PANEL S6B
# EGFR/KRAS mutation-context heatmap
# ============================================================
# ============================================================

df_s6b <- df_human

# ============================================================
# Mutation status
#
# WT and MUT are informative.
# UK / NV / ND and other values become NA.
# ============================================================

df_s6b$EGFR_status <- clean_mutation_status(
  df_s6b$EGFR
)

df_s6b$KRAS_status <- clean_mutation_status(
  df_s6b$KRAS
)

# ============================================================
# Combined mutation context
# ============================================================

df_s6b$EGFR_KRAS_context <- dplyr::case_when(

  df_s6b$EGFR_status == "WT" &
    df_s6b$KRAS_status == "WT" ~
    "EGFR/KRAS WT",

  df_s6b$EGFR_status == "WT" &
    df_s6b$KRAS_status == "MUT" ~
    "KRAS MUT",

  df_s6b$EGFR_status == "MUT" &
    df_s6b$KRAS_status == "WT" ~
    "EGFR MUT",

  df_s6b$EGFR_status == "MUT" &
    df_s6b$KRAS_status == "MUT" ~
    "EGFR/KRAS MUT",

  df_s6b$KRAS_status == "MUT" &
    is.na(
      df_s6b$EGFR_status
    ) ~
    "KRAS MUT",

  df_s6b$EGFR_status == "MUT" &
    is.na(
      df_s6b$KRAS_status
    ) ~
    "EGFR MUT",

  TRUE ~ NA_character_
)

cat(
  "\nEGFR/KRAS contexts in input:\n"
)

print(
  table(
    df_s6b$EGFR_KRAS_context,
    useNA = "ifany"
  )
)

# ============================================================
# Groups displayed in S6B
# ============================================================

s6b_context_order <- c(
  "EGFR/KRAS WT",
  "KRAS MUT",
  "EGFR MUT"
)

df_s6b <- df_s6b %>%
  filter(
    EGFR_KRAS_context %in%
      s6b_context_order
  ) %>%
  mutate(
    EGFR_KRAS_context =
      factor(
        EGFR_KRAS_context,
        levels =
          s6b_context_order
      )
  ) %>%
  arrange(
    EGFR_KRAS_context,
    stage_num,
    Sample
  )

if (nrow(df_s6b) == 0) {
  stop(
    "No samples remain for Supplementary Figure 6B."
  )
}

cat(
  "\nSamples displayed in S6B:\n"
)

print(
  table(
    df_s6b$EGFR_KRAS_context
  )
)

# ============================================================
# Sample labels
# ============================================================

df_s6b$Sample_stage_label <- ifelse(

  is.na(
    df_s6b$stage_num
  ),

  as.character(
    df_s6b$Sample
  ),

  paste0(
    df_s6b$Sample,
    " | S",
    df_s6b$stage_num
  )
)

df_s6b$Sample_stage_label <- factor(
  df_s6b$Sample_stage_label,
  levels = unique(
    df_s6b$Sample_stage_label
  )
)

# ============================================================
# S6B source data
# ============================================================

s6b_source <- df_s6b %>%
  select(
    Sample,
    Sample_stage_label,
    subtype,
    stage_num,
    EGFR,
    KRAS,
    EGFR_status,
    KRAS_status,
    EGFR_KRAS_context,
    Universal_Activated_balance,
    Activated_CAF_score,
    Fib_Univ_z,
    MoMacs_z,
    Em_Ams_z,
    TRM_MDM_balance
  ) %>%
  pivot_longer(
    cols = c(
      Universal_Activated_balance,
      Activated_CAF_score,
      Fib_Univ_z,
      MoMacs_z,
      Em_Ams_z,
      TRM_MDM_balance
    ),
    names_to =
      "feature_internal",
    values_to =
      "value"
  ) %>%
  mutate(

    feature =
      dplyr::recode(

        feature_internal,

        "Universal_Activated_balance" =
          "Universal / Activated CAF ratio",

        "Activated_CAF_score" =
          "Activated CAF score",

        "Fib_Univ_z" =
          "Fibroblasts Universal",

        "MoMacs_z" =
          "MDMs",

        "Em_Ams_z" =
          "TRMs",

        "TRM_MDM_balance" =
          "TRMs / MDMs ratio"
      )
  )

s6b_feature_order <- c(
  "Universal / Activated CAF ratio",
  "Activated CAF score",
  "Fibroblasts Universal",
  "MDMs",
  "TRMs",
  "TRMs / MDMs ratio"
)

s6b_source$feature <- factor(
  s6b_source$feature,
  levels =
    s6b_feature_order
)

write_tsv(
  s6b_source,
  file.path(
    table_dir,
    "SuppFigure6B_source_data.tsv"
  )
)

# ============================================================
# Facet labels
# ============================================================

s6b_facet_labels <- c(

  "EGFR/KRAS WT" =
    "EGFR*'/'*KRAS^WT",

  "KRAS MUT" =
    "KRAS^MUT",

  "EGFR MUT" =
    "EGFR^MUT"
)

# ============================================================
# S6B heatmap
# ============================================================

p_s6b <- ggplot(
  s6b_source,
  aes(
    x = Sample_stage_label,
    y = feature,
    fill = value
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    limits = c(
      -2,
      2
    ),
    breaks = c(
      -2,
      -1,
      0,
      1,
      2
    ),
    oob = scales::squish,
    name = "Score"
  ) +
  facet_grid(
    . ~ EGFR_KRAS_context,
    scales = "free_x",
    space = "free_x",
    drop = TRUE,
    labeller = labeller(
      EGFR_KRAS_context =
        as_labeller(
          s6b_facet_labels,
          label_parsed
        )
    )
  ) +
  scale_y_discrete(
    limits =
      rev(
        s6b_feature_order
      )
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(

    panel.grid =
      element_blank(),

    strip.background =
      element_rect(
        fill = "white",
        color = "grey35",
        linewidth = 0.5
      ),

    strip.text =
      element_text(
        face = "bold",
        size = 10
      ),

    axis.text.x =
      element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 6.5
      ),

    axis.text.y =
      element_text(
        size = 8
      ),

    axis.ticks =
      element_blank(),

    legend.title =
      element_text(
        size = 8
      ),

    legend.text =
      element_text(
        size = 7
      ),

    legend.position =
      "right"
  )

ggsave(
  filename = file.path(
    plot_dir,
    "SuppFigure6B_EGFR_KRAS_heatmap.pdf"
  ),
  plot = p_s6b,
  width = 12,
  height = 4.4,
  units = "in"
)

cat(
  "\nSupplementary Figure 6B completed.\n"
)

# ============================================================
# ============================================================
# PANEL S6C
# TRMs vs fibroblast signatures
#
# Early-stage tumors only:
# stages 1 and 2
# ============================================================
# ============================================================

df_s6c <- df_human %>%
  filter(
    stage_num %in%
      c(
        1,
        2
      )
  )

if (nrow(df_s6c) == 0) {
  stop(
    "No stage 1-2 samples found for Supplementary Figure 6C."
  )
}

df_s6c$stage_factor <- factor(
  df_s6c$stage_num,
  levels = c(
    1,
    2
  )
)

cat(
  "\nSamples used in S6C:\n"
)

print(
  table(
    df_s6c$stage_factor
  )
)

# ============================================================
# Spearman correlation helper
# ============================================================

run_spearman <- function(
    df,
    yvar
) {

  tmp <- data.frame(

    x =
      df$Em_Ams,

    y =
      df[[yvar]]
  )

  tmp <- tmp[
    complete.cases(
      tmp
    ),
    ,
    drop = FALSE
  ]

  if (nrow(tmp) < 4) {

    return(
      data.frame(

        variable =
          yvar,

        rho =
          NA_real_,

        p_value =
          NA_real_,

        n =
          nrow(tmp),

        stringsAsFactors =
          FALSE
      )
    )
  }

  ct <- suppressWarnings(
    stats::cor.test(
      tmp$x,
      tmp$y,
      method = "spearman",
      exact = FALSE
    )
  )

  data.frame(

    variable =
      yvar,

    rho =
      as.numeric(
        ct$estimate
      ),

    p_value =
      ct$p.value,

    n =
      nrow(tmp),

    stringsAsFactors =
      FALSE
  )
}

# ============================================================
# Three comparisons displayed in S6C
# ============================================================

s6c_variables <- c(
  "Fib_Univ_z",
  "fib_act1_z",
  "fib_act2_z"
)

s6c_stats <- bind_rows(
  lapply(
    s6c_variables,
    function(v) {

      run_spearman(
        df_s6c,
        v
      )
    }
  )
)

# ============================================================
# FDR correction
#
# Only the three comparisons shown in the figure are included.
# ============================================================

s6c_stats$p_adj_fdr <- p.adjust(
  s6c_stats$p_value,
  method = "fdr"
)

s6c_stats$comparison <- dplyr::recode(

  s6c_stats$variable,

  "Fib_Univ_z" =
    "TRMs vs Fibroblasts Universal",

  "fib_act1_z" =
    "TRMs vs Activated 1",

  "fib_act2_z" =
    "TRMs vs Activated 2"
)

s6c_stats$variable <- factor(
  s6c_stats$variable,
  levels =
    s6c_variables
)

s6c_stats <- s6c_stats %>%
  arrange(
    variable
  )

write_tsv(
  s6c_stats,
  file.path(
    table_dir,
    "SuppFigure6C_statistics.tsv"
  )
)

cat(
  "\nS6C statistics:\n"
)

print(
  s6c_stats
)

# ============================================================
# S6C source data
# ============================================================

s6c_source <- df_s6c %>%
  transmute(

    Sample =
      Sample,

    subtype =
      subtype,

    stage =
      stage_num,

    TRMs_score =
      Em_Ams,

    Fibroblasts_Universal =
      Fib_Univ_z,

    Activated_1 =
      fib_act1_z,

    Activated_2 =
      fib_act2_z
  )

write_tsv(
  s6c_source,
  file.path(
    table_dir,
    "SuppFigure6C_source_data.tsv"
  )
)

# ============================================================
# Stage colors
# ============================================================

stage_cols <- c(
  "1" = "#4C9BE8",
  "2" = "#4DAF4A"
)

# ============================================================
# Statistics label
# ============================================================

get_stat_label <- function(yvar) {

  stat_row <- s6c_stats[
    as.character(
      s6c_stats$variable
    ) == yvar,
    ,
    drop = FALSE
  ]

  if (nrow(stat_row) != 1) {
    stop(
      "Statistics not found for: ",
      yvar
    )
  }

  paste0(
    "rho = ",
    format(
      round(
        stat_row$rho,
        3
      ),
      nsmall = 3
    ),
    "; FDR = ",
    signif(
      stat_row$p_adj_fdr,
      3
    ),
    "; n = ",
    stat_row$n
  )
}

# ============================================================
# Scatter helper
# ============================================================

make_scatter <- function(
    yvar,
    y_label
) {

  plot_df <- df_s6c[
    complete.cases(
      df_s6c$Em_Ams,
      df_s6c[[yvar]],
      df_s6c$stage_factor
    ),
    ,
    drop = FALSE
  ]

  ggplot(
    plot_df,
    aes(
      x = Em_Ams,
      y = .data[[yvar]],
      color = stage_factor
    )
  ) +
    geom_point(
      size = 1.6,
      alpha = 0.8
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      color = "grey25",
      fill = "grey80",
      linewidth = 0.8
    ) +
    scale_color_manual(
      values =
        stage_cols,
      name =
        "Stage"
    ) +
    labs(
      title =
        get_stat_label(
          yvar
        ),
      x =
        "TRMs score",
      y =
        y_label
    ) +
    theme_classic(
      base_size = 10.5
    ) +
    theme(

      plot.title =
        element_text(
          size = 10,
          face = "plain",
          hjust = 0
        ),

      legend.position =
        "right"
    )
}

# ============================================================
# Three final S6C scatters
# ============================================================

p_s6c_universal <- make_scatter(
  yvar =
    "Fib_Univ_z",
  y_label =
    "Universal fibroblasts"
)

p_s6c_act1 <- make_scatter(
  yvar =
    "fib_act1_z",
  y_label =
    "Activated fibroblast 1 score"
)

p_s6c_act2 <- make_scatter(
  yvar =
    "fib_act2_z",
  y_label =
    "Activated fibroblast 2 score"
)

# ============================================================
# Combined S6C
# ============================================================

p_s6c <- (
  p_s6c_universal |
    p_s6c_act1 |
    p_s6c_act2
) +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )

ggsave(
  filename = file.path(
    plot_dir,
    "SuppFigure6C_TRMs_Fibroblast_correlations.pdf"
  ),
  plot = p_s6c,
  width = 16,
  height = 5,
  units = "in"
)

cat(
  "\nSupplementary Figure 6C completed.\n"
)

# ============================================================
# Final summary
# ============================================================

cat(
  "\n",
  "============================================\n",
  "SUPPLEMENTARY FIGURE 6 COMPLETE\n",
  "============================================\n",

  "Input:\n",
  human_scores_file,
  "\n\n",

  "S6A:\n",
  "  External Kaplan-Meier plots from KMplot\n",
  "  Not generated by this script\n\n",

  "S6B:\n",
  "  SuppFigure6B_EGFR_KRAS_heatmap.pdf\n",
  "  SuppFigure6B_source_data.tsv\n\n",

  "S6C:\n",
  "  SuppFigure6C_TRMs_Fibroblast_correlations.pdf\n",
  "  SuppFigure6C_source_data.tsv\n",
  "  SuppFigure6C_statistics.tsv\n\n",

  "Plots directory:\n",
  plot_dir,
  "\n\n",

  "Tables directory:\n",
  table_dir,
  "\n",

  "============================================\n",
  sep = ""
)