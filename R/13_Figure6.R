# ============================================================
# Figure 6 - Human validation
# ============================================================

# ============================================================
# Libraries
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ============================================================
# Config
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

config_values <- setNames(config$value, config$key)
required_config_keys <- c("rsdir", "outdir")
missing_config_keys <- setdiff(required_config_keys, names(config_values))

if (length(missing_config_keys) > 0L) {
  stop(
    "Faltan claves en config.tsv: ",
    paste(missing_config_keys, collapse = ", ")
  )
}

rsdir <- unname(config_values[["rsdir"]])
outdir <- unname(config_values[["outdir"]])

figure_dir <- file.path(outdir, "Figure6")
table_dir <- file.path(rsdir, "tables", "Figure6")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Input human signature scores
# ============================================================

# Este CSV contiene scores de firmas humanas previamente derivados a nivel de
# paciente y constituye el input reproducible de Figure 6B-E.
human_input_file <- file.path(
  "data",
  "Figure6",
  "firmas_fibros_en_frescos_16abr26_v3.csv"
)

if (!file.exists(human_input_file)) {
  stop("No se encuentra el input humano: ", human_input_file)
}

df_all <- read.csv(
  human_input_file,
  sep = ";",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

# Guard against BOM-encoded headers from the original source file.
colnames(df_all) <- gsub("^\\xef\\xbb\\xbf", "", colnames(df_all))

required_columns <- c(
  "Sample",
  "Fib_Univ",
  "MoMacs",
  "fib_act1",
  "fib_act2",
  "Em_Ams",
  "stage"
)
missing_columns <- setdiff(required_columns, colnames(df_all))
if (length(missing_columns) > 0L) {
  stop(
    "Faltan columnas requeridas en el input humano: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ============================================================
# Prepare patient-level scores
# ============================================================

zscore_manual <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)

  if (is.na(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }

  (x - mean(x, na.rm = TRUE)) / s
}

# Ensure numeric conversion for the required score columns.
df_final <- df_all
for (colname in c("Fib_Univ", "MoMacs", "fib_act1", "fib_act2", "Em_Ams", "stage")) {
  df_final[[colname]] <- as.numeric(df_final[[colname]])
}

df_final$MDM_score <- zscore_manual(df_final$MoMacs)
df_final$TRM_score <- zscore_manual(df_final$Em_Ams)

df_final$Fibroblasts_Activated1_score <-
  zscore_manual(df_final$fib_act1)

df_final$Fibroblasts_Universal_score <-
  zscore_manual(df_final$Fib_Univ)

df_final$Fibroblasts_Activated2_score <-
  zscore_manual(df_final$fib_act2)

# Fibroblast activation is defined historically as the mean of the
# patient-level z-scores for Activated 1 and Activated 2 fibroblast signatures.
df_final$Fibroblast_activation_score <- rowMeans(
  df_final[, c(
    "Fibroblasts_Activated1_score",
    "Fibroblasts_Activated2_score"
  )],
  na.rm = TRUE
)

# TRM / MDM score in the original figure is NOT a mathematical ratio.
# It is the difference between standardized signature scores:
#
#     z(TRM) - z(MDM)
df_final$TRM_MDM_balance <-
  df_final$TRM_score - df_final$MDM_score

df_final$stage_num <- as.numeric(df_final$stage)
df_final$stage_cat <- ifelse(
  df_final$stage_num %in% c(1, 2),
  "Early",
  ifelse(
    df_final$stage_num %in% c(3, 4),
    "Advanced",
    NA
  )
)

df_final$stage_factor <- factor(
  df_final$stage_cat,
  levels = c("Early", "Advanced")
)

df_final$stage_label <- factor(
  df_final$stage_cat,
  levels = c("Early", "Advanced"),
  labels = c("Early (Stage I-II)", "Advanced (Stage III-IV)")
)

# stage individual follows the historical numeric coding.
df_final$stage_individual <- factor(
  df_final$stage_num,
  levels = c(1, 2, 3, 4)
)

print(table(df_final$stage_factor))

# ============================================================
# Helper Spearman
# ============================================================

spearman_complete <- function(x, y) {
  keep <- complete.cases(x, y)
  x_keep <- x[keep]
  y_keep <- y[keep]

  if (length(x_keep) < 5L) {
    return(data.frame(rho = NA_real_, p_value = NA_real_, n = length(x_keep)))
  }

  cor_test <- stats::cor.test(
    x_keep,
    y_keep,
    method = "spearman",
    exact = FALSE
  )

  data.frame(
    rho = as.numeric(cor_test$estimate),
    p_value = as.numeric(cor_test$p.value),
    n = length(x_keep),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# FIGURE 6B - Correlation matrix
# ============================================================

score_order <- c(
  "MDM_score",
  "TRM_score",
  "Fibroblasts_Activated1_score",
  "Fibroblasts_Universal_score",
  "Fibroblasts_Activated2_score"
)

score_label_map <- c(
  "MDM_score" = "MDMs",
  "TRM_score" = "TRMs",
  "Fibroblasts_Activated1_score" = "Fibroblasts Activated 1",
  "Fibroblasts_Universal_score" = "Fibroblasts Universal",
  "Fibroblasts_Activated2_score" = "Fibroblasts Activated 2"
)

pair_df <- expand.grid(
  var1 = score_order,
  var2 = score_order,
  stringsAsFactors = FALSE
) %>%
  filter(var1 != var2) %>%
  mutate(var1 = as.character(var1), var2 = as.character(var2))

pair_df <- pair_df %>%
  rowwise() %>%
  mutate(
    result = list(
      spearman_complete(
        df_final[[var1]],
        df_final[[var2]]
      )
    )
  ) %>%
  ungroup() %>%
  tidyr::unnest_wider(result) %>%
  mutate(
    p_adj_fdr = p.adjust(p_value, method = "fdr")
  )

# Add the symmetric pairs to make the heatmap matrix complete.
all_pairs <- bind_rows(
  pair_df %>% transmute(var1, var2, rho, p_value, p_adj_fdr, n),
  pair_df %>% transmute(
    var1 = var2,
    var2 = var1,
    rho = rho,
    p_value = p_value,
    p_adj_fdr = p_adj_fdr,
    n = n
  )
)

rho_matrix <- matrix(1, nrow = length(score_order), ncol = length(score_order))
rownames(rho_matrix) <- score_order
colnames(rho_matrix) <- score_order

for (i in seq_len(nrow(pair_df))) {
  row_i <- pair_df[i, ]
  rho_matrix[row_i$var1, row_i$var2] <- row_i$rho
  rho_matrix[row_i$var2, row_i$var1] <- row_i$rho
}

rho_long <- as.data.frame(as.table(rho_matrix), stringsAsFactors = FALSE)
colnames(rho_long) <- c("var1", "var2", "rho")

if (anyNA(rho_long$rho)) {
  rho_long$rho[is.na(rho_long$rho)] <- 0
}

rho_long <- rho_long %>%
  mutate(
    var1 = as.character(var1),
    var2 = as.character(var2),
    label1 = score_label_map[var1],
    label2 = score_label_map[var2]
  )

rho_long <- rho_long %>%
  left_join(
    all_pairs %>%
      rename(rho_pair = rho, p_value = p_value, p_adj_fdr = p_adj_fdr, n = n),
    by = c("var1", "var2")
  )

# Rebuild the pair table with the required source-data columns.
figure6b_source <- pair_df %>%
  transmute(
    var1 = var1,
    var2 = var2,
    rho = rho,
    p_value = p_value,
    p_adj_fdr = p_adj_fdr,
    n = n
  ) %>%
  arrange(var1, var2)

write.table(
  figure6b_source,
  file = file.path(table_dir, "Figure6B_source_data.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

rho_plot <- rho_long %>%
  mutate(
    var1 = factor(var1, levels = score_order),
    var2 = factor(var2, levels = rev(score_order))
  ) %>%
  ggplot(aes(x = var1, y = var2, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3.5) +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman rho",
    na.value = "grey90"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_blank(),
    legend.position = "right"
  ) +
  scale_x_discrete(labels = score_label_map[score_order]) +
  scale_y_discrete(labels = rev(score_label_map[score_order]))

ggsave(
  filename = file.path(figure_dir, "Figure6B_Human_State_Correlations.pdf"),
  plot = rho_plot,
  width = 7,
  height = 6,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

# ============================================================
# Figure 6C - TRM/MDM balance by tumor stage
# ============================================================

figure6c_df <- df_final %>%
  filter(!is.na(stage_cat)) %>%
  transmute(
    Sample = Sample,
    stage = stage,
    stage_cat = stage_cat,
    TRM_score = TRM_score,
    MDM_score = MDM_score,
    TRM_MDM_balance = TRM_MDM_balance
  )

wilcox_result <- stats::wilcox.test(
  TRM_MDM_balance ~ stage_cat,
  data = figure6c_df,
  exact = FALSE
)

figure6c_plot_df <- figure6c_df %>%
  mutate(
    stage_plot = factor(
      stage_cat,
      levels = c("Early", "Advanced"),
      labels = c("Early (Stage I-II)", "Advanced (Stage III-IV)")
    )
  )

figure6c_plot <- ggplot(
  figure6c_plot_df,
  aes(x = stage_plot, y = TRM_MDM_balance, fill = stage_plot)
) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(
    values = c(
      "Early (Stage I-II)" = "#E5E5E5",
      "Advanced (Stage III-IV)" = "#D73027"
    )
  ) +
  labs(x = NULL, y = "TRM / MDM score") +
  theme_classic(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  filename = file.path(figure_dir, "Figure6C_TRM_MDM_by_Stage.pdf"),
  plot = figure6c_plot,
  width = 5.5,
  height = 4.5,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

write.table(
  figure6c_df,
  file = file.path(table_dir, "Figure6C_source_data.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

figure6c_stats <- data.frame(
  method = "Wilcoxon rank-sum test",
  statistic = as.numeric(wilcox_result$statistic),
  p_value = as.numeric(wilcox_result$p.value),
  alternative = as.character(wilcox_result$alternative),
  stringsAsFactors = FALSE
)

write.table(
  figure6c_stats,
  file = file.path(table_dir, "Figure6C_Wilcoxon.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# ============================================================
# Figure 6D - TRM vs fibroblast activation
# ============================================================

spearman_d <- stats::cor.test(
  df_final$TRM_score,
  df_final$Fibroblast_activation_score,
  method = "spearman",
  exact = FALSE
)

p_adj_fdr <- p.adjust(spearman_d$p.value, method = "fdr")
# With a single p-value, the FDR adjustment is numerically identical to the
# original p-value, but the label is kept for historical consistency.

figure6d_label <- sprintf(
  "rho = %.2f; FDR = %.3g; n = %d",
  unname(spearman_d$estimate),
  p_adj_fdr,
  sum(complete.cases(df_final$TRM_score, df_final$Fibroblast_activation_score))
)

figure6d_plot <- ggplot(
  df_final,
  aes(x = TRM_score, y = Fibroblast_activation_score)
) +
  geom_point(color = "black", size = 2) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "red",
    linetype = 2,
    fill = "#F4A3A3",
    alpha = 0.35
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = figure6d_label,
    hjust = 1.02,
    vjust = 1.2,
    size = 4
  ) +
  labs(x = "TRM score", y = "Fibroblast activation score") +
  theme_classic(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  filename = file.path(figure_dir, "Figure6D_TRM_FibroblastActivation_Correlation.pdf"),
  plot = figure6d_plot,
  width = 7,
  height = 6,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

figure6d_source <- df_final %>%
  transmute(
    Sample = Sample,
    TRM_score = TRM_score,
    Fibroblast_activation_score = Fibroblast_activation_score
  )

write.table(
  figure6d_source,
  file = file.path(table_dir, "Figure6D_source_data.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

figure6d_stats <- data.frame(
  method = "spearman",
  rho = as.numeric(spearman_d$estimate),
  p_value = as.numeric(spearman_d$p.value),
  p_adj_fdr = as.numeric(p_adj_fdr),
  n = sum(complete.cases(df_final$TRM_score, df_final$Fibroblast_activation_score)),
  stringsAsFactors = FALSE
)

write.table(
  figure6d_stats,
  file = file.path(table_dir, "Figure6D_Spearman.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# ============================================================
# FIGURE 6E - Dominant fibroblast state by stage
# ============================================================

score_state_order <- c(
  "Fibroblasts Universal",
  "Fibroblasts Activated 1",
  "Fibroblasts Activated 2"
)

score_state_map <- c(
  "Fibroblasts Universal" = "Fibroblasts_Universal_score",
  "Fibroblasts Activated 1" = "Fibroblasts_Activated1_score",
  "Fibroblasts Activated 2" = "Fibroblasts_Activated2_score"
)

Figure6_fibro_colors <- c(
  "Fibroblasts Universal"   = "#FBCDC7",
  "Fibroblasts Activated 1" = "#F94144",
  "Fibroblasts Activated 2" = "#6D1A1A"
)

df_final$fibroblast_dominant_state <- NA_character_

for (i in seq_len(nrow(df_final))) {
  current_scores <- as.numeric(
    df_final[i, c(
      "Fibroblasts_Universal_score",
      "Fibroblasts_Activated1_score",
      "Fibroblasts_Activated2_score"
    )]
  )

  if (all(is.na(current_scores))) {
    df_final$fibroblast_dominant_state[i] <- NA_character_
  } else {
    df_final$fibroblast_dominant_state[i] <- score_state_order[which.max(current_scores)]
  }
}

df_final$fibroblast_dominant_state <- factor(
  df_final$fibroblast_dominant_state,
  levels = score_state_order
)

figure6e_df <- df_final %>%
  filter(!is.na(stage_num), !is.na(fibroblast_dominant_state)) %>%
  count(stage_num, fibroblast_dominant_state, name = "n") %>%
  group_by(stage_num) %>%
  mutate(
    total_n = sum(n),
    pct = 100 * n / total_n
  ) %>%
  ungroup() %>%
  mutate(
    stage = as.integer(stage_num),
    fibroblast_state = as.character(fibroblast_dominant_state)
  )

stage_total_counts <- table(df_final$stage_num)
figure6e_df$stage_plot <- factor(
  figure6e_df$stage_num,
  levels = c(1, 2, 3, 4),
  labels = c(
    paste0("1\n(n=", stage_total_counts[["1"]], ")"),
    paste0("2\n(n=", stage_total_counts[["2"]], ")"),
    paste0("3\n(n=", stage_total_counts[["3"]], ")"),
    paste0("4\n(n=", stage_total_counts[["4"]], ")")
  )
)

figure6e_df$fibroblast_state <- factor(
  figure6e_df$fibroblast_state,
  levels = c(
    "Fibroblasts Universal",
    "Fibroblasts Activated 1",
    "Fibroblasts Activated 2"
  )
)

figure6e_plot <- ggplot(
  figure6e_df,
  aes(x = stage_plot, y = pct, fill = fibroblast_state)
) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(
    values = Figure6_fibro_colors,
    name = "Fibroblast state"
  ) +
  labs(x = "Tumor stage", y = "% of patients") +
  theme_classic(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100))

ggsave(
  filename = file.path(figure_dir, "Figure6E_FibroblastStates_by_Stage.pdf"),
  plot = figure6e_plot,
  width = 6.5,
  height = 4.5,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

write.table(
  figure6e_df %>%
    transmute(
      stage = as.integer(stage_num),
      fibroblast_state = fibroblast_state,
      n = n,
      pct = pct,
      total_n = total_n
    ),
  file = file.path(table_dir, "Figure6E_source_data.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

patient_score_df <- df_final %>%
  transmute(
    Sample = Sample,
    stage = stage,
    stage_cat = stage_cat,
    MDM_score = MDM_score,
    TRM_score = TRM_score,
    Fibroblasts_Activated1_score = Fibroblasts_Activated1_score,
    Fibroblasts_Universal_score = Fibroblasts_Universal_score,
    Fibroblasts_Activated2_score = Fibroblasts_Activated2_score,
    Fibroblast_activation_score = Fibroblast_activation_score,
    TRM_MDM_balance = TRM_MDM_balance,
    fibroblast_dominant_state = as.character(fibroblast_dominant_state)
  )

write.table(
  patient_score_df,
  file = file.path(table_dir, "Figure6_patient_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# ============================================================
# 8. NO generar figura combinada
# ============================================================

# ============================================================
# 9. Validaciones de reproducción
# ============================================================

cat("Total pacientes:", nrow(df_final), "\n")
print(table(df_final$stage_num))
print(table(df_final$stage_cat, useNA = "ifany"))
print(figure6b_source)
print(wilcox_result)
print(spearman_d)
print(p_adj_fdr)
print(figure6e_df %>%
  group_by(stage_num) %>%
  summarise(total_n = unique(total_n), .groups = "drop") %>%
  mutate(pct = 100 * total_n / sum(total_n)))

# ============================================================
# Session information
# ============================================================

sessionInfo()
