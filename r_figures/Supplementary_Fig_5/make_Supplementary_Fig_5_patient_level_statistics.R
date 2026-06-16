#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggtext)
  library(svglite)
  library(ragg)
  library(scales)
})

FIG_TEXT_PT <- 10.5
FIG_PANEL_PT <- 12.5
FIG_HEATMAP_TEXT_PT <- 8.0
FIG_GEOM_TEXT_SIZE <- FIG_TEXT_PT / 2.845276
eta_italic_sym <- "\U0001D702"
eta_caption_sym <- "\u03b7"

phase_short_order <- c("Pre", "Early", "Mid", "Late", "Post")
phase_test_order <- c("Early", "Mid", "Late", "Post")
region_order <- c("SOZ", "PZ", "NIZ")
macro_order <- c("SOZ_only", "PZ_only", "NIZ_only", "SOZ_PZ", "SOZ_NIZ", "PZ_NIZ", "SOZ_PZ_NIZ")
role_order <- c("source-like", "balanced", "sink-like")

out_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
source_csv <- file.path(out_dir, "Supplementary_Fig_5_patient_source_data.csv")
stats_csv <- file.path(out_dir, "Supplementary_Fig_5_patient_summary_stats.csv")
checks_csv <- file.path(out_dir, "Supplementary_Fig_5_preplot_consistency_checks.csv")

required_files <- c(source_csv, stats_csv, checks_csv)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

source_data <- read.csv(source_csv, check.names = FALSE)
stats_data <- read.csv(stats_csv, check.names = FALSE)
checks <- read.csv(checks_csv, check.names = FALSE)

if (any(checks$status != "PASS")) {
  stop("Pre-plot consistency checks did not all pass.")
}

# Figure contract:
# Core conclusion: Regional and hyperedge-level eta-burden patterns observed at
# the seizure level persist after repeated seizures are collapsed to patients.
# Archetype: quantitative grid with patient-level trajectory summaries and a
# paired-test matrix.
# Backend: R only for plot generation, export and visual QA.
# Final size: 183 x 150 mm; PNG, SVG, PDF and high-resolution TIFF outputs.

source_data <- source_data %>%
  mutate(
    phase_short = factor(phase_short, levels = phase_short_order),
    panel = factor(panel, levels = c("a", "b", "c"))
  )

stats_data <- stats_data %>%
  mutate(
    phase_short = factor(phase_short, levels = phase_short_order),
    panel = factor(panel, levels = c("a", "b", "c"))
  )

format_group <- function(x) {
  x <- gsub("_", "-", x, fixed = TRUE)
  x
}

region_palette <- c(SOZ = "#D55E00", PZ = "#0072B2", NIZ = "#009E73")
macro_palette <- c(
  SOZ_only = "#202020",
  PZ_only = "#686868",
  NIZ_only = "#9A9A9A",
  SOZ_PZ = "#D6603D",
  SOZ_NIZ = "#4C78A8",
  PZ_NIZ = "#4E9A4B",
  SOZ_PZ_NIZ = "#8056B3"
)
role_palette <- c(`source-like` = "#CC79A7", balanced = "#5A5A5A", `sink-like` = "#56B4E9")

theme_nature <- function(base_size = FIG_TEXT_PT, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#1E1E1E"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#1E1E1E"),
      axis.title = element_text(size = FIG_TEXT_PT),
      axis.title.y = element_text(margin = margin(r = 1.5)),
      axis.text = element_text(size = FIG_TEXT_PT, colour = "#252525"),
      plot.title = element_text(size = FIG_TEXT_PT, face = "bold", hjust = 0),
      legend.title = element_blank(),
      legend.text = element_text(size = FIG_TEXT_PT),
      legend.key.height = unit(3.2, "mm"),
      legend.key.width = unit(5.2, "mm"),
      legend.spacing.y = unit(0.3, "mm"),
      panel.grid.major.y = element_line(linewidth = 0.18, colour = "#E8EDF2"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(3, 3, 3, 2)
    )
}

summary_line_panel <- function(panel_id, group_levels, palette, title, y_label,
                               y_limits = NULL, show_patient = TRUE,
                               legend_position = c(0.02, 0.98),
                               legend_justification = c(0, 1),
                               legend_ncol = 1) {
  summary_df <- stats_data %>%
    filter(panel == panel_id, group %in% group_levels) %>%
    mutate(group = factor(group, levels = group_levels))
  raw_df <- source_data %>%
    filter(panel == panel_id, group %in% group_levels) %>%
    mutate(group = factor(group, levels = group_levels))

  p <- ggplot(summary_df, aes(x = phase_short, y = mean, group = group, colour = group))

  if (show_patient) {
    p <- p +
      geom_line(
        data = raw_df,
        aes(y = value, group = interaction(patient_id, group), colour = group),
        linewidth = 0.25, alpha = 0.12, show.legend = FALSE
      ) +
      geom_point(
        data = raw_df,
        aes(y = value, colour = group),
        size = 0.75, alpha = 0.16, show.legend = FALSE
      )
  }

  p +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.08, linewidth = 0.38, alpha = 0.95) +
    geom_line(linewidth = 0.62) +
    geom_point(shape = 21, size = 1.85, stroke = 0.35, fill = "white") +
    scale_colour_manual(values = palette, labels = format_group) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = y_limits, labels = label_number(accuracy = 0.001), expand = expansion(mult = c(0.04, 0.12))) +
    labs(title = title, x = NULL, y = y_label) +
    theme_nature() +
    theme(
      legend.position = legend_position,
      legend.justification = legend_justification,
      legend.background = element_rect(fill = alpha("white", 0.78), colour = NA),
      legend.margin = margin(0, 1, 0, 1),
      legend.box.margin = margin(0, 0, 0, 0)
    ) +
    guides(colour = guide_legend(ncol = legend_ncol, override.aes = list(linewidth = 0.8, alpha = 1)))
}

p_region <- summary_line_panel(
  panel_id = "a",
  group_levels = region_order,
  palette = region_palette,
  title = "Regional burden per node",
  y_label = paste0("Patient mean ", eta_italic_sym, " burden per node"),
  y_limits = c(0, 0.038),
  show_patient = TRUE,
  legend_position = c(0.98, 0.98),
  legend_justification = c(1, 1),
  legend_ncol = 1
)

p_macro <- summary_line_panel(
  panel_id = "b",
  group_levels = macro_order,
  palette = macro_palette,
  title = "Macro-state burden",
  y_label = paste0("Patient mean total ", eta_italic_sym, " burden"),
  y_limits = c(0, 0.16),
  show_patient = FALSE,
  legend_position = "bottom",
  legend_justification = c(0.5, 0),
  legend_ncol = 2
) +
  theme(
    legend.key.width = unit(4.3, "mm"),
    legend.text = element_text(size = FIG_TEXT_PT - 1.0),
    legend.margin = margin(0, 0, 0, 0)
  )

p_role <- summary_line_panel(
  panel_id = "c",
  group_levels = role_order,
  palette = role_palette,
  title = "Source/sink role burden",
  y_label = paste0("Patient mean total ", eta_italic_sym, " burden"),
  y_limits = c(0, 0.30),
  show_patient = TRUE,
  legend_position = c(0.98, 0.98),
  legend_justification = c(1, 1),
  legend_ncol = 1
)

heat_order <- c(region_order, macro_order, role_order)
heat_labels <- c(
  SOZ = "SOZ", PZ = "PZ", NIZ = "NIZ",
  SOZ_only = "SOZ-only", PZ_only = "PZ-only", NIZ_only = "NIZ-only",
  SOZ_PZ = "SOZ-PZ", SOZ_NIZ = "SOZ-NIZ", PZ_NIZ = "PZ-NIZ",
  SOZ_PZ_NIZ = "SOZ-PZ-NIZ",
  `source-like` = "source-like", balanced = "balanced", `sink-like` = "sink-like"
)

heat_df <- stats_data %>%
  filter(phase_short != "Pre") %>%
  mutate(
    group = factor(group, levels = rev(heat_order)),
    phase_short = factor(as.character(phase_short), levels = phase_test_order),
    family = case_when(
      panel == "a" ~ "Region",
      panel == "b" ~ "Macro-state",
      TRUE ~ "Role"
    ),
    p = p_signrank_delta_vs_pre,
    p_bucket = case_when(
      is.na(p) ~ "not tested",
      p < 0.001 ~ "p < 0.001",
      p < 0.01 ~ "p < 0.01",
      p < 0.05 ~ "p < 0.05",
      TRUE ~ "n.s."
    ),
    direction = case_when(
      is.na(median_delta_vs_pre) ~ "",
      median_delta_vs_pre > 0 ~ "+",
      median_delta_vs_pre < 0 ~ "-",
      TRUE ~ "0"
    ),
    sig = case_when(
      is.na(p) ~ "",
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    ),
    cell_label = ifelse(sig == "", direction, paste0(direction, sig))
  )

p_heat <- ggplot(heat_df, aes(x = phase_short, y = group, fill = p_bucket)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = cell_label), size = FIG_HEATMAP_TEXT_PT / 2.845276, family = "Arial", colour = "#202020") +
  facet_grid(family ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(labels = function(x) unname(heat_labels[x])) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(
    values = c("p < 0.001" = "#B43B2E", "p < 0.01" = "#D97A49", "p < 0.05" = "#F0B46A", "n.s." = "#E8E8E8", "not tested" = "#F6F6F6"),
    breaks = c("p < 0.001", "p < 0.01", "p < 0.05", "n.s.", "not tested"),
    drop = FALSE
  ) +
  labs(title = "Paired tests vs Pre", x = NULL, y = NULL, fill = NULL) +
  theme_nature() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = FIG_TEXT_PT - 0.2),
    axis.text.y = element_text(size = FIG_TEXT_PT - 0.5),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(size = FIG_TEXT_PT - 1.0, face = "bold", angle = 90, margin = margin(r = 1)),
    legend.position = "none",
    legend.text = element_text(size = FIG_TEXT_PT - 0.8),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(4, "mm"),
    plot.margin = margin(3, 2, 3, 2)
  ) +
  guides(fill = "none")

fig <- (p_region | p_macro) / (p_role | p_heat) +
  plot_layout(widths = c(1, 1.18), heights = c(1.0, 1.08)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = FIG_PANEL_PT, face = "bold", family = "Arial"))

write.csv(
  stats_data %>%
    group_by(panel) %>%
    summarise(
      n_patient_min = min(n_patient, na.rm = TRUE),
      n_patient_max = max(n_patient, na.rm = TRUE),
      n_sig_vs_pre = sum(!is.na(p_signrank_delta_vs_pre) & p_signrank_delta_vs_pre < 0.05),
      .groups = "drop"
    ),
  file.path(out_dir, "Supplementary_Fig_5_postplot_export_checks.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "Patient-level statistics confirm that the regional and hyperedge-level eta-burden patterns are not driven by repeated seizures from a small subset of patients.",
    "",
    paste0("After collapsing seizures within each patient (n = 14 patients), regional per-node ", eta_caption_sym, " burden decreased from the pre-ictal baseline across SOZ, PZ and NIZ during ictal and post-ictal stages, while macro-state and source/sink role summaries preserved the distributed cross-regional and balanced-role pattern observed in the seizure-level analysis."),
    "",
    "Boundary: patient-level values are patient means after within-patient seizure averaging; panels b and c summarize total absolute-proxy burden rather than per-node normalized burden."
  ),
  con = file.path(out_dir, "Supplementary_Fig_5_results_description_draft.txt"),
  useBytes = TRUE
)

writeLines(
  c(
    paste0("Supplementary Fig. 5 | Patient-level statistics for ", eta_caption_sym, " burden."),
    paste0("a, Patient-level regional per-node ", eta_caption_sym, " burden across seizure stages after averaging repeated seizures within each patient. Thin lines and faint points denote individual patient means; thick lines and points denote the cohort mean across patients, with error bars showing s.e.m. b, Patient-level total ", eta_caption_sym, " burden grouped by SOZ, PZ and NIZ macro-state combinations. c, Patient-level total ", eta_caption_sym, " burden grouped by source-like, balanced and sink-like hyperedge roles. d, Paired patient-level Wilcoxon signed-rank tests comparing each non-pre-ictal stage against the pre-ictal stage. Cell signs denote the median direction of change relative to pre-ictal baseline; *, ** and *** denote P < 0.05, P < 0.01 and P < 0.001, respectively. Source data are patient-level means from n = 14 patients, derived from n = 24 seizures.")
  ),
  con = file.path(out_dir, "Supplementary_Fig_5_legend_draft.txt"),
  useBytes = TRUE
)

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 150, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot)
  dev.off()

  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()

  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()

  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = dpi)
  print(plot)
  dev.off()
}

save_pub_r(fig, file.path(out_dir, "Supplementary_Fig_5"))

message("Supplementary Fig. 5 written to: ", out_dir)
