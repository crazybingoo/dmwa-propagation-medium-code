suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(svglite)
  library(ragg)
  library(scales)
  library(grid)
})

FIG_TEXT_PT <- 10.5
FIG_PANEL_PT <- 12.5
FIG_GEOM_TEXT_SIZE <- FIG_TEXT_PT / 2.845276

base_candidates <- Sys.glob("example_project/*0514-/nature_fig/Supplementary_Fig_3")
if (length(base_candidates) < 1) {
  stop("Cannot locate Supplementary_Fig_3 directory under example_project/*0514-/nature_fig")
}
base_dir <- normalizePath(base_candidates[1], winslash = "/", mustWork = TRUE)
out_base <- file.path(base_dir, "Supplementary_Fig_3")

phase_levels <- c("early", "mid", "late", "post-ictal")
phase_labels <- c("Early", "Mid", "Late", "Post")
control_levels <- c("DensityMatched", "WeightShuffled", "TopologyShuffled")
control_labels <- c("Density matched", "Weight shuffled", "Topology shuffled")
model_levels <- c("Original", control_levels)
model_labels <- c("Original", control_labels)

fig1_palette <- function() {
  c(
    ink = "#202124",
    neutral_dark = "#4A4F55",
    neutral_mid = "#9AA3AD",
    neutral_light = "#D9DEE7",
    neutral_fill = "#F2F5F8",
    signal = "#2B7A78",
    signal_light = "#D7ECEA",
    accent = "#C4524A",
    accent_light = "#F6D7D2"
  )
}

pal <- fig1_palette()

model_fill_cols <- c(
  "Original" = unname(pal["neutral_fill"]),
  "Density matched" = unname(pal["signal_light"]),
  "Weight shuffled" = unname(pal["signal_light"]),
  "Topology shuffled" = unname(pal["signal_light"])
)

model_line_cols <- c(
  "Original" = unname(pal["ink"]),
  "Density matched" = unname(pal["signal"]),
  "Weight shuffled" = unname(pal["signal"]),
  "Topology shuffled" = unname(pal["signal"])
)

model_violin_fill_cols <- c(
  "Original" = alpha(unname(pal["neutral_fill"]), 0.30),
  "Density matched" = alpha(unname(pal["signal_light"]), 0.24),
  "Weight shuffled" = alpha(unname(pal["signal_light"]), 0.24),
  "Topology shuffled" = alpha(unname(pal["signal_light"]), 0.24)
)

model_box_fill_cols <- c(
  "Original" = alpha("white", 0.82),
  "Density matched" = alpha("white", 0.78),
  "Weight shuffled" = alpha("white", 0.78),
  "Topology shuffled" = alpha("white", 0.78)
)

model_violin_line_cols <- c(
  "Original" = alpha(unname(pal["neutral_dark"]), 0.58),
  "Density matched" = alpha(unname(pal["signal"]), 0.54),
  "Weight shuffled" = alpha(unname(pal["signal"]), 0.54),
  "Topology shuffled" = alpha(unname(pal["signal"]), 0.54)
)

model_point_cols <- c(
  "Original" = unname(pal["neutral_dark"]),
  "Density matched" = unname(pal["signal"]),
  "Weight shuffled" = unname(pal["signal"]),
  "Topology shuffled" = unname(pal["signal"])
)

model_mean_fill_cols <- c(
  "Original" = unname(pal["neutral_fill"]),
  "Density matched" = unname(pal["signal_light"]),
  "Weight shuffled" = unname(pal["signal_light"]),
  "Topology shuffled" = unname(pal["signal_light"])
)

model_mean_line_cols <- c(
  "Original" = unname(pal["ink"]),
  "Density matched" = unname(pal["signal"]),
  "Weight shuffled" = unname(pal["signal"]),
  "Topology shuffled" = unname(pal["signal"])
)

DODGE_WIDTH <- 0.84
VIOLIN_WIDTH <- 0.58
BOX_WIDTH <- 0.13
POINT_SIZE <- 1.45
POINT_ALPHA <- 0.68
POINT_JITTER_WIDTH <- 0.085

theme_nature_fig1 <- function(base_size = FIG_TEXT_PT) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = pal["ink"]),
      axis.ticks = element_line(linewidth = 0.30, colour = pal["ink"]),
      axis.text = element_text(size = FIG_TEXT_PT, colour = pal["ink"]),
      axis.title = element_text(size = FIG_TEXT_PT, colour = pal["ink"]),
      axis.title.y = element_text(margin = margin(r = 1.5)),
      plot.title = element_text(size = FIG_TEXT_PT, face = "bold", hjust = 0, colour = pal["ink"]),
      legend.title = element_text(size = FIG_TEXT_PT, colour = pal["ink"]),
      legend.text = element_text(size = FIG_TEXT_PT, colour = pal["ink"]),
      legend.key.height = unit(3.0, "mm"),
      legend.key.width = unit(5.0, "mm"),
      legend.background = element_blank(),
      panel.grid.major = element_line(linewidth = 0.18, colour = "#ECEFF3"),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = FIG_TEXT_PT, face = "bold", colour = pal["ink"]),
      plot.margin = margin(3.0, 3.0, 3.0, 3.0)
    )
}
theme_set(theme_nature_fig1())

stars_from_p <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", "")
      )
    )
  )
}

delta <- read_csv(
  file.path(base_dir, "Supplementary_Fig_3_delta_eta_source.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    model = factor(model, levels = model_levels),
    model_label = factor(model_label, levels = model_labels),
    phase5 = factor(phase5, levels = phase_levels),
    phase_label = factor(phase_label, levels = phase_labels)
  )

summary_tbl <- read_csv(
  file.path(base_dir, "Supplementary_Fig_3_delta_eta_summary.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    model = factor(model, levels = model_levels),
    model_label = factor(model_label, levels = model_labels),
    phase5 = factor(phase5, levels = phase_levels),
    phase_label = factor(phase_label, levels = phase_labels)
  )

pair_stats <- read_csv(
  file.path(base_dir, "Supplementary_Fig_3_delta_original_vs_controls_stats.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    control_model = factor(control_model, levels = control_levels),
    control_label = factor(control_label, levels = control_labels),
    phase5 = factor(phase5, levels = phase_levels),
    phase_label = factor(phase_label, levels = phase_labels),
    stars = stars_from_p(p_holm_by_control)
  )

plot_data <- bind_rows(lapply(seq_along(control_levels), function(i) {
  ctrl <- control_levels[i]
  ctrl_label <- control_labels[i]
  delta %>%
    filter(model %in% c("Original", ctrl)) %>%
    mutate(
      control_model = ctrl,
      control_label = ctrl_label,
      control_label = factor(control_label, levels = control_labels)
    )
})) %>%
  mutate(
    model_label_chr = as.character(model_label),
    violin_fill = unname(model_violin_fill_cols[model_label_chr]),
    box_fill = unname(model_box_fill_cols[model_label_chr]),
    violin_line_col = unname(model_violin_line_cols[model_label_chr]),
    line_col = unname(model_line_cols[model_label_chr]),
    point_col = unname(model_point_cols[model_label_chr])
  )

mean_data <- bind_rows(lapply(seq_along(control_levels), function(i) {
  ctrl <- control_levels[i]
  ctrl_label <- control_labels[i]
  summary_tbl %>%
    filter(model %in% c("Original", ctrl)) %>%
    mutate(
      control_model = ctrl,
      control_label = ctrl_label,
      control_label = factor(control_label, levels = control_labels)
    )
})) %>%
  mutate(
    model_label_chr = as.character(model_label),
    line_col = unname(model_line_cols[model_label_chr]),
    mean_fill = unname(model_mean_fill_cols[model_label_chr]),
    mean_line = unname(model_mean_line_cols[model_label_chr])
  )

annot_data <- plot_data %>%
  group_by(control_label, phase_label) %>%
  summarise(y = max(delta_eta, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    pair_stats %>%
      select(control_label, phase_label, stars, p_holm_by_control),
    by = c("control_label", "phase_label")
  ) %>%
  mutate(
    y = y + 0.0045,
    label = stars
  ) %>%
  filter(label != "")

y_limits <- range(plot_data$delta_eta, na.rm = TRUE)
y_pad <- diff(y_limits) * 0.12
y_min <- min(-0.004, y_limits[1] - y_pad * 0.2)
y_max <- max(annot_data$y, y_limits[2], na.rm = TRUE) + y_pad * 0.35

set.seed(260615)

make_control_panel <- function(ctrl_label, show_y = FALSE) {
  panel_data <- plot_data %>% filter(control_label == ctrl_label)
  panel_mean <- mean_data %>% filter(control_label == ctrl_label)
  panel_annot <- annot_data %>% filter(control_label == ctrl_label)

  ggplot(panel_data, aes(x = phase_label, y = delta_eta)) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = pal["neutral_dark"]) +
    geom_violin(
      aes(
        group = interaction(phase_label, model_label),
        fill = violin_fill,
        colour = violin_line_col
      ),
      position = position_dodge(width = DODGE_WIDTH),
      width = VIOLIN_WIDTH,
      linewidth = 0.20,
      trim = TRUE,
      scale = "width"
    ) +
    geom_boxplot(
      aes(
        group = interaction(phase_label, model_label),
        fill = box_fill,
        colour = line_col
      ),
      position = position_dodge(width = DODGE_WIDTH),
      width = BOX_WIDTH,
      linewidth = 0.28,
      outlier.shape = NA
    ) +
    geom_point(
      aes(
        group = model_label,
        colour = point_col
      ),
      position = position_jitterdodge(
        jitter.width = POINT_JITTER_WIDTH,
        jitter.height = 0,
        dodge.width = DODGE_WIDTH,
        seed = 260615
      ),
      shape = 16,
      size = POINT_SIZE,
      alpha = POINT_ALPHA,
      stroke = 0
    ) +
    geom_point(
      data = panel_mean,
      aes(
        x = phase_label,
        y = mean_delta_eta,
        group = model_label,
        fill = mean_fill,
        colour = mean_line
      ),
      inherit.aes = FALSE,
      position = position_dodge(width = DODGE_WIDTH),
      shape = 23,
      size = 2.15,
      stroke = 0.40,
    ) +
    geom_text(
      data = panel_annot,
      aes(x = phase_label, y = y, label = label),
      inherit.aes = FALSE,
      size = FIG_GEOM_TEXT_SIZE,
      fontface = "bold",
      colour = pal["ink"],
      vjust = 0
    ) +
    scale_fill_identity() +
    scale_colour_identity() +
    scale_x_discrete(expand = expansion(add = 0.38)) +
    scale_y_continuous(
      labels = label_number(accuracy = 0.005),
      breaks = pretty_breaks(n = 6),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
    labs(
      title = as.character(ctrl_label),
      x = NULL,
      y = if (show_y) expression(Delta * italic(eta) * " relative to pre-ictal") else NULL
    ) +
    guides(fill = "none", colour = "none") +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.margin = margin(0, 0, 1, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
}

p <- (
  make_control_panel("Density matched", show_y = TRUE) |
    make_control_panel("Weight shuffled", show_y = FALSE) |
    make_control_panel("Topology shuffled", show_y = FALSE)
) +
  plot_layout(guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "none",
    plot.tag = element_text(size = FIG_PANEL_PT, face = "bold", colour = pal["ink"])
  )

save_pub_r <- function(plot, filename, width_mm = 180, height_mm = 82, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = dpi, background = "white")
  print(plot)
  dev.off()

  svglite::svglite(paste0(filename, ".svg"), width = w, height = h, bg = "white")
  print(plot)
  dev.off()

  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial", bg = "white")
  print(plot)
  dev.off()

  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi, compression = "lzw", background = "white")
  print(plot)
  dev.off()
}

save_pub_r(p, out_base, width_mm = 183, height_mm = 82, dpi = 600)

write_sidecar_files <- FALSE
if (write_sidecar_files) {
export_checks <- tibble(
  file = paste0("Supplementary_Fig_3.", c("png", "svg", "pdf", "tiff")),
  exists = file.exists(file.path(base_dir, file)),
  bytes = file.info(file.path(base_dir, file))$size
)
write_csv(export_checks, file.path(base_dir, "Supplementary_Fig_3_postplot_export_checks.csv"))

legend_text <- paste0(
  "Supplementary Fig. 3 | Null-control distributions for \u0394\u03b7. ",
  "Distributions of case-level \u0394\u03b7, defined as the model-specific stage-wise mean \u03b7 in each seizure phase minus the pre-ictal mean \u03b7 from the same seizure, are shown for the original DMWA and each null-control model. ",
  "a-c, Original versus density-matched, weight-shuffled and topology-shuffled controls, respectively. ",
  "Grey/black distributions denote the original DMWA; coloured distributions denote the indicated control model in each panel. ",
  "Points denote individual seizures (n = 24); boxes show the median and interquartile range with 1.5\u00d7 IQR whiskers; white diamonds show the mean. ",
  "P values compare original and control \u0394\u03b7 within seizures using two-sided Wilcoxon signed-rank tests with Holm correction within each control family; *P < 0.05, **P < 0.01, ***P < 0.001. ",
  "Random controls used 30 randomizations per window."
)
writeLines(legend_text, file.path(base_dir, "Supplementary_Fig_3_legend_draft.txt"), useBytes = TRUE)

results_text <- paste0(
  "Supplementary Fig. 3 shows the case-level distributions underlying the null-control \u0394\u03b7 analysis. ",
  "The density-matched control retained positive but smaller \u0394\u03b7 than the original DMWA across all post-pre-ictal phases, whereas the weight-shuffled control produced broadly similar positive \u0394\u03b7 distributions. ",
  "Topology shuffling collapsed the phase-related \u0394\u03b7 distribution toward zero, supporting the interpretation that the DMWA topology contributes to the observed propagation-medium state."
)
writeLines(results_text, file.path(base_dir, "Supplementary_Fig_3_results_description_draft.txt"), useBytes = TRUE)

qa_text <- c(
  "Supplementary Fig. 3 QA notes",
  "Core conclusion: null-control distributions expose how density matching, weight shuffling and topology shuffling alter case-level delta eta relative to pre-ictal baseline.",
  "Evidence chain: each facet compares the original DMWA distribution with one control model across early, mid, late and post-ictal phases.",
  "Source data: MATLAB-derived Supplementary_Fig_3_delta_eta_source.csv; delta_eta = eta_mean(phase) - eta_mean(pre-ictal), paired within case and model.",
  "Statistics: two-sided Wilcoxon signed-rank tests comparing original versus control delta_eta within seizures; Holm correction within each control family.",
  "Visual QA checklist: panel titles are bold; delta eta axis uses italic eta; zero line visible; legend does not overlap data; dots, boxes and mean diamonds are visible."
)
writeLines(qa_text, file.path(base_dir, "Supplementary_Fig_3_QA_notes.txt"), useBytes = TRUE)
}

message("Wrote Supplementary Fig. 3 outputs to: ", base_dir)
