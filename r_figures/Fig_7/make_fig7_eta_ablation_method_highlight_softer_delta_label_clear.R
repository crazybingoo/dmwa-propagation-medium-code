#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(scales)
  library(grid)
  library(tibble)
})

out_dir <- file.path("F:/6", "\u521d\u7a3f0514-", "nature_fig", "Fig_7")

needed <- c(
  "Fig7_stage_source_data.csv",
  "Fig7_stage_summary.csv",
  "Fig7_paired_stage_stats.csv",
  "Fig7_variance_source_data.csv",
  "Fig7_roc_auc_summary.csv",
  "Fig7_roc_coordinates.csv",
  "Fig7_preplot_consistency_checks.csv"
)
missing <- needed[!file.exists(file.path(out_dir, needed))]
if (length(missing) > 0) {
  stop("Missing previously audited Fig. 7 data files: ", paste(missing, collapse = ", "))
}

pre_checks <- read.csv(file.path(out_dir, "Fig7_preplot_consistency_checks.csv"), stringsAsFactors = FALSE)
if (any(pre_checks$status != "PASS")) {
  stop("Pre-plot checks are not all PASS; inspect Fig7_preplot_consistency_checks.csv")
}

stage_long <- read.csv(file.path(out_dir, "Fig7_stage_source_data.csv"), stringsAsFactors = FALSE)
stage_summary <- read.csv(file.path(out_dir, "Fig7_stage_summary.csv"), stringsAsFactors = FALSE)
paired_stats <- read.csv(file.path(out_dir, "Fig7_paired_stage_stats.csv"), stringsAsFactors = FALSE)
variance_calc <- read.csv(file.path(out_dir, "Fig7_variance_source_data.csv"), stringsAsFactors = FALSE)
roc_stats <- read.csv(file.path(out_dir, "Fig7_roc_auc_summary.csv"), stringsAsFactors = FALSE)
roc_coords <- read.csv(file.path(out_dir, "Fig7_roc_coordinates.csv"), stringsAsFactors = FALSE)

phase_labels <- c("Pre", "Early", "Mid", "Late", "Post")
model_levels <- c("Original", "Coverage only", "DegreeBias only")

stage_long <- stage_long %>%
  mutate(
    phase_label = factor(phase_label, levels = phase_labels),
    model = factor(model, levels = model_levels)
  )
stage_summary <- stage_summary %>%
  mutate(
    phase_label = factor(phase_label, levels = phase_labels),
    model = factor(model, levels = model_levels)
  )
paired_stats <- paired_stats %>% mutate(model = factor(model, levels = model_levels))
variance_calc <- variance_calc %>%
  mutate(
    phase_label = factor(phase_label, levels = phase_labels),
    model = factor(model, levels = model_levels)
  )
roc_stats <- roc_stats %>% mutate(model = factor(model, levels = model_levels))
roc_coords <- roc_coords %>% mutate(model = factor(model, levels = model_levels))

# Method-highlighting palette:
# Original is the accent but in a softer coral-red so the main method is
# foregrounded without feeling too saturated or too warm.
model_palette <- c(
  "Original" = "#D96661",
  "Coverage only" = "#4778A8",
  "DegreeBias only" = "#7A828C"
)

theme_nature <- function(base_size = 6.6, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(colour = "#202124"),
      axis.line = element_line(linewidth = 0.35, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#222222"),
      axis.title = element_text(size = base_size + 0.1),
      axis.text = element_text(size = base_size - 0.2, colour = "#333333"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.3),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 0.1, face = "bold", colour = "#222222"),
      plot.title = element_text(size = base_size + 0.8, face = "bold", hjust = 0, colour = "#202124"),
      plot.margin = margin(4, 5, 4, 5),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

sig_short <- paired_stats %>%
  filter(comparison == "pre-ictal vs early") %>%
  mutate(
    label = paste0(
      "\u0394=", formatC(mean_delta_eta, format = "f", digits = 4),
      "\nP=", formatC(p_value, format = "f", digits = 3)
    ),
    x = 1.05,
    y = 1.0125
  )

p_a <- ggplot(stage_long, aes(x = phase_label, y = eta)) +
  geom_boxplot(
    aes(colour = model),
    width = 0.48,
    outlier.shape = NA,
    linewidth = 0.32,
    fill = "white",
    alpha = 0.92
  ) +
  geom_point(
    aes(colour = model),
    position = position_jitter(width = 0.075, height = 0, seed = 7),
    size = 1.0,
    alpha = 0.52,
    stroke = 0
  ) +
  geom_line(
    data = stage_summary,
    aes(x = phase_label, y = mean_eta, group = model, colour = model),
    linewidth = 0.45,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = stage_summary,
    aes(x = phase_label, y = mean_eta, colour = model),
    size = 1.45,
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = stage_summary,
    aes(x = phase_label, ymin = mean_eta - ci95_eta, ymax = mean_eta + ci95_eta, colour = model),
    width = 0.13,
    linewidth = 0.32,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = sig_short,
    aes(x = x, y = y, label = label, colour = model),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.0,
    lineheight = 0.88,
    family = "Arial"
  ) +
  facet_wrap(~ model, nrow = 1) +
  scale_colour_manual(values = model_palette) +
  coord_cartesian(ylim = c(0.89, 1.014), clip = "off") +
  labs(
    x = NULL,
    y = expression(eta),
    title = expression("Stage-resolved " * eta * " after component ablation")
  ) +
  guides(colour = "none") +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.spacing.x = unit(4.0, "mm"),
    plot.margin = margin(5, 4, 4, 5)
  )

p_b <- ggplot(variance_calc, aes(x = phase_label, y = var_eta, colour = model, group = model)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.55) +
  scale_colour_manual(values = model_palette) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.00005),
    expand = expansion(mult = c(0.02, 0.14))
  ) +
  labs(
    x = NULL,
    y = expression("Cross-seizure variance of " * eta),
    title = "Stability across seizures"
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = 0.55, size = 1.35))) +
  theme_nature() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.box.margin = margin(-4, 0, 0, 0),
    legend.spacing.x = unit(2.2, "mm"),
    axis.text.x = element_text(size = 6.0)
  )

auc_label_df <- roc_stats %>%
  mutate(
    label = paste0(as.character(model), " AUC=", formatC(auc, format = "f", digits = 3)),
    x = 0.97,
    y = c(0.23, 0.16, 0.09)
  )

p_c <- ggplot(roc_coords, aes(x = fpr, y = tpr, colour = model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.35, colour = "#9A9A9A") +
  geom_step(linewidth = 0.62, direction = "hv") +
  geom_text(
    data = auc_label_df,
    aes(x = x, y = y, label = label, colour = model),
    inherit.aes = FALSE,
    hjust = 1,
    size = 1.95,
    family = "Arial"
  ) +
  scale_colour_manual(values = model_palette, drop = FALSE) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  labs(
    x = "False-positive rate",
    y = "True-positive rate",
    title = "Early-stage discrimination"
  ) +
  theme_nature() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 6.0),
    axis.text.y = element_text(size = 6.0)
  )

fig <- p_a / (p_b | p_c) +
  plot_layout(heights = c(1.42, 1), widths = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 8.2, face = "bold", family = "Arial", colour = "#111111"),
    plot.tag.position = c(0.012, 0.985)
  )

save_pub <- function(plot, filename, width_mm = 183, height_mm = 128, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h, bg = "white")
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial", bg = "white")
  print(plot)
  dev.off()
  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = dpi, background = "white")
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi,
                 compression = "lzw", background = "white")
  print(plot)
  dev.off()
}

base_file <- file.path(out_dir, "Fig7_eta_ablation_components_method_highlight_softer_delta_label_clear")
save_pub(fig, base_file)

export_files <- paste0(base_file, c(".svg", ".pdf", ".png", ".tiff"))
post_checks <- tibble(
  file = basename(export_files),
  exists = file.exists(export_files),
  size_bytes = ifelse(file.exists(export_files), file.info(export_files)$size, NA_real_),
  status = ifelse(file.exists(export_files) & file.info(export_files)$size > 10000, "PASS", "FAIL")
)
write.csv(post_checks, file.path(out_dir, "Fig7_postplot_export_checks_method_highlight_softer_delta_label_clear.csv"), row.names = FALSE)

writeLines(
  c(
    "Fig. 7 method-highlighting note",
    "",
    "Old versions preserved: Fig7_eta_ablation_components.*, Fig7_eta_ablation_components_aligned_palette.*, Fig7_eta_ablation_components_nc_orangered.*, Fig7_eta_ablation_components_method_highlight.* and Fig7_eta_ablation_components_method_highlight_soft.*",
    "Softer method-highlighting delta-label-clear version: Fig7_eta_ablation_components_method_highlight_softer_delta_label_clear.*",
    "Palette: Original #D96661, Coverage only #4778A8, DegreeBias only #7A828C.",
    "Rationale: the main method is foregrounded with a softer coral-red accent, while the ablated variants are visually pushed back with blue and gray. Panel-a labels use the same display content as the current Fig. 8 (Delta and P) and are placed in the upper-left whitespace of each facet to avoid overlap with boxplots or whiskers.",
    "Data and statistics are unchanged from the audited Fig. 7 source tables."
  ),
  file.path(out_dir, "Fig7_palette_alignment_note_method_highlight_softer_delta_label_clear.txt"),
  useBytes = TRUE
)

if (any(post_checks$status == "FAIL")) {
  stop("Export checks failed. Inspect Fig7_postplot_export_checks_method_highlight_softer_delta_label_clear.csv")
}

message("Softer method-highlight delta-label-clear Fig. 7 complete: ", base_file)
