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

base_candidates <- Sys.glob("example_project/*0514-/nature_fig/Supplementary_Fig_2")
if (length(base_candidates) < 1) {
  stop("Cannot locate Supplementary_Fig_2 directory under example_project/*0514-/nature_fig")
}
base_dir <- normalizePath(base_candidates[1], winslash = "/", mustWork = TRUE)
out_base <- file.path(base_dir, "Supplementary_Fig_2")

window_csv <- file.path(base_dir, "Supplementary_Fig_2_window_eta_source.csv")
stage_csv <- file.path(base_dir, "Supplementary_Fig_2_stage_eta_source.csv")
if (!file.exists(window_csv) || !file.exists(stage_csv)) {
  stop("Missing MATLAB-exported source CSV files for Supplementary Fig. 2")
}

phase_levels <- c("pre-ictal", "early", "mid", "late", "post-ictal")
phase_labels <- c("Pre", "Early", "Mid", "Late", "Post")
phase_cols <- c(
  "Pre" = "#5D83B5",
  "Early" = "#6EAD67",
  "Mid" = "#E6A34A",
  "Late" = "#D96661",
  "Post" = "#8D7BB8"
)

pal <- c(
  ink = "#202124",
  neutral_dark = "#4A4F55",
  neutral_mid = "#9AA3AD",
  neutral_light = "#D9DEE7",
  grid = "#ECEFF3"
)

theme_nature <- function(base_size = FIG_TEXT_PT) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = pal["ink"]),
      axis.ticks = element_line(linewidth = 0.30, colour = pal["ink"]),
      axis.text = element_text(size = FIG_TEXT_PT, colour = pal["ink"]),
      axis.title = element_text(size = FIG_TEXT_PT, colour = pal["ink"]),
      axis.title.y = element_text(margin = margin(r = 1.5)),
      legend.title = element_text(size = FIG_TEXT_PT, colour = pal["ink"]),
      legend.text = element_text(size = FIG_TEXT_PT, colour = pal["ink"]),
      legend.key.height = unit(3.0, "mm"),
      legend.key.width = unit(5.5, "mm"),
      legend.background = element_blank(),
      panel.grid.major = element_line(linewidth = 0.18, colour = pal["grid"]),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", colour = NA),
      strip.text = element_text(size = FIG_TEXT_PT, face = "bold", colour = pal["ink"]),
      plot.margin = margin(3.0, 4.8, 3.0, 3.0)
    )
}
theme_set(theme_nature())

save_pub_r <- function(plot, filename, width_mm = 245, height_mm = 160, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot)
  dev.off()

  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()

  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = dpi)
  print(plot)
  dev.off()

  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()
}

window_raw <- read_csv(window_csv, show_col_types = FALSE)
stage_raw <- read_csv(stage_csv, show_col_types = FALSE)

case_order <- stage_raw %>%
  distinct(case_id) %>%
  mutate(
    seizure_index = row_number(),
    seizure_label = sprintf("Seizure %02d", seizure_index)
  )

plot_df <- window_raw %>%
  select(case_id, window_idx, phase5, eta, num_hyperedges) %>%
  mutate(
    phase5 = factor(phase5, levels = phase_levels, labels = phase_labels)
  ) %>%
  left_join(case_order, by = "case_id") %>%
  group_by(case_id) %>%
  arrange(window_idx, .by_group = TRUE) %>%
  mutate(
    rel_window = if (max(window_idx) > min(window_idx)) {
      100 * (window_idx - min(window_idx)) / (max(window_idx) - min(window_idx))
    } else {
      0
    }
  ) %>%
  ungroup() %>%
  mutate(seizure_label = factor(seizure_label, levels = case_order$seizure_label))

baseline_df <- plot_df %>%
  filter(phase5 == "Pre") %>%
  group_by(seizure_label) %>%
  summarise(pre_eta = mean(eta, na.rm = TRUE), .groups = "drop")

case_key <- case_order %>%
  select(seizure_index, seizure_label, case_id)

clean_source <- plot_df %>%
  select(seizure_index, seizure_label, case_id, window_idx, rel_window, phase5, eta, num_hyperedges)

stage_summary <- stage_raw %>%
  mutate(phase5 = factor(phase5, levels = phase_levels, labels = phase_labels)) %>%
  group_by(phase5) %>%
  summarise(
    n_seizures = n(),
    mean_eta = mean(eta_mean, na.rm = TRUE),
    sd_eta = sd(eta_mean, na.rm = TRUE),
    median_eta = median(eta_mean, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(case_key, file.path(base_dir, "Supplementary_Fig_2_case_key.csv"))
write_csv(clean_source, file.path(base_dir, "Supplementary_Fig_2_plot_source.csv"))
write_csv(stage_summary, file.path(base_dir, "Supplementary_Fig_2_stage_summary.csv"))

eta_range <- range(plot_df$eta, na.rm = TRUE)
eta_limits <- c(floor(eta_range[1] * 100) / 100, ceiling(eta_range[2] * 100) / 100)
eta_breaks <- pretty(eta_limits, n = 4)
eta_breaks <- eta_breaks[eta_breaks >= eta_limits[1] & eta_breaks <= eta_limits[2]]

trajectory_plot <- ggplot(plot_df, aes(x = rel_window, y = eta)) +
  geom_hline(
    data = baseline_df,
    aes(yintercept = pre_eta),
    inherit.aes = FALSE,
    linewidth = 0.26,
    linetype = "22",
    colour = pal["neutral_mid"]
  ) +
  geom_line(aes(colour = phase5, group = seizure_label), linewidth = 0.34, alpha = 0.92, lineend = "round") +
  facet_wrap(~seizure_label, ncol = 6) +
  scale_colour_manual(values = phase_cols, drop = FALSE, name = "Stage") +
  scale_x_continuous(
    breaks = c(0, 50, 100),
    labels = c("0", "50", "100"),
    expand = expansion(mult = c(0.018, 0.045))
  ) +
  scale_y_continuous(
    breaks = eta_breaks,
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.04, 0.05))
  ) +
  coord_cartesian(ylim = eta_limits, clip = "on") +
  labs(
    x = "Normalized window position (%)",
    y = expression(italic(eta))
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.1, alpha = 1), nrow = 1)) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.spacing = unit(1.0, "mm"),
    panel.spacing.x = unit(3.1, "mm"),
    panel.spacing.y = unit(3.2, "mm")
  )

save_pub_r(trajectory_plot, out_base, width_mm = 245, height_mm = 160, dpi = 600)

export_files <- paste0(out_base, c(".png", ".svg", ".pdf", ".tiff"))
export_checks <- tibble(
  file = basename(export_files),
  exists = file.exists(export_files),
  bytes = if_else(file.exists(export_files), as.numeric(file.info(export_files)$size), NA_real_)
)
write_csv(export_checks, file.path(base_dir, "Supplementary_Fig_2_postplot_export_checks.csv"))

legend_text <- paste(
  "Supplementary Fig. 2 | Seizure-level \u03b7 trajectories.",
  "Window-resolved effective refractive index (\u03b7) trajectories are shown for all 24 seizures.",
  "Each panel denotes one seizure, with the x axis normalized to the window sequence of that recording.",
  "Line colour indicates the five seizure stages: pre-ictal, early ictal, mid-ictal, late ictal and post-ictal.",
  "Dashed grey lines mark the corresponding pre-ictal mean \u03b7 for each seizure.",
  "Source data are window-level \u03b7 estimates derived from 3-s windows after applying one seizure-level PLV threshold selected from the density-threshold elbow.",
  sep = " "
)
writeLines(legend_text, file.path(base_dir, "Supplementary_Fig_2_legend_draft.txt"), useBytes = TRUE)

results_text <- paste(
  "Window-level \u03b7 trajectories were available for all 24 seizures and showed seizure-specific temporal profiles across pre-ictal, ictal and post-ictal windows.",
  "These traces provide the seizure-level basis for the stage-wise cohort summaries reported in the main text.",
  sep = " "
)
writeLines(results_text, file.path(base_dir, "Supplementary_Fig_2_results_description_draft.txt"), useBytes = TRUE)

qa_text <- c(
  "Supplementary Fig. 2 QA notes",
  sprintf("Core conclusion: all %d seizures have traceable window-level eta trajectories supporting the stage-wise eta summaries.", n_distinct(plot_df$case_id)),
  "Archetype: quantitative grid.",
  "Backend: R only (ggplot2 + svglite/cairo_pdf/ragg).",
  sprintf("Source rows: %d window-level observations; %d seizure-stage rows.", nrow(window_raw), nrow(stage_raw)),
  sprintf("Eta range displayed: %.4f to %.4f; observed range: %.4f to %.4f.", eta_limits[1], eta_limits[2], eta_range[1], eta_range[2]),
  "Visual checks to perform: facet labels readable, eta symbol italic, phase legend clear, no line/legend overlap."
)
writeLines(qa_text, file.path(base_dir, "Supplementary_Fig_2_QA_notes.txt"), useBytes = TRUE)

message("Supplementary Fig. 2 exports written to: ", base_dir)
