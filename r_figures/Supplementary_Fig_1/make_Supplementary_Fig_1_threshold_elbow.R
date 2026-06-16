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
FIG_HEATMAP_TEXT_PT <- 8.0
FIG_HEATMAP_TEXT_SIZE <- FIG_HEATMAP_TEXT_PT / 2.845276

base_candidates <- Sys.glob("example_project/*0514-/nature_fig/Supplementary_Fig_1")
if (length(base_candidates) < 1) {
  stop("Cannot locate Supplementary_Fig_1 directory under example_project/*0514-/nature_fig")
}
base_dir <- normalizePath(base_candidates[1], winslash = "/", mustWork = TRUE)
out_base <- file.path(base_dir, "Supplementary_Fig_1")

phase_levels <- c("pre-ictal", "early-ictal", "mid-ictal", "late-ictal", "post-ictal")
phase_labels <- c("Pre", "Early", "Mid", "Late", "Post")
phase_cols <- c(
  "pre-ictal" = "#5D83B5",
  "early-ictal" = "#6EAD67",
  "mid-ictal" = "#E6A34A",
  "late-ictal" = "#D96661",
  "post-ictal" = "#8D7BB8"
)
phase_cols_plot <- stats::setNames(unname(phase_cols), phase_labels)

pal <- c(
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

theme_nature <- function(base_size = FIG_TEXT_PT) {
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
theme_set(theme_nature())

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 150, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()
  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = dpi, background = "white")
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi, background = "white", compression = "lzw")
  print(plot)
  dev.off()
}

curve_data <- read_csv(file.path(base_dir, "Supplementary_Fig_1_threshold_curves.csv"), show_col_types = FALSE)
case_summary <- read_csv(file.path(base_dir, "Supplementary_Fig_1_case_summary.csv"), show_col_types = FALSE)
window_density <- read_csv(file.path(base_dir, "Supplementary_Fig_1_window_density_at_elbow.csv"), show_col_types = FALSE)

stopifnot(n_distinct(case_summary$case_id) == 24)
stopifnot(all(c("case_id", "plv_threshold", "mean_pair_density") %in% names(curve_data)))
stopifnot(all(c("case_id", "selected_threshold", "selected_density") %in% names(case_summary)))

case_summary <- case_summary %>%
  arrange(selected_threshold) %>%
  mutate(
    threshold_rank = row_number(),
    case_label = paste0(case_id, "\n", "t = ", number(selected_threshold, accuracy = 0.001))
  )

curve_data <- curve_data %>%
  mutate(case_id = factor(case_id, levels = case_summary$case_id))

cohort_curve <- curve_data %>%
  group_by(plv_threshold) %>%
  summarise(
    mean_density = mean(mean_pair_density, na.rm = TRUE),
    q25_density = quantile(mean_pair_density, 0.25, na.rm = TRUE),
    q75_density = quantile(mean_pair_density, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

selected_points <- case_summary %>%
  select(case_id, selected_threshold, selected_density, threshold_rank) %>%
  mutate(case_id = factor(case_id, levels = levels(curve_data$case_id)))

closest_case <- function(target) {
  case_summary$case_id[which.min(abs(case_summary$selected_threshold - target))]
}
rep_cases <- unique(c(
  case_summary$case_id[1],
  closest_case(stats::median(case_summary$selected_threshold)),
  case_summary$case_id[nrow(case_summary)]
))
rep_roles <- c("Low", "Median", "High")[seq_along(rep_cases)]
rep_label_tbl <- tibble(
  case_id = rep_cases,
  role = rep_roles
) %>%
  left_join(case_summary %>% select(case_id, selected_threshold), by = "case_id") %>%
  mutate(case_label = paste0(role, "\n", "t = ", number(selected_threshold, accuracy = 0.001)))

rep_data <- curve_data %>%
  filter(as.character(case_id) %in% rep_cases) %>%
  left_join(rep_label_tbl %>% select(case_id, case_label), by = "case_id") %>%
  mutate(case_label = factor(case_label, levels = rep_label_tbl$case_label))

rep_points <- selected_points %>%
  filter(as.character(case_id) %in% rep_cases) %>%
  left_join(rep_label_tbl %>% select(case_id, case_label), by = "case_id") %>%
  mutate(case_label = factor(case_label, levels = levels(rep_data$case_label)))

rep_chords <- rep_data %>%
  group_by(case_id, case_label) %>%
  summarise(
    x0 = min(plv_threshold),
    y0 = mean_pair_density[which.min(plv_threshold)],
    x1 = max(plv_threshold),
    y1 = mean_pair_density[which.max(plv_threshold)],
    .groups = "drop"
  )

stage_density <- window_density %>%
  filter(stage %in% phase_levels) %>%
  mutate(stage = factor(stage, levels = phase_levels, labels = phase_labels)) %>%
  group_by(case_id, stage) %>%
  summarise(
    mean_density = mean(pair_density_at_elbow, na.rm = TRUE),
    selected_threshold = first(selected_threshold),
    .groups = "drop"
  )

stage_summary <- stage_density %>%
  group_by(stage) %>%
  summarise(
    n = sum(is.finite(mean_density)),
    sd = sd(mean_density, na.rm = TRUE),
    mean_density = mean(mean_density, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se = sd / sqrt(n),
    ci = qt(0.975, df = n - 1) * se
  )

write_csv(stage_summary, file.path(base_dir, "Supplementary_Fig_1_stage_density_summary.csv"))
write_csv(case_summary, file.path(base_dir, "Supplementary_Fig_1_case_summary_ranked.csv"))
write_csv(
  tibble(
    role = c("low elbow threshold", "median elbow threshold", "high elbow threshold"),
    case_id = rep_cases
  ) %>%
    left_join(case_summary %>% select(case_id, selected_threshold, selected_density, n_channels, n_windows), by = "case_id"),
  file.path(base_dir, "Supplementary_Fig_1_representative_cases.csv")
)

threshold_median <- median(case_summary$selected_threshold, na.rm = TRUE)
threshold_iqr <- quantile(case_summary$selected_threshold, c(0.25, 0.75), na.rm = TRUE)
threshold_range <- range(case_summary$selected_threshold, na.rm = TRUE)

p_a <- ggplot() +
  geom_line(
    data = curve_data,
    aes(plv_threshold, mean_pair_density, group = case_id),
    linewidth = 0.34,
    colour = alpha(pal["neutral_mid"], 0.42)
  ) +
  geom_ribbon(
    data = cohort_curve,
    aes(plv_threshold, ymin = q25_density, ymax = q75_density),
    fill = pal["signal_light"],
    alpha = 0.65,
    colour = NA
  ) +
  geom_line(
    data = cohort_curve,
    aes(plv_threshold, mean_density),
    linewidth = 0.82,
    colour = pal["signal"]
  ) +
  geom_point(
    data = selected_points,
    aes(selected_threshold, selected_density),
    size = 1.9,
    shape = 21,
    stroke = 0.35,
    fill = pal["accent"],
    colour = "white"
  ) +
  annotate(
    "text",
    x = 0.97,
    y = 0.96,
    hjust = 1,
    vjust = 1,
    label = "n = 24 seizures",
    size = FIG_GEOM_TEXT_SIZE,
    colour = pal["neutral_dark"]
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = number_format(accuracy = 0.01)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = number_format(accuracy = 0.01)) +
  labs(
    title = "Seizure-specific density-threshold curves",
    x = "PLV threshold",
    y = "Mean pair-edge density"
  )

p_b <- ggplot(rep_data, aes(plv_threshold, mean_pair_density)) +
  geom_segment(
    data = rep_chords,
    aes(x = x0, y = y0, xend = x1, yend = y1),
    inherit.aes = FALSE,
    linewidth = 0.34,
    linetype = "22",
    colour = pal["neutral_mid"]
  ) +
  geom_line(linewidth = 0.74, colour = pal["signal"]) +
  geom_point(
    data = rep_points,
    aes(selected_threshold, selected_density),
    inherit.aes = FALSE,
    size = 2.1,
    shape = 21,
    stroke = 0.35,
    fill = pal["accent"],
    colour = "white"
  ) +
  facet_wrap(~ case_label, nrow = 1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = number_format(accuracy = 0.1)) +
  labs(
    title = "Representative elbow calls",
    x = "PLV threshold",
    y = "Mean density"
  ) +
  theme(
    strip.text = element_text(size = FIG_HEATMAP_TEXT_PT, face = "bold"),
    axis.text = element_text(size = FIG_HEATMAP_TEXT_PT),
    axis.title = element_text(size = FIG_TEXT_PT),
    panel.spacing.x = unit(3.2, "mm")
  )

p_c <- ggplot(case_summary, aes(x = 1, y = selected_threshold)) +
  geom_violin(
    width = 0.70,
    fill = pal["neutral_fill"],
    colour = pal["neutral_mid"],
    linewidth = 0.35,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    fill = "white",
    colour = pal["ink"],
    linewidth = 0.35
  ) +
  geom_point(
    position = position_jitter(width = 0.055, height = 0, seed = 9),
    shape = 21,
    size = 2.15,
    stroke = 0.35,
    fill = pal["signal"],
    colour = "white"
  ) +
  annotate(
    "text",
    x = 1.39,
    y = threshold_median,
    hjust = 0,
    label = paste0("median\n", number(threshold_median, accuracy = 0.001)),
    size = FIG_GEOM_TEXT_SIZE,
    colour = pal["neutral_dark"]
  ) +
  scale_y_continuous(breaks = seq(0.2, 0.9, 0.2), labels = number_format(accuracy = 0.1)) +
  coord_cartesian(xlim = c(0.62, 1.86), ylim = c(0.2, 0.9), clip = "off") +
  labs(
    title = "Elbow thresholds",
    x = NULL,
    y = "Elbow threshold"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

p_d <- ggplot(stage_density, aes(stage, mean_density, colour = stage, fill = stage)) +
  geom_line(aes(group = case_id), colour = alpha(pal["neutral_mid"], 0.33), linewidth = 0.32) +
  geom_boxplot(
    width = 0.52,
    outlier.shape = NA,
    alpha = 0.22,
    linewidth = 0.35,
    colour = pal["ink"]
  ) +
  geom_point(
    position = position_jitter(width = 0.07, height = 0, seed = 12),
    shape = 21,
    size = 1.95,
    stroke = 0.25,
    colour = "white"
  ) +
  geom_point(
    data = stage_summary,
    aes(stage, mean_density),
    inherit.aes = FALSE,
    shape = 23,
    size = 3.0,
    stroke = 0.45,
    fill = "white",
    colour = pal["ink"]
  ) +
  geom_errorbar(
    data = stage_summary,
    aes(stage, ymin = mean_density - ci, ymax = mean_density + ci),
    inherit.aes = FALSE,
    width = 0.10,
    linewidth = 0.38,
    colour = pal["ink"]
  ) +
  scale_colour_manual(values = phase_cols_plot, guide = "none") +
  scale_fill_manual(values = phase_cols_plot, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = number_format(accuracy = 0.01)) +
  labs(
    title = "Density after applying the seizure-level elbow",
    x = NULL,
    y = "Stage mean density"
  ) +
  theme(plot.margin = margin(5.0, 3.0, 3.0, 8.0))

design <- "
AAB
AAC
DDD
"

fig <- p_a + p_b + p_c + p_d +
  plot_layout(design = design, heights = c(1.05, 1.05, 0.95), widths = c(1.10, 1.10, 1.0), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = FIG_PANEL_PT, face = "bold", colour = pal["ink"]),
    plot.tag.position = c(0.012, 0.985)
  )

save_pub_r(fig, out_base, width_mm = 183, height_mm = 160, dpi = 600)

export_checks <- tibble(
  check = c(
    "seizures represented",
    "threshold grid",
    "selected threshold median",
    "selected threshold IQR",
    "selected threshold range",
    "stage summaries",
    "export formats"
  ),
  value = c(
    as.character(n_distinct(case_summary$case_id)),
    paste0(min(curve_data$plv_threshold), " to ", max(curve_data$plv_threshold), " by 0.0025"),
    number(threshold_median, accuracy = 0.0001),
    paste0(number(threshold_iqr[[1]], accuracy = 0.0001), " to ", number(threshold_iqr[[2]], accuracy = 0.0001)),
    paste0(number(threshold_range[[1]], accuracy = 0.0001), " to ", number(threshold_range[[2]], accuracy = 0.0001)),
    paste(stage_summary$stage, stage_summary$n, sep = " n=", collapse = "; "),
    "SVG, PDF, PNG and TIFF"
  )
)
write_csv(export_checks, file.path(base_dir, "Supplementary_Fig_1_postplot_export_checks.csv"))

legend_text <- paste0(
  "Supplementary Fig. 1 | Threshold selection from seizure-specific density-threshold curves. ",
  "a, Mean pair-edge density as a function of the PLV threshold for each seizure. Grey lines show individual seizures, ",
  "the teal line shows the cohort mean, the shaded band shows the interquartile range, and red points mark the elbow-derived threshold for each seizure. ",
  "b, Representative low-, median- and high-threshold examples showing the distance-to-chord elbow call on the mean density-threshold curve. ",
  "c, Distribution of selected seizure-level PLV thresholds across n = 24 seizures; points denote individual seizures. ",
  "d, Stage-wise mean pair-edge density after applying the same selected threshold to all windows from each seizure. ",
  "Boxes show the median and interquartile range across seizures, lines connect seizure-level stage summaries, and white diamonds with error bars show the mean and 95% confidence interval. ",
  "All source data were generated from 3-s windows advanced in 1-s steps."
)
writeLines(legend_text, con = file.path(base_dir, "Supplementary_Fig_1_legend_draft.txt"), useBytes = TRUE)

description_text <- paste0(
  "To document threshold selection for the PLV-to-hyperedge step, we computed a density-threshold curve for each seizure using the same 3-s windows and 1-s step used in the main analysis. ",
  "A single PLV threshold was selected per seizure from the elbow of the mean density-threshold curve and was then applied to all windows from that seizure. ",
  "Across 24 seizures, the selected thresholds had a median of ", number(threshold_median, accuracy = 0.001),
  " (IQR, ", number(threshold_iqr[[1]], accuracy = 0.001), " to ", number(threshold_iqr[[2]], accuracy = 0.001),
  "; range, ", number(threshold_range[[1]], accuracy = 0.001), " to ", number(threshold_range[[2]], accuracy = 0.001), "). ",
  "The selected thresholds therefore reflect seizure-specific PLV distributions rather than a fixed global cutoff, while preserving a consistent threshold within each seizure for subsequent window-level hyperedge construction."
)
writeLines(description_text, con = file.path(base_dir, "Supplementary_Fig_1_results_description_draft.txt"), useBytes = TRUE)

qa_text <- c(
  "Supplementary Fig. 1 QA notes",
  "Core conclusion: seizure-level PLV thresholds were selected from density-threshold elbows and then applied consistently across windows from the same seizure.",
  "Archetype: quantitative grid with one hero curve panel plus representative, distributional and stage-application diagnostics.",
  "Backend: R / ggplot2 / patchwork only for drawing, preview export and final export.",
  paste0("Source data: ", n_distinct(case_summary$case_id), " seizures; threshold grid ", min(curve_data$plv_threshold), " to ", max(curve_data$plv_threshold), " by 0.0025."),
  paste0("Selected threshold median: ", number(threshold_median, accuracy = 0.0001), "; IQR: ", number(threshold_iqr[[1]], accuracy = 0.0001), "-", number(threshold_iqr[[2]], accuracy = 0.0001), "; range: ", number(threshold_range[[1]], accuracy = 0.0001), "-", number(threshold_range[[2]], accuracy = 0.0001), "."),
  "Statistics: Fig. 1 is a threshold-selection diagnostic; no inferential test is used.",
  "Export check: SVG, PDF, PNG and 600-dpi TIFF produced from the same R object."
)
writeLines(qa_text, con = file.path(base_dir, "Supplementary_Fig_1_QA_notes.txt"), useBytes = TRUE)

message("Supplementary Fig. 1 exported to: ", base_dir)
