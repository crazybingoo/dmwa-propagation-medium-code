library(ggplot2)
library(patchwork)
library(ggtext)
library(dplyr)
library(tidyr)
library(readr)
library(scales)
library(svglite)
library(ragg)

FIG_TEXT_PT <- 10.5
FIG_PANEL_PT <- 12.5
FIG_DENSE_PT <- 8.0

eta_sym <- "\u03b7"
beta_sym <- "\u03b2"
rho_sym <- "\u03c1"
lambda_sym <- "\u03bb"
eta_md <- paste0("<i>", eta_sym, "</i>")
beta_md <- paste0("<i>", beta_sym, "</i>")
R_md <- "<i>R</i>"
lambda1_md <- paste0("<i>", lambda_sym, "</i><sub>1</sub>")
c_md <- "<i>c</i>"
h_md <- "<i>h</i>"

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
out_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(file_arg[1], winslash = "/", mustWork = TRUE))
} else {
  file.path("figures", "Supplementary_Fig_4")
}
src_dir <- file.path(out_dir, "source_data")

beta_df <- read_csv(file.path(src_dir, "SuppFig4_HR_beta_summary.csv"), show_col_types = FALSE)
trial_df <- read_csv(file.path(src_dir, "SuppFig4_HR_trial_metrics.csv"), show_col_types = FALSE)
stats_df <- read_csv(file.path(src_dir, "SuppFig4_HR_trend_stats.csv"), show_col_types = FALSE)
grid_df <- read_csv(file.path(src_dir, "SuppFig4_mechanism_grid.csv"), show_col_types = FALSE)
summary_df <- read_csv(file.path(src_dir, "SuppFig4_numeric_summary.csv"), show_col_types = FALSE)

theme_set(
  theme_classic(base_size = FIG_TEXT_PT, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.30, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      axis.text = element_text(colour = "black", size = FIG_TEXT_PT),
      axis.title = element_text(colour = "black", size = FIG_TEXT_PT),
      axis.title.y = element_text(margin = margin(r = 1.5)),
      plot.title = ggtext::element_markdown(face = "bold", size = FIG_TEXT_PT, hjust = 0),
      plot.subtitle = ggtext::element_markdown(
        size = FIG_DENSE_PT,
        colour = "grey28",
        margin = margin(t = 0.4, b = 1.2)
      ),
      plot.tag = element_text(face = "bold", size = FIG_PANEL_PT, colour = "black"),
      legend.title = element_text(size = FIG_DENSE_PT),
      legend.text = ggtext::element_markdown(size = FIG_DENSE_PT),
      legend.key.height = unit(2.6, "mm"),
      legend.key.width = unit(5.0, "mm"),
      strip.background = element_blank(),
      strip.text = ggtext::element_markdown(face = "bold", size = FIG_TEXT_PT),
      panel.grid = element_blank(),
      plot.margin = margin(3, 3, 3, 3)
    )
)

blue <- "#2B6CB0"
teal <- "#1F9A8A"
red <- "#C44E52"
navy <- "#2F4154"
orange <- "#E66101"
gold <- "#E6AB02"
purple <- "#7B3294"
green <- "#66A61E"
grey_mid <- "#6B7280"
heat_low <- "#1D3F8E"
heat_mid <- "#12AFC7"
heat_high <- "#F7CB35"

metric_levels <- c("eta", "R", "lambda1")
metric_labels <- c(eta_md, R_md, lambda1_md)
metric_cols <- setNames(c(navy, teal, red), metric_labels)

beta_long <- beta_df %>%
  transmute(
    beta,
    eta = eta_mean_avg,
    eta_sem = eta_mean_sem,
    R = avgW_mean_avg,
    R_sem = avgW_mean_sem,
    lambda1 = lambda1_mean_avg,
    lambda1_sem = lambda1_mean_sem
  ) %>%
  pivot_longer(
    cols = all_of(metric_levels),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    sem = case_when(
      metric == "eta" ~ eta_sem,
      metric == "R" ~ R_sem,
      metric == "lambda1" ~ lambda1_sem
    ),
    metric = factor(metric, levels = metric_levels, labels = metric_labels)
  )

trial_long <- trial_df %>%
  transmute(
    beta,
    eta = eta_mean,
    R = avgW_mean,
    lambda1 = lambda1_mean
  ) %>%
  pivot_longer(
    cols = all_of(metric_levels),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = metric_levels, labels = metric_labels))

trend_labels <- stats_df %>%
  filter(metric %in% c("eta_mean", "avgW_mean", "lambda1_mean")) %>%
  mutate(
    metric = recode(metric, eta_mean = "eta", avgW_mean = "R", lambda1_mean = "lambda1"),
    metric = factor(metric, levels = metric_levels, labels = metric_labels),
    label = sprintf("%s = %.3f; P = %.3g", rho_sym, spearman_rho_vs_beta, p_value)
  )

p_hr <- ggplot(beta_long, aes(beta, value, colour = metric, fill = metric)) +
  geom_point(
    data = trial_long,
    aes(beta, value, colour = metric),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.012, height = 0, seed = 42),
    size = 0.72,
    alpha = 0.28
  ) +
  geom_line(linewidth = 0.52) +
  geom_errorbar(aes(ymin = value - sem, ymax = value + sem), width = 0.025, linewidth = 0.30) +
  geom_point(size = 1.55, shape = 21, stroke = 0.22, colour = "white") +
  geom_text(
    data = trend_labels,
    aes(x = -0.02, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.25,
    size = FIG_DENSE_PT / 2.845276,
    colour = "grey25"
  ) +
  facet_grid(metric ~ ., scales = "free_y", switch = "y", labeller = label_value) +
  scale_colour_manual(values = metric_cols, guide = "none") +
  scale_fill_manual(values = metric_cols, guide = "none") +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(-0.04, 1.04)) +
  scale_y_continuous(n.breaks = 3) +
  labs(
    title = "HR beta scan",
    subtitle = paste0("Trial-level diagnostics of ", eta_md, ", ", R_md, " and ", lambda1_md),
    x = paste0("Epileptiform drive ", beta_md),
    y = NULL,
    tag = "a"
  ) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = ggtext::element_markdown(angle = 0, margin = margin(r = 1.4)),
    axis.title.x = ggtext::element_markdown(size = FIG_TEXT_PT),
    panel.spacing.y = unit(1.2, "mm")
  )

norm_df <- beta_long %>%
  group_by(metric) %>%
  mutate(
    baseline = value[beta == min(beta)][1],
    value_norm = value / baseline,
    sem_norm = sem / abs(baseline)
  ) %>%
  ungroup()

p_norm <- ggplot(norm_df, aes(beta, value_norm, colour = metric, fill = metric)) +
  geom_hline(yintercept = 1, linewidth = 0.25, linetype = "22", colour = "grey55") +
  geom_line(linewidth = 0.58) +
  geom_errorbar(aes(ymin = value_norm - sem_norm, ymax = value_norm + sem_norm), width = 0.025, linewidth = 0.30) +
  geom_point(size = 1.75, shape = 21, stroke = 0.25, colour = "white") +
  scale_colour_manual(values = metric_cols, name = NULL) +
  scale_fill_manual(values = metric_cols, name = NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(-0.04, 1.04)) +
  labs(
    title = "Ratio decomposition",
    subtitle = paste0(R_md, " and ", lambda1_md, " rise together; ", eta_md, " remains bounded"),
    x = paste0("Epileptiform drive ", beta_md),
    y = paste0("Relative to ", beta_md, " = 0"),
    tag = "b"
  ) +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.86), colour = NA),
    axis.title.x = ggtext::element_markdown(size = FIG_TEXT_PT),
    axis.title.y = ggtext::element_markdown(size = FIG_TEXT_PT, margin = margin(r = 1.5)),
    plot.subtitle = ggtext::element_markdown(size = FIG_DENSE_PT, colour = "grey28")
  )

map_df <- grid_df %>%
  transmute(c, h, eta, R, lambda1) %>%
  pivot_longer(
    cols = c(eta, R, lambda1),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(metric) %>%
  mutate(value_scaled = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup() %>%
  mutate(metric = factor(metric, levels = metric_levels, labels = metric_labels))

p_maps <- ggplot(map_df, aes(c, h, fill = value_scaled)) +
  geom_raster(interpolate = FALSE) +
  facet_grid(. ~ metric, labeller = label_value) +
  coord_equal(expand = FALSE) +
  scale_fill_gradientn(
    colours = c(heat_low, "#244C9A", heat_mid, heat_high),
    values = rescale(c(0, 0.28, 0.62, 1)),
    limits = c(0, 1),
    oob = squish,
    guide = "none"
  ) +
  scale_x_continuous(breaks = c(0, 1), labels = c("0", "1"), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), expand = c(0, 0)) +
  labs(
    title = "Structural mechanism grid",
    subtitle = paste0("Matched ", c_md, "-", h_md, " simulations; colours are scaled within each metric"),
    x = paste0("Core concentration ", c_md),
    y = paste0("Homogenization ", h_md),
    tag = "c"
  ) +
  theme(
    panel.spacing.x = unit(4.0, "mm"),
    axis.title.x = ggtext::element_markdown(size = FIG_TEXT_PT, margin = margin(t = 1.0)),
    axis.title.y = ggtext::element_markdown(size = FIG_TEXT_PT, margin = margin(r = 0.8)),
    strip.text = ggtext::element_markdown(face = "bold", size = FIG_TEXT_PT)
  )

cbar_df <- tibble(x = seq(0, 1, length.out = 200), y = 1, value_scaled = x)
p_cbar <- ggplot(cbar_df, aes(x, y, fill = value_scaled)) +
  geom_raster() +
  scale_fill_gradientn(
    colours = c(heat_low, "#244C9A", heat_mid, heat_high),
    values = rescale(c(0, 0.28, 0.62, 1)),
    limits = c(0, 1),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = sprintf("%.2f", c(0, 0.25, 0.50, 0.75, 1.00)),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Within-metric scaled value", y = NULL) +
  theme_void(base_family = "Arial") +
  theme(
    axis.text.x = element_text(size = FIG_DENSE_PT, colour = "black", margin = margin(t = 1)),
    axis.ticks.x = element_line(linewidth = 0.20, colour = "white"),
    axis.title.x = element_text(size = FIG_DENSE_PT, colour = "black", margin = margin(t = 2)),
    plot.margin = margin(-1, 0, 1, 0)
  )

p_maps <- p_maps / p_cbar + plot_layout(heights = c(1, 0.16))

slice_levels <- c(0.10, 0.30, 0.50, 0.70, 0.90)
slice_cols <- c(blue, orange, gold, purple, green)
nearest <- function(x, target) x[which.min(abs(x - target))]

fixed_h <- sort(unique(grid_df$h))
fixed_c <- sort(unique(grid_df$c))
closest_h <- vapply(slice_levels, function(v) nearest(fixed_h, v), numeric(1))
closest_c <- vapply(slice_levels, function(v) nearest(fixed_c, v), numeric(1))

fixed_h_df <- bind_rows(lapply(seq_along(closest_h), function(i) {
  grid_df %>%
    filter(abs(h - closest_h[i]) < 1e-9) %>%
    transmute(
      x = c,
      fixed_level = sprintf("%.2f", slice_levels[i]),
      slice_type = "Fixed h",
      eta, R, lambda1
    )
}))

fixed_c_df <- bind_rows(lapply(seq_along(closest_c), function(i) {
  grid_df %>%
    filter(abs(c - closest_c[i]) < 1e-9) %>%
    transmute(
      x = h,
      fixed_level = sprintf("%.2f", slice_levels[i]),
      slice_type = "Fixed c",
      eta, R, lambda1
    )
}))

slice_df <- bind_rows(fixed_h_df, fixed_c_df) %>%
  pivot_longer(cols = c(eta, R, lambda1), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = metric_levels, labels = metric_labels),
    fixed_level = factor(fixed_level, levels = sprintf("%.2f", slice_levels)),
    slice_type = factor(
      slice_type,
      levels = c("Fixed h", "Fixed c")
    )
  )

p_slices <- ggplot(slice_df, aes(x, value, colour = fixed_level, group = fixed_level)) +
  geom_line(linewidth = 0.42) +
  geom_point(size = 0.68, stroke = 0.20, shape = 21, fill = "white") +
  facet_grid(slice_type ~ metric, scales = "free_y", labeller = label_value) +
  scale_colour_manual(values = setNames(slice_cols, sprintf("%.2f", slice_levels)), name = "Fixed level") +
  scale_x_continuous(breaks = c(0, 1), limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "Fixed-parameter slices",
    subtitle = paste0("Raw values across matched ", c_md, "-", h_md, " trajectories; top varies ", c_md, ", bottom varies ", h_md),
    x = "Varied structural parameter",
    y = "Raw metric value",
    tag = "d"
  ) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -1, b = -2),
    legend.key.width = unit(4.8, "mm"),
    strip.text.x = ggtext::element_markdown(face = "bold", size = FIG_TEXT_PT),
    strip.text.y.right = element_text(face = "bold", angle = 270, size = FIG_DENSE_PT, margin = margin(l = 1.5)),
    panel.spacing.x = unit(5.0, "mm"),
    panel.spacing.y = unit(1.4, "mm")
  )

top_row <- p_hr + p_norm + plot_layout(widths = c(1.25, 0.95))
fig <- top_row / p_maps / p_slices
fig <- fig + plot_layout(heights = c(1.00, 0.92, 1.18))

save_pub <- function(plot, filename, width_mm = 183, height_mm = 178, dpi = 600) {
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

  ragg::agg_tiff(
    paste0(filename, ".tiff"),
    width = w,
    height = h,
    units = "in",
    res = dpi,
    background = "white",
    compression = "lzw"
  )
  print(plot)
  dev.off()
}

base <- file.path(out_dir, "Supplementary_Fig_4_simulation_diagnostics")
save_pub(fig, base)

eta_trend <- stats_df %>% filter(metric == "eta_mean")
R_trend <- stats_df %>% filter(metric == "avgW_mean")
lambda_trend <- stats_df %>% filter(metric == "lambda1_mean")
range_row <- function(metric_name) {
  summary_df %>% filter(source == "structural_grid", metric == metric_name) %>% slice(1)
}
eta_range <- range_row("eta")
R_range <- range_row("R")
lambda_range <- range_row("lambda1")

legend_text <- c(
  paste0(
    "Supplementary Fig. 4 | Simulation diagnostics for ", eta_sym, ", R and ", lambda_sym, "1."
  ),
  paste0(
    "a, Hindmarsh-Rose signal simulations were passed through the same DMW-HLG construction and ",
    eta_sym, " readout used for the empirical SEEG analyses. Points denote trial-level means, and lines show mean +/- s.e.m. across six trials per beta level."
  ),
  paste0(
    "b, Normalized decomposition of ", eta_sym, " = R/", lambda_sym, "1 across beta levels. Trial-level Spearman correlations against beta were ",
    eta_sym, ": rho = ", sprintf("%.3f", eta_trend$spearman_rho_vs_beta), ", P = ", sprintf("%.3g", eta_trend$p_value),
    "; R: rho = ", sprintf("%.3f", R_trend$spearman_rho_vs_beta), ", P = ", sprintf("%.3g", R_trend$p_value),
    "; and ", lambda_sym, "1: rho = ", sprintf("%.3f", lambda_trend$spearman_rho_vs_beta), ", P = ", sprintf("%.3g", lambda_trend$p_value), "."
  ),
  paste0(
    "c, Structural mechanism grid showing within-metric scaled values of ", eta_sym, ", R and ", lambda_sym, "1 across core concentration c and homogenization h."
  ),
  paste0(
    "d, Fixed-parameter slices from the same grid show raw metric values along matched c-h trajectories. The grid contained 41 x 41 settings and 30 Monte Carlo realizations per setting; raw ranges were ",
    eta_sym, " = ", sprintf("%.3f", eta_range$min_or_rho), "-", sprintf("%.3f", eta_range$max_or_p),
    ", R = ", sprintf("%.2f", R_range$min_or_rho), "-", sprintf("%.2f", R_range$max_or_p),
    " and ", lambda_sym, "1 = ", sprintf("%.2f", lambda_range$min_or_rho), "-", sprintf("%.2f", lambda_range$max_or_p), "."
  )
)
writeLines(legend_text, file.path(out_dir, "Supplementary_Fig_4_legend_draft.txt"), useBytes = TRUE)

qa_note <- c(
  "Core conclusion: simulation diagnostics show that eta is a bounded ratio between the resource term R and the leading spectral term lambda1, rather than a simple monotonic surrogate of either term.",
  "Evidence hierarchy: panel a shows HR trial-level raw metrics; panel b shows the ratio decomposition; panel c shows the structural grid; panel d shows raw fixed-parameter slices.",
  "Statistics: HR beta scan uses six trials per beta level, with trial-level Spearman correlations against beta reported in the legend.",
  "Source data: HR summaries were exported from the native MATLAB DMW-HLG eta simulation; structural maps were exported from eta_phase_mechanism_data.mat.",
  "Export check: PNG, SVG, PDF and TIFF were generated by R only."
)
writeLines(qa_note, file.path(out_dir, "Supplementary_Fig_4_QA_notes.txt"), useBytes = TRUE)

exports <- tibble(
  file = paste0(base, c(".png", ".svg", ".pdf", ".tiff")),
  exists = file.exists(file),
  bytes = if_else(exists, file.info(file)$size, as.numeric(NA))
)
write_csv(exports, file.path(out_dir, "Supplementary_Fig_4_postplot_export_checks.csv"))
message("Saved Supplementary Fig. 4 exports to: ", out_dir)
