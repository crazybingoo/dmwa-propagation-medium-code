library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(readr)
library(scales)
library(svglite)
library(ragg)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
out_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(file_arg[1], winslash = "/", mustWork = TRUE))
} else {
  file.path("F:/6", "\u521d\u7a3f0514-", "nature_fig", "Fig_9")
}
src_dir <- file.path(out_dir, "source_data")

trace_file <- file.path(src_dir, "HR_representative_traces_full.csv")
plv_file <- file.path(src_dir, "HR_representative_plv_full.csv")
if (!file.exists(trace_file)) {
  trace_file <- file.path(src_dir, "HR_representative_traces.csv")
}
if (!file.exists(plv_file)) {
  plv_file <- file.path(src_dir, "HR_representative_plv.csv")
}

trace_df <- read_csv(trace_file, show_col_types = FALSE)
plv_df <- read_csv(plv_file, show_col_types = FALSE)
beta_df <- read_csv(file.path(src_dir, "HR_eta_beta_summary.csv"), show_col_types = FALSE)
trial_df <- read_csv(file.path(src_dir, "HR_eta_trial_metrics.csv"), show_col_types = FALSE)
stats_df <- read_csv(file.path(src_dir, "HR_eta_trend_stats.csv"), show_col_types = FALSE)
mech_df <- read_csv(file.path(src_dir, "mechanism_eta_phase_grid.csv"), show_col_types = FALSE)

theme_set(
  theme_classic(base_size = 6.4, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.30, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      axis.text = element_text(colour = "black", size = 5.7),
      axis.title = element_text(colour = "black", size = 6.2),
      plot.title = element_text(face = "bold", size = 7.0, hjust = 0),
      plot.subtitle = element_text(size = 5.7, colour = "grey28", margin = margin(t = 0.5, b = 1.5)),
      plot.tag = element_text(face = "bold", size = 9.2, colour = "black"),
      plot.tag.position = c(0.005, 0.995),
      legend.title = element_text(size = 5.7),
      legend.text = element_text(size = 5.3),
      legend.key.height = unit(2.5, "mm"),
      legend.key.width = unit(3.4, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 5.8),
      panel.grid = element_blank(),
      plot.margin = margin(2.5, 3, 2.5, 3)
    )
)

blue <- "#2B6CB0"
red <- "#C44E52"
navy <- "#2F4154"
brown <- "#8C6D31"
teal <- "#1F9A8A"
grey_mid <- "#6B7280"
orange <- "#E66101"
gold <- "#E6AB02"
purple <- "#7B3294"
green <- "#66A61E"
heat_low <- "#244C9A"
heat_mid <- "#24A9C7"
heat_high <- "#F2C94C"
phase_low <- "#1D3F8E"
phase_mid <- "#12AFC7"
phase_high <- "#F7CB35"
trace_low <- "#2F4154"
trace_mid <- "#F7F7F2"
trace_high <- "#C44E52"

drive_labeller <- function(beta_label) {
  factor(
    beta_label,
    levels = c("beta=0", "beta=1"),
    labels = c("\u03b2 = 0", "\u03b2 = 1")
  )
}

# ---------- Panel a: HR signal rasters and matched PLV matrices ----------
trace_channels <- sort(unique(trace_df$channel))
trace_plot_df <- trace_df %>%
  filter(channel %in% trace_channels) %>%
  arrange(beta, channel, time) %>%
  group_by(beta_label, channel) %>%
  mutate(
    draw_i = row_number(),
    value_z = pmax(pmin(as.numeric(scale(value)), 2.5), -2.5),
    drive = drive_labeller(beta_label)
  ) %>%
  filter(draw_i %% 2 == 1) %>%
  ungroup()

p_trace <- ggplot(trace_plot_df, aes(time, channel, fill = value_z)) +
  geom_raster(interpolate = FALSE) +
  facet_grid(drive ~ ., switch = "y") +
  scale_fill_gradient2(
    low = trace_low,
    mid = trace_mid,
    high = trace_high,
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = squish,
    name = "z-scored\nx(t)",
    guide = guide_colorbar(
      barheight = unit(14, "mm"),
      barwidth = unit(2.3, "mm"),
      ticks.linewidth = 0.25
    )
  ) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0, 5, 10)) +
  scale_y_continuous(
    breaks = c(1, 9, 18),
    limits = c(0.5, max(trace_channels) + 0.5),
    expand = c(0, 0)
  ) +
  labs(
    title = "Neural-signal simulation",
    subtitle = "Representative HR membrane-potential activity",
    x = "Time (a.u.)",
    y = "Channel",
    tag = "a"
  ) +
  theme(
    legend.position = "right",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, margin = margin(r = 2)),
    plot.tag = element_text(margin = margin(r = 2, b = 0)),
    panel.spacing.y = unit(1.0, "mm")
  )

plv_plot_df <- plv_df %>%
  mutate(drive = drive_labeller(beta_label))

p_plv <- ggplot(plv_plot_df, aes(source, target, fill = plv)) +
  geom_raster() +
  facet_grid(drive ~ .) +
  coord_equal(expand = FALSE) +
  scale_fill_gradientn(
    colours = c("#F8FBFE", "#BAD5EA", "#5E98CB", "#173B6D"),
    limits = c(0.45, 1),
    oob = squish,
    name = "PLV",
    guide = guide_colorbar(
      barheight = unit(15, "mm"),
      barwidth = unit(2.2, "mm"),
      ticks.linewidth = 0.25
    )
  ) +
  scale_x_continuous(breaks = c(1, 9, 18), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(1, 9, 18), expand = c(0, 0)) +
  labs(
    title = "Matched PLV input",
    subtitle = "Same PLV-to-DMW-HLG pipeline",
    x = "Channel",
    y = "Channel"
  ) +
  theme(
    legend.position = "right",
    legend.box.margin = margin(0, -1, 0, -7),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(0.1, "mm"),
    plot.margin = margin(2.5, 0, 2.5, 3),
    panel.spacing.y = unit(1.0, "mm")
  )

# ---------- Panel b: beta scan under the fixed construction ----------
metric_levels <- c("mean_plv", "pair_density", "eta")
metric_labels <- c("Mean PLV", "Pair-edge density", "\u03b7")

beta_long <- beta_df %>%
  transmute(
    beta,
    eta = eta_mean_avg,
    eta_sem = eta_mean_sem,
    mean_plv = mean_plv_mean_avg,
    mean_plv_sem = mean_plv_mean_sem,
    pair_density = pair_density_mean_avg,
    pair_density_sem = pair_density_mean_sem
  ) %>%
  pivot_longer(
    cols = all_of(metric_levels),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    sem = case_when(
      metric == "eta" ~ eta_sem,
      metric == "mean_plv" ~ mean_plv_sem,
      metric == "pair_density" ~ pair_density_sem
    ),
    metric = factor(metric, levels = metric_levels, labels = metric_labels)
  )

trial_long <- trial_df %>%
  transmute(
    beta,
    eta = eta_mean,
    mean_plv = mean_plv_mean,
    pair_density = pair_density_mean
  ) %>%
  pivot_longer(
    cols = all_of(metric_levels),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = metric_levels, labels = metric_labels))

metric_cols <- c("Mean PLV" = blue, "Pair-edge density" = brown, "\u03b7" = navy)

p_beta <- ggplot(beta_long, aes(beta, value, colour = metric, fill = metric)) +
  geom_point(
    data = trial_long,
    aes(beta, value, colour = metric),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.012, height = 0, seed = 42),
    size = 0.55,
    alpha = 0.24
  ) +
  geom_line(linewidth = 0.50) +
  geom_errorbar(aes(ymin = value - sem, ymax = value + sem), width = 0.025, linewidth = 0.30) +
  geom_point(size = 1.45, stroke = 0.22, colour = "white", shape = 21) +
  facet_grid(metric ~ ., scales = "free_y", switch = "y") +
  scale_colour_manual(values = metric_cols, guide = "none") +
  scale_fill_manual(values = metric_cols, guide = "none") +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(-0.03, 1.03)) +
  labs(
    title = "Fixed DMW-HLG readout",
    subtitle = "\u03b7 diverges from synchrony and density",
    x = "Epileptiform drive \u03b2",
    y = NULL,
    tag = "b"
  ) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, margin = margin(r = 2)),
    panel.spacing.y = unit(1.2, "mm")
  )

# ---------- Panels c-d: structural mechanism maps and fixed-parameter slices ----------
slice_levels <- c(0.10, 0.30, 0.50, 0.70, 0.90)
h_targets <- slice_levels
closest_h <- vapply(h_targets, function(target) {
  mech_df$h[which.min(abs(mech_df$h - target))]
}, numeric(1))
h_names <- sprintf("h = %.2f", h_targets)
c_targets <- slice_levels
closest_c <- vapply(c_targets, function(target) {
  mech_df$c[which.min(abs(mech_df$c - target))]
}, numeric(1))
c_names <- sprintf("c = %.2f", c_targets)
slice_cols <- c("#2B6CB0", orange, gold, purple, green)
h_cols <- setNames(slice_cols, h_names)
c_cols <- setNames(slice_cols, c_names)

map_df <- mech_df %>%
  transmute(
    c,
    h,
    eta = eta,
    resource = avgW,
    lambda1 = lambda1
  ) %>%
  pivot_longer(
    cols = c(eta, resource, lambda1),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(metric) %>%
  mutate(value_scaled = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup() %>%
  mutate(
    metric = factor(
      metric,
      levels = c("eta", "resource", "lambda1"),
      labels = c("\u03b7", "R", "\u03bb1")
    )
  )

extract_legend <- function(plot) {
  grob <- ggplotGrob(plot)
  guide_idx <- which(grepl("^guide-box", vapply(grob$grobs, function(x) x$name, character(1))))
  guide_idx <- guide_idx[!vapply(grob$grobs[guide_idx], inherits, logical(1), what = "zeroGrob")]
  if (length(guide_idx) == 0) {
    return(grid::nullGrob())
  }
  grob$grobs[[guide_idx[1]]]
}

p_maps_base <- ggplot(map_df, aes(c, h, fill = value_scaled)) +
  geom_raster(interpolate = FALSE) +
  geom_hline(
    yintercept = closest_h,
    linewidth = 0.16,
    linetype = "11",
    colour = "grey12",
    alpha = 0.35
  ) +
  geom_vline(
    xintercept = closest_c,
    linewidth = 0.16,
    linetype = "11",
    colour = "grey12",
    alpha = 0.35
  ) +
  facet_grid(. ~ metric) +
  coord_equal(expand = FALSE) +
  scale_fill_gradientn(
    colours = c(phase_low, heat_low, phase_mid, phase_high),
    values = rescale(c(0, 0.28, 0.60, 1)),
    limits = c(0, 1),
    oob = squish,
    name = "Within-metric\nscaled value",
    guide = guide_colorbar(
      barheight = unit(17, "mm"),
      barwidth = unit(2.2, "mm"),
      ticks.linewidth = 0.25
    )
  ) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  labs(
    title = "Structural mechanism maps",
    subtitle = "\u03b7, R and \u03bb1 respond differently across the same c-h parameter space",
    x = "Core concentration c",
    y = "Homogenization h",
    tag = "c"
  ) +
  theme(
    legend.position = "right",
    legend.box.spacing = unit(0.4, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    panel.spacing.x = unit(4.0, "mm"),
    strip.text = element_text(face = "bold", size = 6.4),
    axis.text.x = element_text(size = 5.2),
    axis.title.x = element_text(margin = margin(t = 1.0)),
    axis.title.y = element_text(size = 6.2, margin = margin(r = 1.0)),
    plot.margin = margin(2.5, 2, 2.5, 3)
  )

p_maps_legend <- wrap_elements(
  full = extract_legend(
    p_maps_base +
      theme(
        legend.position = "right",
        legend.box.spacing = unit(0, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(0, 0, 0, 0)
      )
  )
)

p_maps <- (
  p_maps_base +
    guides(fill = "none") +
    theme(plot.margin = margin(2.5, 1, 2.5, 3))
) + p_maps_legend + plot_spacer()
p_maps <- p_maps + plot_layout(nrow = 1, widths = c(1, 0.085, 0.23))

fixed_h_df <- bind_rows(lapply(seq_along(closest_h), function(i) {
  mech_df %>%
    filter(abs(h - closest_h[i]) < 1e-9) %>%
    arrange(c) %>%
    mutate(fixed_level = sprintf("%.2f", h_targets[i]))
})) %>%
  select(x = c, fixed_level, resource = avgW, lambda1, eta) %>%
  pivot_longer(
    cols = c(resource, lambda1, eta),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("resource", "lambda1", "eta"),
      labels = c("R", "\u03bb1", "\u03b7")
    ),
    slice_type = factor("Fixed h, varying c", levels = c("Fixed h, varying c", "Fixed c, varying h")),
    fixed_level = factor(fixed_level, levels = sprintf("%.2f", slice_levels))
  )

fixed_c_df <- bind_rows(lapply(seq_along(closest_c), function(i) {
  mech_df %>%
    filter(abs(c - closest_c[i]) < 1e-9) %>%
    arrange(h) %>%
    mutate(fixed_level = sprintf("%.2f", c_targets[i]))
})) %>%
  select(x = h, fixed_level, resource = avgW, lambda1, eta) %>%
  pivot_longer(
    cols = c(resource, lambda1, eta),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("resource", "lambda1", "eta"),
      labels = c("R", "\u03bb1", "\u03b7")
    ),
    slice_type = factor("Fixed c, varying h", levels = c("Fixed h, varying c", "Fixed c, varying h")),
    fixed_level = factor(fixed_level, levels = sprintf("%.2f", slice_levels))
  )

p_slice_df <- bind_rows(fixed_h_df, fixed_c_df)

make_slice_panel <- function(slice_label, metric_label, row_label = NULL, show_title = FALSE, show_x = FALSE) {
  plot_df <- p_slice_df %>%
    filter(slice_type == slice_label, metric == metric_label)

  ggplot(plot_df, aes(x, value, colour = fixed_level)) +
    geom_line(linewidth = 0.42) +
    geom_point(size = 0.70, stroke = 0.22, shape = 21, fill = "white") +
    scale_colour_manual(values = setNames(slice_cols, sprintf("%.2f", slice_levels)), name = "Fixed level") +
    scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
    labs(
      title = if (show_title) as.character(metric_label) else NULL,
      x = if (show_x) "Varied structural parameter" else NULL,
      y = row_label
    ) +
    theme(
      legend.position = "bottom",
      legend.margin = margin(t = -0.5, b = -2),
      legend.key.width = unit(4.6, "mm"),
      legend.spacing.x = unit(1.2, "mm"),
      plot.title = element_text(face = "bold", size = 6.4, hjust = 0.5, margin = margin(b = 1.2)),
      axis.title.x = element_text(size = 6.1, margin = margin(t = 1.0)),
      axis.title.y = element_text(face = "bold", size = 5.9, angle = 0, margin = margin(r = 2.2)),
      axis.text.x = element_text(size = 5.2),
      plot.margin = margin(1.0, 2.3, 1.0, 2.3)
    )
}

p_slice_header <- ggplot() +
  labs(
    title = "Fixed-parameter slices",
    subtitle = "Raw metric values across matched c-h trajectories",
    tag = "d"
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 7.0, hjust = 0, colour = "black"),
    plot.subtitle = element_text(size = 5.7, colour = "grey28", margin = margin(t = 0.5, b = 1.5)),
    plot.tag = element_text(face = "bold", size = 9.2, colour = "black"),
    plot.tag.position = c(0.005, 0.995),
    plot.margin = margin(1.5, 3, 0.5, 3)
  )

p_slice_grid <- (
  make_slice_panel("Fixed h, varying c", "R", "Fixed h,\nvarying c", TRUE, FALSE) +
    make_slice_panel("Fixed h, varying c", "\u03bb1", NULL, TRUE, FALSE) +
    make_slice_panel("Fixed h, varying c", "\u03b7", NULL, TRUE, FALSE)
) / (
  make_slice_panel("Fixed c, varying h", "R", "Fixed c,\nvarying h", FALSE, FALSE) +
    make_slice_panel("Fixed c, varying h", "\u03bb1", NULL, FALSE, TRUE) +
    make_slice_panel("Fixed c, varying h", "\u03b7", NULL, FALSE, FALSE)
)
p_slice_grid <- p_slice_grid +
  plot_layout(guides = "collect", widths = c(1, 1, 1), heights = c(1, 1)) &
  theme(legend.position = "bottom")

p_slices <- p_slice_header / p_slice_grid
p_slices <- p_slices + plot_layout(heights = c(0.06, 1))

top_row <- wrap_plots(p_trace, p_plv, p_beta, nrow = 1, widths = c(1.02, 0.70, 0.86))

fig <- top_row / p_maps / p_slices
fig <- fig + plot_layout(heights = c(1.05, 0.58, 1.08))

save_pub <- function(plot, filename, width_mm = 183, height_mm = 190, dpi = 600) {
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

base <- file.path(out_dir, "Fig9_HR_mechanism_eta_not_density")
save_pub(fig, base)

beta0 <- beta_df %>% filter(beta == min(beta))
beta1 <- beta_df %>% filter(beta == max(beta))
beta_peak <- beta_df %>% slice_max(eta_mean_avg, n = 1, with_ties = FALSE)

eta_rho <- stats_df %>% filter(metric == "eta_mean")
plv_rho <- stats_df %>% filter(metric == "mean_plv_mean")
density_rho <- stats_df %>% filter(metric == "pair_density_mean")

qa_note <- c(
  "Core conclusion: HR neural-signal simulations and structural mechanism simulations jointly support that eta is not reducible to density or synchrony.",
  "Figure layout: top row shows the neural-signal simulation; middle row compares eta, R and lambda1 in the same c-h structural parameter space; bottom row shows fixed-h and fixed-c slices.",
  sprintf(
    "HR beta scan: mean PLV %.3f -> %.3f; pair-edge density %.3f -> %.3f; eta %.3f -> %.3f, peaking at %.3f when beta = %.2f.",
    beta0$mean_plv_mean_avg, beta1$mean_plv_mean_avg,
    beta0$pair_density_mean_avg, beta1$pair_density_mean_avg,
    beta0$eta_mean_avg, beta1$eta_mean_avg,
    beta_peak$eta_mean_avg, beta_peak$beta
  ),
  sprintf(
    "Trend check across 30 trial-level observations: Spearman rho vs beta = %.3f for eta, %.3f for mean PLV and %.3f for pair-edge density.",
    eta_rho$spearman_rho_vs_beta, plv_rho$spearman_rho_vs_beta, density_rho$spearman_rho_vs_beta
  ),
  "Control-variable note: HR signals are converted to PLV input and then passed through the same DMW-HLG construction used for SEEG; the construction itself is not changed.",
  "Structural-map note: eta, R and lambda1 heatmaps use within-metric min-max scaling for pattern comparison; raw ranges are eta 0.765-0.978, R 3.57-17.0 and lambda1 4.15-18.8.",
  "Fixed-slice note: panel d uses five fixed levels, 0.10, 0.30, 0.50, 0.70 and 0.90, for both h and c; curves show raw metric values rather than normalized values.",
  "Review-risk note: the pair-edge density is nearly controlled but shows a small high-beta drift, consistent with PLV ties in the quantile-thresholded construction.",
  "Exports: PNG, SVG, PDF and high-resolution TIFF."
)
writeLines(qa_note, file.path(out_dir, "Fig9_QA_notes.txt"))

message("Saved Fig. 9 exports to: ", out_dir)
