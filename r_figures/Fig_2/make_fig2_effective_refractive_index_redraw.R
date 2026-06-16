suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggbeeswarm)
  library(R.matlab)
  library(svglite)
  library(ragg)
  library(grid)
  library(scales)
})

FIG_TEXT_PT <- 6.2
FIG_PANEL_PT <- 8.0
FIG_GEOM_TEXT_SIZE <- FIG_TEXT_PT / 2.845276

base_candidates <- Sys.glob("example_project/*0514-/nature_fig/Fig_2")
if (length(base_candidates) < 1) {
  stop("Cannot locate Fig_2 directory under example_project/*0514-/nature_fig/Fig_2")
}
base_dir <- normalizePath(base_candidates[1], winslash = "/", mustWork = TRUE)
out_base <- file.path(base_dir, "Fig2_effective_refractive_index_increases_redraw")

phase_levels <- c("pre-ictal", "early", "mid", "late", "post-ictal")
phase_labels <- c("Pre", "Early", "Mid", "Late", "Post")
effect_levels <- phase_levels[-1]
effect_labels <- phase_labels[-1]

phase_cols <- c(
  "pre-ictal" = "#5D83B5",
  "early" = "#6EAD67",
  "mid" = "#E6A34A",
  "late" = "#D96661",
  "post-ictal" = "#8D7BB8"
)
phase3_levels <- c("pre-ictal", "ictal", "post-ictal")
phase3_labels <- c("Pre", "Ictal", "Post")
phase3_cols <- c(
  "pre-ictal" = phase_cols[["pre-ictal"]],
  "ictal" = "#D98245",
  "post-ictal" = phase_cols[["post-ictal"]]
)

eta_sym <- "\u03b7"
delta_sym <- "\u0394"
density_sym <- "\u03b4\u00b2"

theme_nature <- function(base_size = FIG_TEXT_PT) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "#202124"),
      axis.ticks = element_line(linewidth = 0.28, colour = "#202124"),
      axis.text = element_text(size = base_size, colour = "#202124"),
      axis.title = element_text(size = base_size, colour = "#202124"),
      plot.title = element_text(size = base_size, face = "bold", hjust = 0),
      legend.position = "none",
      panel.grid.major = element_line(linewidth = 0.16, colour = "#ECEFF3"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(2.8, 3.2, 2.8, 3.2)
    )
}
theme_set(theme_nature())

stars_from_p <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", "")
      )
    )
  )
}

format_p <- function(p) {
  ifelse(p < 0.001, "P<0.001", paste0("P=", signif(p, 2)))
}

ci_one_sample <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || stats::sd(x) == 0) {
    return(c(mean(x), mean(x)))
  }
  stats::t.test(x)$conf.int
}

eta_stage <- read_csv(file.path(base_dir, "ALL_CASES_stage5_eta.csv"), show_col_types = FALSE) %>%
  mutate(
    phase5 = factor(phase5, levels = phase_levels),
    phase_label = factor(phase_labels[as.integer(phase5)], levels = phase_labels)
  )

eta_window <- read_csv(file.path(base_dir, "seizure_01_window_level_eta.csv"), show_col_types = FALSE) %>%
  mutate(
    phase5 = factor(phase5, levels = phase_levels),
    phase_label = factor(phase_labels[as.integer(phase5)], levels = phase_labels)
  )

stage_summary <- eta_stage %>%
  group_by(phase5, phase_label) %>%
  summarise(
    n = n(),
    mean = mean(eta_mean, na.rm = TRUE),
    sd = sd(eta_mean, na.rm = TRUE),
    se = sd / sqrt(n),
    ci = qt(0.975, df = n - 1) * se,
    .groups = "drop"
  )

wide_eta <- eta_stage %>%
  select(case_id, phase5, eta_mean) %>%
  pivot_wider(names_from = phase5, values_from = eta_mean)

delta_long <- wide_eta %>%
  select(case_id, all_of(phase_levels)) %>%
  pivot_longer(
    cols = all_of(effect_levels),
    names_to = "phase5",
    values_to = "eta_stage"
  ) %>%
  mutate(
    delta_eta = eta_stage - `pre-ictal`,
    phase5 = factor(phase5, levels = effect_levels),
    phase_label = factor(effect_labels[as.integer(phase5)], levels = rev(effect_labels))
  )

paired_stats <- delta_long %>%
  group_by(phase5, phase_label) %>%
  summarise(
    n = sum(is.finite(delta_eta)),
    mean_diff = mean(delta_eta, na.rm = TRUE),
    median_diff = median(delta_eta, na.rm = TRUE),
    ci_low = ci_one_sample(delta_eta)[1],
    ci_high = ci_one_sample(delta_eta)[2],
    p_wilcox = suppressWarnings(wilcox.test(delta_eta, mu = 0, exact = FALSE)$p.value),
    dz = mean(delta_eta, na.rm = TRUE) / sd(delta_eta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    p_holm = p.adjust(p_wilcox, method = "holm"),
    stars = stars_from_p(p_holm),
    p_label = paste0(stars, "  ", format_p(p_holm))
  )

write_csv(stage_summary, file.path(base_dir, "Fig2_redraw_stage_summary.csv"))
write_csv(paired_stats, file.path(base_dir, "Fig2_redraw_paired_statistics.csv"))

stage_bounds <- eta_window %>%
  filter(!is.na(phase5)) %>%
  group_by(phase5, phase_label) %>%
  summarise(
    xmin = min(window_idx),
    xmax = max(window_idx),
    xmid = (xmin + xmax) / 2,
    .groups = "drop"
  )

transition_x <- stage_bounds %>%
  arrange(factor(phase5, levels = phase_levels)) %>%
  slice(seq_len(nrow(stage_bounds) - 1)) %>%
  pull(xmax)

eta_ymax <- max(eta_window$eta, na.rm = TRUE)
eta_ymin <- min(eta_window$eta, na.rm = TRUE)
eta_pad <- 0.10 * (eta_ymax - eta_ymin)

p_a <- ggplot(eta_window, aes(window_idx, eta)) +
  geom_rect(
    data = stage_bounds,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = phase5),
    inherit.aes = FALSE, alpha = 0.12, colour = NA
  ) +
  geom_line(linewidth = 0.30, colour = "#202124") +
  geom_vline(xintercept = transition_x, linetype = "22", linewidth = 0.28, colour = "#8A8F98") +
  geom_text(
    data = stage_bounds,
    aes(x = xmid, y = eta_ymax + eta_pad * 0.34, label = phase_label),
    inherit.aes = FALSE, size = FIG_GEOM_TEXT_SIZE, family = "Arial",
    fontface = "bold", colour = "#202124"
  ) +
  scale_fill_manual(values = phase_cols) +
  coord_cartesian(ylim = c(eta_ymin - eta_pad, eta_ymax + eta_pad * 0.70), clip = "off") +
  labs(
    title = paste0("Representative seizure: time-resolved ", eta_sym),
    x = "Time windows",
    y = paste0("Effective refractive index, ", eta_sym)
  )

p_b <- ggplot(eta_stage, aes(phase_label, eta_mean, group = case_id)) +
  geom_line(colour = "#C8CED7", linewidth = 0.22, alpha = 0.24) +
  geom_point(aes(fill = phase5), shape = 21, size = 0.95, stroke = 0.08, colour = "white", alpha = 0.58) +
  geom_errorbar(
    data = stage_summary,
    aes(phase_label, ymin = mean - ci, ymax = mean + ci),
    inherit.aes = FALSE, width = 0.11, linewidth = 0.30, colour = "#202124"
  ) +
  geom_line(
    data = stage_summary,
    aes(phase_label, mean, group = 1),
    inherit.aes = FALSE, colour = "#202124", linewidth = 0.48
  ) +
  geom_point(
    data = stage_summary,
    aes(phase_label, mean, fill = phase5),
    inherit.aes = FALSE, shape = 21, size = 2.30, stroke = 0.34, colour = "#202124"
  ) +
  scale_fill_manual(values = phase_cols) +
  labs(
    title = paste0("Stage-wise ", eta_sym, " across seizures"),
    x = NULL,
    y = paste0("Stage mean ", eta_sym)
  )

delta_xlim <- range(delta_long$delta_eta, paired_stats$ci_low, paired_stats$ci_high, na.rm = TRUE)
delta_pad <- diff(delta_xlim) * 0.24
label_x <- delta_xlim[2] + delta_pad * 0.55

p_c <- ggplot(delta_long, aes(delta_eta, phase_label, fill = phase5)) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.28, colour = "#7A828C") +
  geom_point(
    position = position_jitter(width = 0, height = 0.08, seed = 4),
    shape = 21, size = 1.18, stroke = 0.08, colour = "white", alpha = 0.58
  ) +
  geom_segment(
    data = paired_stats,
    aes(x = ci_low, xend = ci_high, y = phase_label, yend = phase_label),
    inherit.aes = FALSE, linewidth = 0.36, colour = "#202124"
  ) +
  geom_point(
    data = paired_stats,
    aes(mean_diff, phase_label, fill = phase5),
    inherit.aes = FALSE, shape = 23, size = 2.35, stroke = 0.34, colour = "#202124"
  ) +
  geom_text(
    data = paired_stats,
    aes(x = label_x, y = phase_label, label = p_label),
    inherit.aes = FALSE, family = "Arial", size = FIG_GEOM_TEXT_SIZE * 0.92,
    hjust = 0, colour = "#202124"
  ) +
  scale_fill_manual(values = phase_cols[effect_levels]) +
  coord_cartesian(xlim = c(delta_xlim[1] - delta_pad * 0.25, label_x + delta_pad * 0.95), clip = "off") +
  labs(
    title = paste0("Paired effect relative to Pre"),
    x = paste0(delta_sym, eta_sym, " relative to Pre"),
    y = NULL
  ) +
  theme(
    plot.margin = margin(2.8, 9.5, 2.8, 3.2)
  )

mat_dir <- file.path(base_dir, "seizure_data")
mat_files <- list.files(mat_dir, pattern = "\\.mat$", full.names = TRUE)
exclude_files <- c("example_excluded_state_01", "example_excluded_state_02",
                   "example_excluded_state_03", "example_excluded_state_04")

read_stage_mat <- function(path) {
  name <- tools::file_path_sans_ext(basename(path))
  if (name %in% exclude_files) return(NULL)
  obj <- readMat(path)
  pid <- strsplit(name, "_", fixed = TRUE)[[1]][1]
  pre <- as.numeric(obj[["stage.pre"]][1, ])
  ictal <- as.numeric(obj[["stage.ictal"]][1, ])
  post <- as.numeric(obj[["stage.post"]][1, ])
  tibble(
    file_id = name,
    patient = pid,
    phase3 = factor(phase3_levels, levels = phase3_levels),
    density = c(pre[1], ictal[1], post[1]),
    eta = c(pre[2], ictal[2], post[2])
  )
}

state_raw <- bind_rows(lapply(mat_files, read_stage_mat))
if (nrow(state_raw) < 1) {
  stop("No state-space .mat files were available for panel d")
}

state_norm <- state_raw %>%
  group_by(patient) %>%
  mutate(
    density_sd = sd(density, na.rm = TRUE),
    eta_sd = sd(eta, na.rm = TRUE),
    density_z = (density - mean(density, na.rm = TRUE)) / ifelse(density_sd == 0, 1, density_sd),
    eta_z = (eta - mean(eta, na.rm = TRUE)) / ifelse(eta_sd == 0, 1, eta_sd)
  ) %>%
  ungroup() %>%
  select(-density_sd, -eta_sd) %>%
  arrange(file_id, phase3)

state_centres <- state_norm %>%
  group_by(phase3) %>%
  summarise(
    density_z = mean(density_z, na.rm = TRUE),
    eta_z = mean(eta_z, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    phase3 = factor(phase3, levels = phase3_levels),
    phase_label = factor(phase3_labels[as.integer(phase3)], levels = phase3_labels)
  ) %>%
  arrange(phase3)

state_labels <- state_centres %>%
  mutate(
    label_x = case_when(
      phase3 == "pre-ictal" ~ density_z - 0.30,
      phase3 == "ictal" ~ density_z + 0.60,
      TRUE ~ density_z - 0.78
    ),
    label_y = case_when(
      phase3 == "pre-ictal" ~ eta_z + 0.24,
      phase3 == "ictal" ~ eta_z + 0.34,
      TRUE ~ eta_z + 0.08
    )
  )

state_xlim <- range(state_norm$density_z, na.rm = TRUE)
state_ylim <- range(state_norm$eta_z, na.rm = TRUE)
state_pad <- max(diff(state_xlim), diff(state_ylim)) * 0.09
state_xlim <- state_xlim + c(-state_pad, state_pad)
state_ylim <- state_ylim + c(-state_pad, state_pad)

write_csv(state_norm, file.path(base_dir, "Fig2_redraw_state_space_points.csv"))

p_d <- ggplot(state_norm, aes(density_z, eta_z)) +
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.24, colour = "#A8AFB8") +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.24, colour = "#A8AFB8") +
  stat_ellipse(
    aes(colour = phase3),
    type = "norm", level = 0.68, linewidth = 0.30, alpha = 0.70,
    show.legend = FALSE
  ) +
  geom_point(
    aes(fill = phase3, shape = phase3),
    size = 1.34, alpha = 0.48, stroke = 0.10, colour = "white"
  ) +
  geom_point(
    data = state_centres,
    aes(density_z, eta_z, fill = phase3, shape = phase3),
    inherit.aes = FALSE,
    size = 2.75, stroke = 0.40, colour = "#202124"
  ) +
  geom_label(
    data = state_labels,
    aes(label_x, label_y, label = phase_label),
    inherit.aes = FALSE, family = "Arial", size = FIG_GEOM_TEXT_SIZE * 0.95,
    fontface = "bold", colour = "#202124", fill = alpha("white", 0.88),
    linewidth = 0, label.padding = unit(0.10, "lines")
  ) +
  scale_fill_manual(values = phase3_cols) +
  scale_colour_manual(values = phase3_cols) +
  scale_shape_manual(values = c("pre-ictal" = 21, "ictal" = 22, "post-ictal" = 24)) +
  coord_cartesian(xlim = state_xlim, ylim = state_ylim, clip = "off") +
  labs(
    title = paste0("Propagation-medium state space"),
    x = paste0("Generalized density, ", density_sym, " (z-score)"),
    y = paste0(eta_sym, " (z-score)")
  )

layout_design <- "
AAAAAA
AAAAAA
BBCCDD
BBCCDD
"

fig <- p_a + p_b + p_c + p_d +
  plot_layout(
    design = layout_design,
    heights = c(0.72, 0.72, 1.00, 1.00),
    widths = rep(1, 6)
  ) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = FIG_PANEL_PT, face = "bold", family = "Arial", colour = "#202124"),
    plot.tag.position = c(0.010, 0.990)
  )

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 118, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  svglite::svglite(paste0(filename, ".svg"), width = w, height = h, bg = "white")
  print(plot)
  dev.off()

  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial", bg = "white")
  print(plot)
  dev.off()

  ragg::agg_tiff(
    paste0(filename, ".tiff"), width = w, height = h, units = "in",
    res = dpi, background = "white", compression = "lzw"
  )
  print(plot)
  dev.off()

  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = 300, background = "white")
  print(plot)
  dev.off()
}

save_pub_r(fig, out_base)

legend_text <- paste0(
  "Fig. 2 | Effective refractive index increases during seizures.\n",
  "a, Time-resolved ", eta_sym, " for a representative seizure, with shaded windows indicating pre-ictal, early, mid, late and post-ictal stages. ",
  "b, Stage-resolved ", eta_sym, " across 24 seizures. Thin grey lines denote seizure-level paired trajectories; small points denote individual seizure-stage means, and large outlined points with error bars denote cohort means and 95% confidence intervals. ",
  "c, Paired stage effects in ", eta_sym, " relative to the pre-ictal baseline. Small points denote seizure-level paired differences; diamonds and horizontal intervals denote cohort mean differences and 95% confidence intervals. Asterisks and P values indicate two-sided paired Wilcoxon signed-rank tests against zero with Holm correction. ",
  "d, Within-patient normalized state space of generalized density and ", eta_sym, ". Colored points denote seizure-stage observations, ellipses show the 68% normal contour for each stage, and large outlined points denote stage centroids. ",
  "Source data are from n = 24 seizures for b and c; d uses ", n_distinct(state_norm$file_id),
  " seizure summaries after excluding predefined exploratory files."
)
writeLines(legend_text, file.path(base_dir, "Fig2_redraw_legend_draft.txt"), useBytes = TRUE)

qa_lines <- c(
  "Core conclusion: the DMWA-derived effective refractive index eta increases from the pre-ictal baseline during seizure and post-ictal stages, and the increase accompanies a generalized-density/eta state-space shift.",
  "Figure archetype: asymmetric quantitative composite with a full-width representative time series and three compact cohort-level evidence panels.",
  "Backend: R only; ggplot2/patchwork/ggbeeswarm/R.matlab/svglite/cairo_pdf/ragg.",
  "Export: SVG/PDF/TIFF/PNG at double-column width.",
  "Layout revision: panels b and c are separated into absolute-stage and paired-effect views; panel d is cropped to data-supported ranges, arrows and individual trajectory lines were removed, first-version stage contours were restored, and stage labels were offset with light label backgrounds.",
  "Statistics: seizure-level summaries; paired deltas relative to pre-ictal baseline; paired Wilcoxon signed-rank tests with Holm correction.",
  paste0("Stage means Pre/Early/Mid/Late/Post: ", paste(sprintf("%.4f", stage_summary$mean), collapse = ", ")),
  paste0("Mean paired deltas Early/Mid/Late/Post: ", paste(sprintf("%.4f", paired_stats$mean_diff), collapse = ", ")),
  paste0("Holm-adjusted P Early/Mid/Late/Post: ", paste(signif(paired_stats$p_holm, 3), collapse = ", ")),
  paste0("State-space centroid density z Pre/Ictal/Post: ", paste(sprintf("%.3f", state_centres$density_z), collapse = ", ")),
  paste0("State-space centroid eta z Pre/Ictal/Post: ", paste(sprintf("%.3f", state_centres$eta_z), collapse = ", "))
)
writeLines(qa_lines, file.path(base_dir, "Fig2_redraw_QA_notes.txt"), useBytes = TRUE)

message("Fig. 2 redraw outputs written to: ", base_dir)
