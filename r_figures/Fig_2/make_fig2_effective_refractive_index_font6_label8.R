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
})

FIG_TEXT_PT <- 6
FIG_PANEL_PT <- 8
FIG_GEOM_TEXT_SIZE <- FIG_TEXT_PT / 2.845276


base_candidates <- Sys.glob("F:/6/*0514-/nature_fig/Fig_2")
if (length(base_candidates) < 1) {
  stop("Cannot locate Fig_2 directory under F:/6/*0514-/nature_fig/Fig_2")
}
base_dir <- normalizePath(base_candidates[1], winslash = "/", mustWork = TRUE)
out_base <- file.path(base_dir, "Fig2_effective_refractive_index_increases_with_state_space_font6_label8")

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

theme_nature <- function(base_size = FIG_TEXT_PT) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#202124"),
      axis.ticks = element_line(linewidth = 0.30, colour = "#202124"),
      axis.text = element_text(size = FIG_TEXT_PT, colour = "#202124"),
      axis.title = element_text(size = FIG_TEXT_PT, colour = "#202124"),
      plot.title = element_text(size = FIG_TEXT_PT, face = "bold", hjust = 0),
      legend.position = "none",
      panel.grid.major = element_line(linewidth = 0.18, colour = "#ECEFF3"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(3.0, 3.0, 3.0, 3.0)
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

eta_window <- read_csv(file.path(base_dir, "lhs_cut07_window_level_eta.csv"), show_col_types = FALSE) %>%
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
    phase_label = factor(effect_labels[as.integer(phase5)], levels = effect_labels)
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
    stars = stars_from_p(p_holm)
  )

write_csv(stage_summary, file.path(base_dir, "Fig2_relayout_stage_summary.csv"))
write_csv(paired_stats, file.path(base_dir, "Fig2_relayout_paired_statistics.csv"))

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
  geom_vline(xintercept = transition_x, linetype = "22", linewidth = 0.30, colour = "#8A8F98") +
  geom_text(
    data = stage_bounds,
    aes(x = xmid, y = eta_ymax + eta_pad * 0.34, label = phase_label),
    inherit.aes = FALSE, size = FIG_GEOM_TEXT_SIZE, family = "Arial", fontface = "bold", colour = "#202124"
  ) +
  scale_fill_manual(values = phase_cols) +
  coord_cartesian(ylim = c(eta_ymin - eta_pad, eta_ymax + eta_pad * 0.70), clip = "off") +
  labs(
    title = paste0("Representative seizure: time-resolved ", eta_sym),
    x = "Time windows",
    y = paste0("Effective refractive index, ", eta_sym)
  )

p_b <- ggplot(eta_stage, aes(phase_label, eta_mean, group = case_id)) +
  geom_line(colour = "#B9C1CC", linewidth = 0.22, alpha = 0.18) +
  geom_point(aes(fill = phase5), shape = 21, size = 0.95, stroke = 0.08, colour = "white", alpha = 0.55) +
  geom_line(
    data = stage_summary,
    aes(phase_label, mean, group = 1),
    inherit.aes = FALSE, colour = "#202124", linewidth = 0.48
  ) +
  geom_errorbar(
    data = stage_summary,
    aes(phase_label, ymin = mean - ci, ymax = mean + ci),
    inherit.aes = FALSE, width = 0.12, linewidth = 0.32, colour = "#202124"
  ) +
  geom_point(
    data = stage_summary,
    aes(phase_label, mean, fill = phase5),
    inherit.aes = FALSE, shape = 21, size = 2.45, stroke = 0.36, colour = "#202124"
  ) +
  scale_fill_manual(values = phase_cols) +
  labs(
    title = paste0("Cohort trajectory of ", eta_sym),
    x = NULL,
    y = paste0("Stage mean ", eta_sym)
  )

star_y <- max(delta_long$delta_eta, na.rm = TRUE) + 0.006

p_c <- ggplot(delta_long, aes(phase_label, delta_eta, fill = phase5)) +
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.28, colour = "#7A828C") +
  ggbeeswarm::geom_quasirandom(
    width = 0.12, shape = 21, size = 1.25, stroke = 0.10,
    colour = "white", alpha = 0.58
  ) +
  geom_errorbar(
    data = paired_stats,
    aes(phase_label, ymin = ci_low, ymax = ci_high),
    inherit.aes = FALSE, width = 0.12, linewidth = 0.34, colour = "#202124"
  ) +
  geom_line(
    data = paired_stats,
    aes(phase_label, mean_diff, group = 1),
    inherit.aes = FALSE, linewidth = 0.45, colour = "#202124"
  ) +
  geom_point(
    data = paired_stats,
    aes(phase_label, mean_diff, fill = phase5),
    inherit.aes = FALSE, shape = 23, size = 2.35, stroke = 0.35, colour = "#202124"
  ) +
  geom_text(
    data = paired_stats,
    aes(phase_label, y = star_y, label = stars),
    inherit.aes = FALSE, family = "Arial", size = FIG_GEOM_TEXT_SIZE, fontface = "bold", colour = "#202124"
  ) +
  scale_fill_manual(values = phase_cols[effect_levels]) +
  coord_cartesian(
    ylim = c(min(delta_long$delta_eta, na.rm = TRUE) - 0.004, star_y + 0.004),
    clip = "off"
  ) +
  labs(
    title = paste0("Paired stage effect relative to Pre"),
    x = NULL,
    y = paste0(delta_sym, eta_sym, " relative to Pre")
  )

mat_dir <- file.path(base_dir, "seizure_data")
mat_files <- list.files(mat_dir, pattern = "\\.mat$", full.names = TRUE)
exclude_files <- c("gwh_s2", "ssh_s3", "ssh_s5", "ssh_s6")

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
  select(-density_sd, -eta_sd)

state_centres <- state_norm %>%
  group_by(phase3) %>%
  summarise(
    density_z = mean(density_z, na.rm = TRUE),
    eta_z = mean(eta_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(phase3 = factor(phase3, levels = phase3_levels))

arrow_df <- state_centres %>%
  arrange(phase3) %>%
  mutate(
    xend = lead(density_z),
    yend = lead(eta_z)
  ) %>%
  filter(!is.na(xend))

write_csv(state_norm, file.path(base_dir, "Fig2_with_state_space_points.csv"))

p_d <- ggplot(state_norm, aes(density_z, eta_z)) +
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.25, colour = "#A8AFB8") +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.25, colour = "#A8AFB8") +
  geom_point(
    aes(fill = phase3, shape = phase3),
    size = 1.45, alpha = 0.44, stroke = 0.12, colour = "white"
  ) +
  geom_segment(
    data = arrow_df,
    aes(x = density_z, y = eta_z, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(1.55, "mm"), type = "closed"),
    linewidth = 0.42, colour = "#202124", alpha = 0.75
  ) +
  geom_point(
    data = state_centres,
    aes(density_z, eta_z, fill = phase3, shape = phase3),
    inherit.aes = FALSE,
    size = 3.0, stroke = 0.42, colour = "#202124"
  ) +
  scale_fill_manual(values = phase3_cols, labels = phase3_labels, name = NULL) +
  scale_shape_manual(values = c("pre-ictal" = 21, "ictal" = 22, "post-ictal" = 24), labels = phase3_labels, name = NULL) +
  coord_fixed(xlim = c(-2.8, 2.8), ylim = c(-2.8, 2.8), clip = "off") +
  labs(
    title = paste0("State space of ", eta_sym, " and ", "\u03b4", "\u00b2"),
    x = "Density delta(2) (z-score)",
    y = paste0(eta_sym, " (z-score)")
  ) +
  theme(
    legend.position = c(0.04, 0.96),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.box.background = element_rect(fill = scales::alpha("white", 0.86), colour = NA),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = FIG_TEXT_PT),
    legend.key.height = unit(2.8, "mm"),
    legend.key.width = unit(3.4, "mm")
  )

layout_design <- "
AAAAAA
AAAAAA
BBBDDD
BBBDDD
CCCDDD
CCCDDD
"

fig <- p_a + p_b + p_c + p_d +
  plot_layout(
    design = layout_design,
    heights = c(0.72, 0.72, 0.58, 0.58, 0.58, 0.58)
  ) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = FIG_PANEL_PT, face = "bold", family = "Arial", colour = "#202124"),
    plot.tag.position = c(0.012, 0.988)
  )

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 158, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  svglite::svglite(paste0(filename, ".svg"), width = w, height = h, bg = "white")
  print(plot)
  dev.off()

  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial", bg = "white")
  print(plot)
  dev.off()

  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi, background = "white", compression = "lzw")
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
  "b, Stage-resolved ", eta_sym, " across 24 seizures. Thin grey lines denote seizure-level paired trajectories; small points denote individual seizure-stage means, and large outlined points denote cohort means with 95% confidence intervals. ",
  "c, Paired change in ", eta_sym, " relative to the pre-ictal baseline. Points denote seizure-level stage effects; diamonds and error bars denote cohort mean differences and 95% confidence intervals. Asterisks indicate paired Wilcoxon signed-rank tests against the pre-ictal baseline with Holm correction. Source data are from n = 24 seizures."
)
writeLines(legend_text, file.path(base_dir, "Fig2_relayout_legend_draft.txt"), useBytes = TRUE)

legend_text_with_state <- paste0(
  "Fig. 2 | Effective refractive index increases during seizures.\n",
  "a, Time-resolved ", eta_sym, " for a representative seizure, with shaded windows indicating pre-ictal, early, mid, late and post-ictal stages. ",
  "b, Stage-resolved ", eta_sym, " across 24 seizures. Thin grey lines denote seizure-level paired trajectories; small points denote individual seizure-stage means, and large outlined points denote cohort means with 95% confidence intervals. ",
  "c, Paired change in ", eta_sym, " relative to the pre-ictal baseline. Points denote seizure-level stage effects; diamonds and error bars denote cohort mean differences and 95% confidence intervals. Asterisks indicate paired Wilcoxon signed-rank tests against the pre-ictal baseline with Holm correction. ",
  "d, State-space representation of generalized density and ", eta_sym, " after within-patient z-scoring. Small points denote seizure-stage observations; large outlined points denote stage centroids, and arrows indicate the pre-ictal to ictal to post-ictal progression. Source data are from n = 24 seizures for b and c; d uses ",
  n_distinct(state_norm$file_id), " seizure summaries after excluding predefined exploratory files."
)
writeLines(legend_text_with_state, file.path(base_dir, "Fig2_with_state_space_legend_draft.txt"), useBytes = TRUE)

qa_lines <- c(
  "Core conclusion: the DMWA-derived effective refractive index eta increases from the pre-ictal baseline during seizure and post-ictal stages.",
  "Figure archetype: asymmetric quantitative layout with a full-width representative time series, two compact cohort-level effect panels and a larger state-space context panel.",
  "Backend: R only; ggplot2/patchwork/ggbeeswarm/R.matlab/svglite/cairo_pdf/ragg.",
  "Export: SVG/PDF/TIFF/PNG at double-column width.",
  "Statistics: seizure-level summaries; paired deltas relative to pre-ictal baseline; Wilcoxon/Holm statistics exported.",
  paste0("Stage means Pre/Early/Mid/Late/Post: ", paste(sprintf("%.4f", stage_summary$mean), collapse = ", ")),
  paste0("Mean paired deltas Early/Mid/Late/Post: ", paste(sprintf("%.4f", paired_stats$mean_diff), collapse = ", ")),
  paste0("Holm-adjusted P Early/Mid/Late/Post: ", paste(signif(paired_stats$p_holm, 3), collapse = ", "))
)
writeLines(qa_lines, file.path(base_dir, "Fig2_relayout_QA_notes.txt"), useBytes = TRUE)

writeLines(qa_lines, file.path(base_dir, "Fig2_with_state_space_QA_notes.txt"), useBytes = TRUE)

message("Fig. 2 relayout outputs written to: ", base_dir)
