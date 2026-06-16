#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(scales)
  library(grid)
  library(tibble)
})

source_dir <- file.path("data", "Fig_8")
out_dir <- file.path("figures", "Fig_8")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stage_file <- file.path(source_dir, "ALL_CASES_stage5_eta_Comparison.csv")
window_file <- file.path(source_dir, "ALL_CASES_window_level_eta_Comparison.csv")
variance_file <- file.path(source_dir, "Stats_Variance_Comparison_eta.csv")

needed <- c(stage_file, window_file, variance_file)
missing <- needed[!file.exists(needed)]
if (length(missing) > 0) {
  stop("Missing Fig. 8 source files: ", paste(missing, collapse = ", "))
}

stage_raw <- read.csv(stage_file, stringsAsFactors = FALSE, check.names = FALSE)
window_raw <- read.csv(window_file, stringsAsFactors = FALSE, check.names = FALSE)
variance_raw <- read.csv(variance_file, stringsAsFactors = FALSE, check.names = FALSE)

required_stage <- c("case_id", "Model", "phase5", "n_windows", "eta_mean", "eta_std", "eta_median")
required_window <- c("case_id", "window_idx", "phase5", "Model", "eta", "num_hyperedges")
required_variance <- c("Model", "Phase", "Variance")
if (!all(required_stage %in% names(stage_raw))) {
  stop("Stage table missing columns: ", paste(setdiff(required_stage, names(stage_raw)), collapse = ", "))
}
if (!all(required_window %in% names(window_raw))) {
  stop("Window table missing columns: ", paste(setdiff(required_window, names(window_raw)), collapse = ", "))
}
if (!all(required_variance %in% names(variance_raw))) {
  stop("Variance table missing columns: ", paste(setdiff(required_variance, names(variance_raw)), collapse = ", "))
}

model_levels <- c("DMW-HLG", "High-Order(3-nodes)", "Low-Order(2-nodes)")
model_labels <- c(
  "DMW-HLG" = "DMW-HLG",
  "High-Order(3-nodes)" = "High-order",
  "Low-Order(2-nodes)" = "Low-order"
)
phase_levels <- c("pre-ictal", "early", "mid", "late", "post-ictal")
phase_labels <- c(
  "pre-ictal" = "Pre",
  "early" = "Early",
  "mid" = "Mid",
  "late" = "Late",
  "post-ictal" = "Post"
)

stage <- stage_raw %>%
  mutate(
    Model = factor(Model, levels = model_levels),
    phase5 = factor(phase5, levels = phase_levels),
    eta_mean = as.numeric(eta_mean),
    eta_std = as.numeric(eta_std),
    eta_median = as.numeric(eta_median),
    n_windows = as.integer(n_windows)
  ) %>%
  filter(!is.na(Model), !is.na(phase5), !is.na(eta_mean))

window_data <- window_raw %>%
  mutate(
    Model = factor(Model, levels = model_levels),
    phase5 = factor(phase5, levels = phase_levels),
    eta = as.numeric(eta),
    num_hyperedges = as.numeric(num_hyperedges)
  ) %>%
  filter(!is.na(Model), !is.na(phase5), !is.na(eta))

variance_tbl <- variance_raw %>%
  mutate(
    Model = factor(Model, levels = model_levels),
    Phase = factor(Phase, levels = phase_levels),
    Variance = as.numeric(Variance)
  ) %>%
  filter(!is.na(Model), !is.na(Phase), !is.na(Variance))

# Recalculate source-data products from audited inputs.
stage_summary <- stage %>%
  group_by(Model, phase5) %>%
  summarise(
    n_seizures = n_distinct(case_id),
    mean_eta = mean(eta_mean, na.rm = TRUE),
    sd_eta = sd(eta_mean, na.rm = TRUE),
    se_eta = sd_eta / sqrt(n_seizures),
    ci95_eta = qt(0.975, df = n_seizures - 1) * se_eta,
    median_eta = median(eta_mean, na.rm = TRUE),
    q1_eta = quantile(eta_mean, 0.25, na.rm = TRUE, names = FALSE),
    q3_eta = quantile(eta_mean, 0.75, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  )

paired_stats <- stage %>%
  select(case_id, Model, phase5, eta_mean) %>%
  pivot_wider(names_from = phase5, values_from = eta_mean) %>%
  filter(!is.na(`pre-ictal`), !is.na(early)) %>%
  group_by(Model) %>%
  summarise(
    n_seizures = n(),
    mean_pre = mean(`pre-ictal`, na.rm = TRUE),
    mean_early = mean(early, na.rm = TRUE),
    delta_early_minus_pre = mean(early - `pre-ictal`, na.rm = TRUE),
    sd_delta = sd(early - `pre-ictal`, na.rm = TRUE),
    paired_effect_dz = delta_early_minus_pre / sd_delta,
    p_value = wilcox.test(early, `pre-ictal`, paired = TRUE, exact = FALSE, correct = FALSE)$p.value,
    .groups = "drop"
  )

calc_roc <- function(scores, labels, positive = "early") {
  valid <- !is.na(scores) & !is.na(labels)
  scores <- scores[valid]
  labels <- labels[valid]
  y <- labels == positive
  pos <- sum(y)
  neg <- sum(!y)
  if (pos == 0 || neg == 0) {
    stop("ROC requires both positive and negative observations.")
  }

  # Rank-sum AUC with average ranks for ties.
  auc <- (sum(rank(scores, ties.method = "average")[y]) - pos * (pos + 1) / 2) / (pos * neg)
  flip_scores <- FALSE
  if (auc < 0.5) {
    scores <- -scores
    auc <- (sum(rank(scores, ties.method = "average")[y]) - pos * (pos + 1) / 2) / (pos * neg)
    flip_scores <- TRUE
  }

  ord <- order(scores, decreasing = TRUE)
  sorted_y <- y[ord]
  tpr <- c(0, cumsum(sorted_y) / pos, 1)
  fpr <- c(0, cumsum(!sorted_y) / neg, 1)
  tibble(fpr = fpr, tpr = tpr, auc = auc, n_positive = pos, n_negative = neg, flipped = flip_scores)
}

roc_coords <- window_data %>%
  filter(phase5 %in% c("pre-ictal", "early")) %>%
  group_by(Model) %>%
  group_modify(~ calc_roc(.x$eta, as.character(.x$phase5), positive = "early")) %>%
  ungroup()

roc_summary <- roc_coords %>%
  group_by(Model) %>%
  summarise(
    auc = first(auc),
    n_positive_early_windows = first(n_positive),
    n_negative_pre_windows = first(n_negative),
    flipped = first(flipped),
    .groups = "drop"
  )

variance_calc <- stage %>%
  group_by(Model, phase5) %>%
  summarise(Variance_recalc = var(eta_mean, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    variance_tbl %>% rename(phase5 = Phase),
    by = c("Model", "phase5")
  ) %>%
  mutate(abs_diff = abs(Variance_recalc - Variance))

expected_stage_rows <- length(model_levels) * length(phase_levels) * 24
expected_window_rows <- 10276 * length(model_levels)
preplot_checks <- tibble(
  check = c(
    "stage file exists",
    "window file exists",
    "variance file exists",
    "stage rows equal 24 seizures x 5 phases x 3 models",
    "window rows equal 10276 windows x 3 models",
    "stage model-phase cells have n=24",
    "variance table matches recalculated seizure-level variance",
    "ROC uses pre-ictal and early window-level observations"
  ),
  value = c(
    as.character(file.exists(stage_file)),
    as.character(file.exists(window_file)),
    as.character(file.exists(variance_file)),
    as.character(nrow(stage)),
    as.character(nrow(window_data)),
    paste(range(stage_summary$n_seizures), collapse = "-"),
    formatC(max(variance_calc$abs_diff, na.rm = TRUE), format = "e", digits = 3),
    paste0(
      "pre=", unique(roc_summary$n_negative_pre_windows),
      "; early=", unique(roc_summary$n_positive_early_windows)
    )
  ),
  status = c(
    ifelse(file.exists(stage_file), "PASS", "FAIL"),
    ifelse(file.exists(window_file), "PASS", "FAIL"),
    ifelse(file.exists(variance_file), "PASS", "FAIL"),
    ifelse(nrow(stage) == expected_stage_rows, "PASS", "FAIL"),
    ifelse(nrow(window_data) == expected_window_rows, "PASS", "FAIL"),
    ifelse(all(stage_summary$n_seizures == 24), "PASS", "FAIL"),
    ifelse(max(variance_calc$abs_diff, na.rm = TRUE) < 1e-15, "PASS", "FAIL"),
    ifelse(all(roc_summary$n_negative_pre_windows == 4139) &&
             all(roc_summary$n_positive_early_windows == 964), "PASS", "FAIL")
  )
)

write.csv(stage, file.path(out_dir, "Fig8_stage_source_data.csv"), row.names = FALSE)
write.csv(window_data, file.path(out_dir, "Fig8_window_source_data.csv"), row.names = FALSE)
write.csv(stage_summary, file.path(out_dir, "Fig8_stage_summary.csv"), row.names = FALSE)
write.csv(paired_stats, file.path(out_dir, "Fig8_paired_stage_stats.csv"), row.names = FALSE)
write.csv(variance_calc, file.path(out_dir, "Fig8_variance_source_data.csv"), row.names = FALSE)
write.csv(roc_coords, file.path(out_dir, "Fig8_roc_coordinates.csv"), row.names = FALSE)
write.csv(roc_summary, file.path(out_dir, "Fig8_roc_auc_summary.csv"), row.names = FALSE)
write.csv(preplot_checks, file.path(out_dir, "Fig8_preplot_consistency_checks.csv"), row.names = FALSE)

if (any(preplot_checks$status != "PASS")) {
  stop("Pre-plot consistency checks failed; inspect Fig8_preplot_consistency_checks.csv")
}

# Same restrained method-family palette used in the recent Fig. 7 revision:
# the proposed DMW-HLG is the coral accent; lower-order baselines use blue/grey.
model_palette <- c(
  "DMW-HLG" = "#D96A63",
  "High-Order(3-nodes)" = "#5A84AF",
  "Low-Order(2-nodes)" = "#8B919A"
)

theme_nature <- function(base_size = 6.5, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(colour = "#202124"),
      axis.line = element_line(linewidth = 0.34, colour = "#232323"),
      axis.ticks = element_line(linewidth = 0.34, colour = "#232323"),
      axis.title = element_text(size = base_size + 0.1),
      axis.text = element_text(size = base_size - 0.2, colour = "#333333"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.2),
      plot.title = element_text(size = base_size + 0.8, face = "bold", hjust = 0, colour = "#202124"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 0.2, face = "bold", colour = "#222222"),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(4, 5, 4, 5)
    )
}

paired_label_df <- paired_stats %>%
  mutate(
    label = paste0(
      "\u0394=", formatC(delta_early_minus_pre, format = "f", digits = 4),
      "\nP=", formatC(p_value, format = "f", digits = 3)
    ),
    x = 1.05,
    y = 1.0125
  )

p_a <- ggplot(stage, aes(x = phase5, y = eta_mean, colour = Model)) +
  geom_boxplot(
    width = 0.50,
    outlier.shape = NA,
    linewidth = 0.30,
    fill = "white",
    alpha = 0.88
  ) +
  geom_point(
    position = position_jitter(width = 0.075, height = 0, seed = 8),
    size = 0.95,
    alpha = 0.48,
    stroke = 0
  ) +
  geom_line(
    data = stage_summary,
    aes(x = phase5, y = mean_eta, group = Model, colour = Model),
    linewidth = 0.44,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = stage_summary,
    aes(x = phase5, y = mean_eta, colour = Model),
    size = 1.35,
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = stage_summary,
    aes(x = phase5, ymin = mean_eta - ci95_eta, ymax = mean_eta + ci95_eta, colour = Model),
    width = 0.13,
    linewidth = 0.30,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = paired_label_df,
    aes(x = x, y = y, label = label, colour = Model),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 1.8,
    lineheight = 0.92,
    family = "Arial"
  ) +
  facet_wrap(~ Model, nrow = 1, labeller = as_labeller(model_labels)) +
  scale_colour_manual(values = model_palette, breaks = model_levels, labels = model_labels) +
  scale_x_discrete(labels = phase_labels) +
  coord_cartesian(ylim = c(0.945, 1.014), clip = "off") +
  labs(
    x = NULL,
    y = expression(eta),
    title = paste0("Stage-resolved ", intToUtf8(0x03b7), " across order representations")
  ) +
  guides(colour = "none") +
  theme_nature() +
  theme(
    panel.spacing.x = unit(4.0, "mm"),
    axis.text.x = element_text(size = 6.0),
    plot.margin = margin(5, 4, 5, 5)
  )

variance_plot_data <- variance_tbl %>%
  mutate(ModelLabel = factor(model_labels[as.character(Model)], levels = unname(model_labels)))

p_b <- ggplot(variance_plot_data, aes(x = Phase, y = Variance, colour = Model, group = Model)) +
  geom_line(linewidth = 0.50) +
  geom_point(size = 1.45) +
  scale_colour_manual(values = model_palette, breaks = model_levels, labels = model_labels) +
  scale_x_discrete(labels = phase_labels) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.00002),
    expand = expansion(mult = c(0.03, 0.12))
  ) +
  labs(
    x = NULL,
    y = expression("Cross-seizure variance of " * eta),
    title = "Cross-seizure variability"
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = 0.50, size = 1.25))) +
  theme_nature() +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.background = element_blank(),
    legend.box.margin = margin(-3, 0, 0, 0),
    legend.spacing.x = unit(2.0, "mm"),
    axis.text.x = element_text(size = 6.0)
  )

auc_label_df <- roc_summary %>%
  arrange(Model) %>%
  mutate(
    label = paste0(model_labels[as.character(Model)], " AUC=", formatC(auc, format = "f", digits = 3)),
    x = 0.97,
    y = c(0.24, 0.16, 0.08)
  )

p_c <- ggplot(roc_coords, aes(x = fpr, y = tpr, colour = Model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.34, colour = "#9A9A9A") +
  geom_step(linewidth = 0.58, direction = "hv") +
  geom_text(
    data = auc_label_df,
    aes(x = x, y = y, label = label, colour = Model),
    inherit.aes = FALSE,
    hjust = 1,
    size = 1.85,
    family = "Arial"
  ) +
  scale_colour_manual(values = model_palette, breaks = model_levels, labels = model_labels) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  labs(
    x = "False-positive rate",
    y = "True-positive rate",
    title = "Pre-ictal vs early discrimination"
  ) +
  guides(colour = "none") +
  theme_nature() +
  theme(
    axis.text.x = element_text(size = 6.0),
    axis.text.y = element_text(size = 6.0),
    plot.margin = margin(4, 5, 4, 5)
  )

fig <- p_a / (p_b | p_c) +
  plot_layout(heights = c(1.42, 1), widths = c(1.08, 0.92)) +
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

base_file <- file.path(out_dir, "Fig8_eta_order_representation_comparison_method_first_label_clear_titlebold")
save_pub(fig, base_file)

export_files <- paste0(base_file, c(".svg", ".pdf", ".png", ".tiff"))
post_checks <- tibble(
  file = basename(export_files),
  exists = file.exists(export_files),
  size_bytes = ifelse(file.exists(export_files), file.info(export_files)$size, NA_real_),
  status = ifelse(file.exists(export_files) & file.info(export_files)$size > 10000, "PASS", "FAIL")
)
write.csv(post_checks, file.path(out_dir, "Fig8_postplot_export_checks_method_first_label_clear.csv"), row.names = FALSE)

figure_numeric_summary <- paired_stats %>%
  left_join(roc_summary, by = "Model") %>%
  left_join(
    variance_tbl %>%
      group_by(Model) %>%
      summarise(mean_cross_seizure_variance = mean(Variance, na.rm = TRUE), .groups = "drop"),
    by = "Model"
  ) %>%
  mutate(ModelLabel = model_labels[as.character(Model)]) %>%
  select(
    Model, ModelLabel, n_seizures, mean_pre, mean_early, delta_early_minus_pre,
    p_value, paired_effect_dz, auc, n_negative_pre_windows, n_positive_early_windows,
    mean_cross_seizure_variance
  )
write.csv(figure_numeric_summary, file.path(out_dir, "Fig8_figure_numeric_summary_method_first_label_clear.csv"), row.names = FALSE)

writeLines(
  c(
    "Fig. 8 figure contract and audit note",
    "",
    "Core conclusion: higher-order representations provide stronger pre-ictal-to-early eta sensitivity than the low-order baseline, with DMW-HLG showing the largest mean early increase and AUC but also the greatest cross-seizure variance.",
    "Archetype: quantitative grid with a dominant stage-resolved comparison panel and two supporting stability/discrimination panels.",
    "Backend: R only (ggplot2 + patchwork + svglite/cairo_pdf/ragg).",
    "Display order: DMW-HLG first, followed by the high-order and low-order baselines, matching the Fig. 7 method-first visual logic.",
    "Panel-a statistics labels are placed in the upper-left whitespace of each facet to avoid overlap with boxplots and whiskers.",
    "Palette: DMW-HLG #D96A63, High-order #5A84AF, Low-order #8B919A, matching the restrained Fig. 7 method-highlight preference.",
    "Data: source tables copied from audited Fig. 8 directory; no values were fabricated or manually edited.",
    "Statistics: stage-level paired pre-ictal-to-early comparisons use two-sided Wilcoxon signed-rank tests on n=24 seizure-level stage means; ROC curves use window-level eta values for pre-ictal and early stages."
  ),
  file.path(out_dir, "Fig8_figure_contract_and_audit_note_method_first_label_clear.txt"),
  useBytes = TRUE
)

if (any(post_checks$status == "FAIL")) {
  stop("Export checks failed; inspect Fig8_postplot_export_checks_method_first_label_clear.csv")
}

message("Fig. 8 complete: ", base_file)
