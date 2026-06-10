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

data_dir <- "E:/wcldematlab/keep/new_idea/8 - n_u_v/2_Fig7"
out_dir <- file.path("F:/6", "\u521d\u7a3f0514-", "nature_fig", "Fig_7")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stage_file <- file.path(data_dir, "ALL_CASES_stage5_eta_ablation.csv")
window_file <- file.path(data_dir, "ALL_CASES_window_eta_ablation.csv")
variance_file <- file.path(data_dir, "variance_table_eta_ablation.csv")
matlab_ablation_script <- file.path(data_dir, "compare_ablation_models_eta.m")
matlab_roc_script <- file.path(data_dir, "plot_ROC_7_eta.m")
gain_script <- "E:/wcldematlab/keep/new_idea/8 - n_u_v/2/gain_hyperEdges_23.m"

required_files <- c(stage_file, window_file, variance_file, matlab_ablation_script, matlab_roc_script, gain_script)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required Fig. 7 source files:\n", paste(missing_files, collapse = "\n"))
}

phase_levels <- c("pre-ictal", "early", "mid", "late", "post-ictal")
phase_labels <- c("Pre", "Early", "Mid", "Late", "Post")
model_cols <- c("eta_orig_mean", "eta_cov_mean", "eta_deg_mean")
model_labels <- c(
  eta_orig_mean = "Original",
  eta_cov_mean = "Coverage only",
  eta_deg_mean = "DegreeBias only"
)
model_levels <- unname(model_labels)
model_palette <- c(
  "Original" = "#232323",
  "Coverage only" = "#477DA0",
  "DegreeBias only" = "#8B6F4E"
)

stage_raw <- read.csv(stage_file, check.names = FALSE, stringsAsFactors = FALSE)
window_raw <- read.csv(window_file, check.names = FALSE, stringsAsFactors = FALSE)
variance_ref <- read.csv(variance_file, check.names = FALSE, stringsAsFactors = FALSE)

stage_long <- stage_raw %>%
  pivot_longer(all_of(model_cols), names_to = "model_col", values_to = "eta") %>%
  mutate(
    model = factor(unname(model_labels[model_col]), levels = model_levels),
    phase5 = factor(phase5, levels = phase_levels),
    phase_label = factor(phase_labels[match(as.character(phase5), phase_levels)], levels = phase_labels)
  ) %>%
  arrange(model, case_id, phase5)

write.csv(stage_long %>% select(case_id, phase5, phase_label, n_windows, model, eta),
          file.path(out_dir, "Fig7_stage_source_data.csv"), row.names = FALSE)

stage_summary <- stage_long %>%
  group_by(model, phase5, phase_label) %>%
  summarise(
    n_seizures = n(),
    mean_eta = mean(eta, na.rm = TRUE),
    sd_eta = sd(eta, na.rm = TRUE),
    sem_eta = sd_eta / sqrt(n_seizures),
    ci95_eta = qt(0.975, df = n_seizures - 1) * sem_eta,
    median_eta = median(eta, na.rm = TRUE),
    q1_eta = quantile(eta, 0.25, na.rm = TRUE),
    q3_eta = quantile(eta, 0.75, na.rm = TRUE),
    min_eta = min(eta, na.rm = TRUE),
    max_eta = max(eta, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(stage_summary, file.path(out_dir, "Fig7_stage_summary.csv"), row.names = FALSE)

paired_stat <- function(df, model_col, a, b) {
  x <- df %>%
    filter(phase5 == a) %>%
    select(case_id, x = all_of(model_col))
  y <- df %>%
    filter(phase5 == b) %>%
    select(case_id, y = all_of(model_col))
  merged <- inner_join(x, y, by = "case_id")
  delta <- merged$y - merged$x
  p_val <- suppressWarnings(wilcox.test(merged$x, merged$y, paired = TRUE, exact = FALSE)$p.value)
  tibble(
    model = unname(model_labels[[model_col]]),
    comparison = paste(a, "vs", b),
    n_seizures = nrow(merged),
    mean_delta_eta = mean(delta, na.rm = TRUE),
    median_delta_eta = median(delta, na.rm = TRUE),
    sd_delta_eta = sd(delta, na.rm = TRUE),
    paired_dz = ifelse(sd(delta, na.rm = TRUE) > 0, mean(delta, na.rm = TRUE) / sd(delta, na.rm = TRUE), NA_real_),
    p_value = p_val
  )
}

comparison_pairs <- list(
  c("pre-ictal", "early"),
  c("early", "mid"),
  c("mid", "late"),
  c("late", "post-ictal"),
  c("pre-ictal", "mid"),
  c("mid", "post-ictal")
)
paired_stats <- bind_rows(lapply(model_cols, function(model_col) {
  bind_rows(lapply(comparison_pairs, function(pair) paired_stat(stage_raw, model_col, pair[1], pair[2])))
})) %>%
  mutate(
    model = factor(model, levels = model_levels),
    p_label = case_when(
      is.na(p_value) ~ "P=NA",
      p_value < 0.001 ~ "P<0.001",
      TRUE ~ paste0("P=", formatC(p_value, format = "f", digits = 3))
    )
  )
write.csv(paired_stats, file.path(out_dir, "Fig7_paired_stage_stats.csv"), row.names = FALSE)

variance_calc <- stage_long %>%
  group_by(model, phase5, phase_label) %>%
  summarise(var_eta = var(eta, na.rm = TRUE), .groups = "drop") %>%
  mutate(model = factor(model, levels = model_levels))
write.csv(variance_calc, file.path(out_dir, "Fig7_variance_source_data.csv"), row.names = FALSE)

variance_calc_wide <- variance_calc %>%
  mutate(
    Phase = as.character(phase5),
    var_col = recode(as.character(model),
                     "Original" = "Var_Orig",
                     "Coverage only" = "Var_Cov",
                     "DegreeBias only" = "Var_Deg")
  ) %>%
  select(Phase, var_col, var_eta) %>%
  pivot_wider(names_from = var_col, values_from = var_eta) %>%
  arrange(match(Phase, phase_levels))

variance_joined <- variance_calc_wide %>%
  inner_join(variance_ref, by = "Phase", suffix = c("_calc", "_ref")) %>%
  mutate(
    diff_orig = abs(Var_Orig_calc - Var_Orig_ref),
    diff_cov = abs(Var_Cov_calc - Var_Cov_ref),
    diff_deg = abs(Var_Deg_calc - Var_Deg_ref)
  )

roc_auc <- function(labels, scores, positive = "early") {
  ok <- is.finite(scores) & !is.na(labels)
  labels <- labels[ok]
  scores <- scores[ok]
  pos <- scores[labels == positive]
  neg <- scores[labels != positive]
  if (length(pos) == 0 || length(neg) == 0) {
    return(NA_real_)
  }
  comp <- outer(pos, neg, "-")
  (sum(comp > 0) + 0.5 * sum(comp == 0)) / (length(pos) * length(neg))
}

roc_curve <- function(labels, scores, positive = "early") {
  ok <- is.finite(scores) & !is.na(labels)
  labels <- labels[ok]
  scores <- scores[ok]
  thresholds <- c(Inf, sort(unique(scores), decreasing = TRUE), -Inf)
  bind_rows(lapply(thresholds, function(thr) {
    pred <- scores >= thr
    tp <- sum(pred & labels == positive)
    fp <- sum(pred & labels != positive)
    fn <- sum(!pred & labels == positive)
    tn <- sum(!pred & labels != positive)
    tibble(
      threshold = thr,
      fpr = ifelse(fp + tn > 0, fp / (fp + tn), NA_real_),
      tpr = ifelse(tp + fn > 0, tp / (tp + fn), NA_real_)
    )
  })) %>%
    distinct(fpr, tpr, .keep_all = TRUE) %>%
    arrange(fpr, tpr)
}

roc_input <- stage_raw %>%
  filter(phase5 %in% c("pre-ictal", "early"))

roc_stats <- bind_rows(lapply(model_cols, function(model_col) {
  labels <- roc_input$phase5
  scores <- roc_input[[model_col]]
  auc <- roc_auc(labels, scores, positive = "early")
  flipped <- FALSE
  if (!is.na(auc) && auc < 0.5) {
    scores <- -scores
    auc <- roc_auc(labels, scores, positive = "early")
    flipped <- TRUE
  }
  tibble(
    model = factor(unname(model_labels[[model_col]]), levels = model_levels),
    n_pre = sum(labels == "pre-ictal"),
    n_early = sum(labels == "early"),
    auc = auc,
    score_flipped = flipped
  )
}))
write.csv(roc_stats, file.path(out_dir, "Fig7_roc_auc_summary.csv"), row.names = FALSE)

roc_coords <- bind_rows(lapply(model_cols, function(model_col) {
  labels <- roc_input$phase5
  scores <- roc_input[[model_col]]
  auc <- roc_auc(labels, scores, positive = "early")
  if (!is.na(auc) && auc < 0.5) {
    scores <- -scores
  }
  roc_curve(labels, scores, positive = "early") %>%
    mutate(model = factor(unname(model_labels[[model_col]]), levels = model_levels))
}))
write.csv(roc_coords, file.path(out_dir, "Fig7_roc_coordinates.csv"), row.names = FALSE)

window_means <- window_raw %>%
  filter(phase5 %in% phase_levels) %>%
  group_by(case_id, phase5) %>%
  summarise(
    n_windows_from_window_csv = n(),
    eta_orig_mean_from_window = mean(eta_orig, na.rm = TRUE),
    eta_cov_mean_from_window = mean(eta_cov, na.rm = TRUE),
    eta_deg_mean_from_window = mean(eta_deg, na.rm = TRUE),
    .groups = "drop"
  )
stage_joined <- stage_raw %>%
  left_join(window_means, by = c("case_id", "phase5")) %>%
  mutate(
    n_windows_diff = n_windows - n_windows_from_window_csv,
    diff_orig = abs(eta_orig_mean - eta_orig_mean_from_window),
    diff_cov = abs(eta_cov_mean - eta_cov_mean_from_window),
    diff_deg = abs(eta_deg_mean - eta_deg_mean_from_window)
  )
write.csv(stage_joined, file.path(out_dir, "Fig7_window_to_stage_aggregation_check.csv"), row.names = FALSE)

read_text_safe <- function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}
ablation_code <- read_text_safe(matlab_ablation_script)
gain_code <- read_text_safe(gain_script)
roc_code <- read_text_safe(matlab_roc_script)

pre_checks <- tibble(
  check = c(
    "required source files present",
    "stage rows",
    "unique seizures",
    "phase coverage per seizure",
    "window rows match summed stage windows",
    "window-to-stage eta aggregation",
    "variance table matches recomputation",
    "ROC sample size",
    "MATLAB eta definition present",
    "MATLAB 3 s window present",
    "PLV 0.55 percentile threshold present",
    "MATLAB ROC source phase pair"
  ),
  status = c(
    ifelse(length(missing_files) == 0, "PASS", "FAIL"),
    ifelse(nrow(stage_raw) == 24 * 5, "PASS", "FAIL"),
    ifelse(length(unique(stage_raw$case_id)) == 24, "PASS", "FAIL"),
    ifelse(all(table(stage_raw$case_id) == 5) && setequal(unique(stage_raw$phase5), phase_levels), "PASS", "FAIL"),
    ifelse(nrow(window_raw) == sum(stage_raw$n_windows), "PASS", "FAIL"),
    ifelse(max(c(stage_joined$diff_orig, stage_joined$diff_cov, stage_joined$diff_deg), na.rm = TRUE) < 1e-12 &&
             all(stage_joined$n_windows_diff == 0), "PASS", "FAIL"),
    ifelse(max(c(variance_joined$diff_orig, variance_joined$diff_cov, variance_joined$diff_deg), na.rm = TRUE) < 1e-15, "PASS", "FAIL"),
    ifelse(all(roc_stats$n_pre == 24 & roc_stats$n_early == 24), "PASS", "FAIL"),
    ifelse(grepl("eta_val = \\(sum_W / n\\) / lambda_1", ablation_code, fixed = FALSE), "PASS", "FAIL"),
    ifelse(grepl("winLenSec = 3", ablation_code, fixed = TRUE), "PASS", "FAIL"),
    ifelse(grepl("0.55", gain_code, fixed = TRUE), "PASS", "FAIL"),
    ifelse(grepl("pre-ictal", roc_code, fixed = TRUE) && grepl("early", roc_code, fixed = TRUE), "PASS", "FAIL")
  ),
  detail = c(
    paste(required_files, collapse = " | "),
    paste0("n=", nrow(stage_raw), "; expected=120"),
    paste0("n=", length(unique(stage_raw$case_id)), "; expected=24"),
    paste0("phase levels=", paste(unique(stage_raw$phase5), collapse = ", ")),
    paste0("window rows=", nrow(window_raw), "; stage n_windows sum=", sum(stage_raw$n_windows)),
    paste0("max eta diff=", formatC(max(c(stage_joined$diff_orig, stage_joined$diff_cov, stage_joined$diff_deg), na.rm = TRUE),
                                    format = "e", digits = 3),
           "; max n_windows diff=", max(abs(stage_joined$n_windows_diff), na.rm = TRUE)),
    paste0("max variance diff=", formatC(max(c(variance_joined$diff_orig, variance_joined$diff_cov, variance_joined$diff_deg), na.rm = TRUE),
                                         format = "e", digits = 3)),
    paste0(paste0(as.character(roc_stats$model), ": ", roc_stats$n_pre, " pre/", roc_stats$n_early, " early"), collapse = "; "),
    "compare_ablation_models_eta.m defines eta as (sum_W/n)/lambda_1",
    "compare_ablation_models_eta.m uses winLenSec = 3",
    "gain_hyperEdges_23.m uses the 0.55 PLV percentile rule",
    "plot_ROC_7_eta.m compares pre-ictal versus early with early as positive class"
  )
)
write.csv(pre_checks, file.path(out_dir, "Fig7_preplot_consistency_checks.csv"), row.names = FALSE)
if (any(pre_checks$status == "FAIL")) {
  stop("Pre-plot checks failed. Inspect Fig7_preplot_consistency_checks.csv")
}

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
      legend.key.width = unit(4.5, "mm"),
      legend.key.height = unit(3.2, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 0.1, face = "bold", colour = "#222222"),
      plot.title = element_text(size = base_size + 0.8, face = "bold", hjust = 0, colour = "#202124"),
      plot.subtitle = element_text(size = base_size - 0.2, hjust = 0, colour = "#4A4A4A", margin = margin(b = 2)),
      plot.margin = margin(4, 5, 4, 5),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

sig_short <- paired_stats %>%
  filter(comparison == "pre-ictal vs early") %>%
  mutate(
    label = case_when(
      p_value < 0.001 ~ "Pre-Early\nP<0.001",
      TRUE ~ paste0("Pre-Early\nP=", formatC(p_value, format = "f", digits = 3))
    )
  ) %>%
  left_join(stage_summary %>%
              group_by(model) %>%
              summarise(y = max(max_eta, na.rm = TRUE) + 0.0045, .groups = "drop"),
            by = "model")

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
    aes(x = 3, y = y, label = label, colour = model),
    inherit.aes = FALSE,
    size = 2.0,
    lineheight = 0.88,
    family = "Arial"
  ) +
  facet_wrap(~ model, nrow = 1) +
  scale_colour_manual(values = model_palette) +
  coord_cartesian(ylim = c(0.89, 1.008), clip = "off") +
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
  theme_nature() +
  theme(
    legend.position = c(0.58, 0.86),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.86), colour = NA),
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

base_file <- file.path(out_dir, "Fig7_eta_ablation_components")
save_pub(fig, base_file)

export_files <- paste0(base_file, c(".svg", ".pdf", ".png", ".tiff"))
post_checks <- tibble(
  file = basename(export_files),
  exists = file.exists(export_files),
  size_bytes = ifelse(file.exists(export_files), file.info(export_files)$size, NA_real_),
  status = ifelse(file.exists(export_files) & file.info(export_files)$size > 10000, "PASS", "FAIL")
)
write.csv(post_checks, file.path(out_dir, "Fig7_postplot_export_checks.csv"), row.names = FALSE)
if (any(post_checks$status == "FAIL")) {
  stop("Export checks failed. Inspect Fig7_postplot_export_checks.csv")
}

fmt <- function(x, digits = 4) formatC(x, format = "f", digits = digits)
fmt_p <- function(p) ifelse(p < 0.001, "P < 0.001", paste0("P = ", formatC(p, format = "f", digits = 3)))
get_stat <- function(model, comparison, col) {
  paired_stats %>%
    filter(as.character(.data$model) == .env$model, .data$comparison == .env$comparison) %>%
    pull(all_of(col)) %>%
    .[1]
}
get_auc <- function(model) roc_stats %>% filter(as.character(.data$model) == .env$model) %>% pull(auc) %>% .[1]
get_var <- function(model, phase) variance_calc %>%
  filter(as.character(.data$model) == .env$model, as.character(.data$phase5) == .env$phase) %>%
  pull(var_eta) %>%
  .[1]

results_text <- paste0(
  "DMWA components provide complementary sensitivity and stability\n",
  "To determine how the two weighting terms in DMWA contribute to the eta signal, we recomputed eta after selectively retaining cross-order Coverage or same-order DegreeBias. All three matrices preserved the seizure-stage increase from the pre-ictal period to early ictal activity, indicating that neither component alone abolishes the propagation-medium signature. The mean pre-ictal-to-early eta increase was ",
  fmt(get_stat("Original", "pre-ictal vs early", "mean_delta_eta"), 4), " for the full DMWA (", fmt_p(get_stat("Original", "pre-ictal vs early", "p_value")), "), ",
  fmt(get_stat("Coverage only", "pre-ictal vs early", "mean_delta_eta"), 4), " for Coverage only (", fmt_p(get_stat("Coverage only", "pre-ictal vs early", "p_value")), ") and ",
  fmt(get_stat("DegreeBias only", "pre-ictal vs early", "mean_delta_eta"), 4), " for DegreeBias only (", fmt_p(get_stat("DegreeBias only", "pre-ictal vs early", "p_value")), "; Fig. 7a). Thus, both cross-order overlap constraints and same-order local dominance contain seizure-stage information.\n\n",
  "The ablations differed more clearly in stability. Cross-seizure variance was lowest for the full DMWA and highest for the DegreeBias-only matrix at every stage (for example, early-stage variance: full DMWA ",
  fmt(get_var("Original", "early"), 6), ", Coverage only ", fmt(get_var("Coverage only", "early"), 6), ", DegreeBias only ", fmt(get_var("DegreeBias only", "early"), 6), "; Fig. 7b). This pattern supports the interpretation that Coverage contributes a stabilizing geometric constraint, whereas DegreeBias amplifies locally dominant same-order structure and increases inter-seizure heterogeneity. In pre-ictal versus early discrimination, performance was similar but slightly highest for the full DMWA (AUC = ",
  fmt(get_auc("Original"), 3), "), followed by Coverage only (AUC = ", fmt(get_auc("Coverage only"), 3), ") and DegreeBias only (AUC = ", fmt(get_auc("DegreeBias only"), 3), "; Fig. 7c). These results argue for the combined DMWA formulation as a balance between sensitivity to seizure-stage reorganization and stability across seizures."
)
writeLines(results_text, file.path(out_dir, "Fig7_results_description_draft.txt"), useBytes = TRUE)

legend_text <- paste0(
  "Fig. 7 | DMWA components provide complementary sensitivity and stability. ",
  "a, Stage-resolved effective refractive index (eta) after component ablation. Original denotes the full DMWA weighting, which combines cross-order Coverage with same-order DegreeBias; Coverage only retains the cross-order overlap/containment term; DegreeBias only retains the same-order local-dominance term. Boxes show the median and interquartile range across n = 24 seizures, whiskers show 1.5 times the interquartile range, dots denote individual seizures, and overlaid points with error bars show the mean and 95% confidence interval. Pre-ictal-to-early comparisons use paired two-sided Wilcoxon signed-rank tests. ",
  "b, Cross-seizure variance of eta for each ablation model and stage, computed from seizure-level stage means. ",
  "c, Receiver-operating characteristic curves for discriminating early ictal from pre-ictal stage summaries using eta. AUC values are computed from n = 24 pre-ictal and n = 24 early seizure-level observations per model; early is the positive class. Source data are seizure-level eta summaries derived from 3-s windows and the 0.55 PLV percentile hyperedge-extraction threshold."
)
writeLines(legend_text, file.path(out_dir, "Fig7_legend_draft.txt"), useBytes = TRUE)

qa_notes <- paste0(
  "Fig. 7 QA notes\n\n",
  "Figure contract:\n",
  "Core conclusion: cross-order Coverage and same-order DegreeBias each preserve the early ictal eta increase, but the full DMWA provides the best balance between early-stage sensitivity and cross-seizure stability.\n",
  "Evidence chain: panel a tests stage sensitivity for full and ablated matrices; panel b tests stability as cross-seizure variance; panel c tests early pre-ictal versus ictal discrimination.\n",
  "Archetype: quantitative grid with one full-width distribution panel and two compact validation panels.\n",
  "Backend: R only; ggplot2 + patchwork + svglite/cairo_pdf/ragg.\n",
  "Export: double-column SVG/PDF/TIFF/PNG with editable vector output.\n",
  "Palette: restrained neutral/blue/brown model palette, shared across all panels and aligned with the low-saturation NC style used in preceding figures.\n\n",
  "Data consistency:\n",
  "- Stage source rows: ", nrow(stage_raw), " rows = 24 seizures x 5 stages.\n",
  "- Window source rows: ", nrow(window_raw), "; summed stage n_windows: ", sum(stage_raw$n_windows), ".\n",
  "- Window-to-stage aggregation check maximum eta difference: ",
  formatC(max(c(stage_joined$diff_orig, stage_joined$diff_cov, stage_joined$diff_deg), na.rm = TRUE), format = "e", digits = 3), ".\n",
  "- Variance table recomputation maximum absolute difference: ",
  formatC(max(c(variance_joined$diff_orig, variance_joined$diff_cov, variance_joined$diff_deg), na.rm = TRUE), format = "e", digits = 3), ".\n",
  "- ROC comparison: pre-ictal versus early; early is the positive class; n = 24 per class per model.\n",
  "- AUC values: Original ", fmt(get_auc("Original"), 3), ", Coverage only ", fmt(get_auc("Coverage only"), 3), ", DegreeBias only ", fmt(get_auc("DegreeBias only"), 3), ".\n",
  "- Source-code audit: compare_ablation_models_eta.m defines eta as (sum_W/n)/lambda_1, uses 3-s windows, and constructs Original, Coverage-only and DegreeBias-only matrices from Coverage and DegreeBias terms. The hyperedge extraction function uses the 0.55 PLV percentile threshold.\n\n",
  "Interpretation boundary:\n",
  "This figure supports component-level sensitivity/stability differences within DMWA. It should not be used to claim that ablation alone proves ordered topology; that conclusion belongs to the density and randomization controls.\n\n",
  "Working memory note:\n",
  "For this manuscript workflow, do not modify Word files without explicit user approval. Provide additions or replacements as text for the user to paste manually. Fig. 7 uses data from the local 2_Fig7 directory and outputs to the manuscript nature_fig/Fig_7 directory."
)
writeLines(qa_notes, file.path(out_dir, "Fig7_QA_notes.txt"), useBytes = TRUE)

message("Fig. 7 complete: ", base_file)
