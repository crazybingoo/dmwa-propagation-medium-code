#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(svglite)
  library(ragg)
})

FIG_TEXT_PT <- 6
FIG_PANEL_PT <- 8
FIG_GEOM_TEXT_SIZE <- FIG_TEXT_PT / 2.845276


out_dir_candidates <- Sys.glob(file.path("F:/6", "*0514-", "nature_fig", "Fig_6"))
if (length(out_dir_candidates) < 1) {
  stop("Cannot locate Fig_6 directory under F:/6/*0514-/nature_fig")
}
out_dir <- normalizePath(out_dir_candidates[1], winslash = "/", mustWork = TRUE)
source_csv <- file.path(out_dir, "Fig6_robustness_source_data.csv")
stage_csv <- file.path(out_dir, "Fig6_robustness_stage_means.csv")

if (!file.exists(source_csv)) {
  stop("Missing source data: ", source_csv)
}
if (!file.exists(stage_csv)) {
  stop("Missing stage means: ", stage_csv)
}

df <- read.csv(source_csv, check.names = FALSE)
stage_df <- read.csv(stage_csv, check.names = FALSE)

df <- df %>%
  mutate(
    window_f = factor(window_s, levels = sort(unique(window_s)), labels = paste0(sort(unique(window_s)), " s")),
    threshold_f = factor(plv_threshold, levels = sort(unique(plv_threshold))),
    p_label = case_when(
      p_value < 0.001 ~ "P<0.001",
      TRUE ~ paste0("P=", formatC(p_value, format = "f", digits = 3))
    ),
    p_num = formatC(neg_log10_p, format = "f", digits = 2),
    diff_num = formatC(mean_ictal_minus_pre, format = "f", digits = 4)
  )

stage_df <- stage_df %>%
  mutate(
    window_f = factor(window_s, levels = sort(unique(window_s)), labels = paste0(sort(unique(window_s)), " s")),
    threshold_f = factor(plv_threshold, levels = sort(unique(plv_threshold)))
  )

phase_order_ok <- all(sort(unique(df$window_s)) == c(1, 2, 3, 4, 5)) &&
  all(abs(sort(unique(df$plv_threshold)) - c(0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70)) < 1e-9)

pre_checks <- tibble::tibble(
  check = c(
    "source rows",
    "parameter grid",
    "positive effects",
    "significance p<0.05",
    "significance p<0.01",
    "screenshot value range"
  ),
  status = c(
    ifelse(nrow(df) == 35, "PASS", "FAIL"),
    ifelse(phase_order_ok, "PASS", "FAIL"),
    ifelse(all(df$mean_ictal_minus_pre > 0, na.rm = TRUE), "PASS", "FAIL"),
    ifelse(all(df$p_value < 0.05, na.rm = TRUE), "PASS", "FAIL"),
    ifelse(all(df$p_value < 0.01, na.rm = TRUE), "PASS", "FAIL"),
    ifelse(min(df$neg_log10_p) > 2.9 && max(df$neg_log10_p) > 4.3 &&
             min(df$mean_ictal_minus_pre) > 0.003 && max(df$mean_ictal_minus_pre) < 0.024, "PASS", "WARN")
  ),
  detail = c(
    paste0("n=", nrow(df)),
    paste0("windows=", paste(sort(unique(df$window_s)), collapse = ","),
           "; thresholds=", paste(formatC(sort(unique(df$plv_threshold)), format = "f", digits = 2), collapse = ",")),
    paste0("effect range=", formatC(min(df$mean_ictal_minus_pre), format = "f", digits = 6),
           " to ", formatC(max(df$mean_ictal_minus_pre), format = "f", digits = 6)),
    paste0(sum(df$p_value < 0.05), " of ", nrow(df), " settings"),
    paste0(sum(df$p_value < 0.01), " of ", nrow(df), " settings"),
    paste0("-log10(P) range=", formatC(min(df$neg_log10_p), format = "f", digits = 3),
           " to ", formatC(max(df$neg_log10_p), format = "f", digits = 3))
  )
)
write.csv(pre_checks, file.path(out_dir, "Fig6_preplot_consistency_checks.csv"), row.names = FALSE)

if (any(pre_checks$status == "FAIL")) {
  stop("Pre-plot checks failed; inspect Fig6_preplot_consistency_checks.csv")
}

theme_heat <- function(base_size = FIG_TEXT_PT, base_family = "Arial") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      axis.title = element_text(size = FIG_TEXT_PT),
      axis.text = element_text(size = FIG_TEXT_PT, colour = "#252525"),
      plot.title = element_text(size = FIG_TEXT_PT, face = "bold", hjust = 0, colour = "#202124"),
      panel.grid = element_blank(),
      legend.title = element_text(size = FIG_TEXT_PT),
      legend.text = element_text(size = FIG_TEXT_PT),
      legend.key.height = unit(14, "mm"),
      plot.margin = margin(5, 6, 5, 5),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

text_colour_for <- function(x, threshold) {
  ifelse(x >= threshold, "white", "#151515")
}

# Low-saturation NC-style sequential palette aligned with the preceding
# DMWA figures. Both heatmaps use the same blue-grey family; the colour-bar
# titles distinguish statistical strength from effect magnitude.
pal_sig <- c("#F3F6F8", "#DDE7EF", "#B9CCD9", "#7EA4BE", "#477DA0", "#245A78")

p_a <- ggplot(df, aes(x = threshold_f, y = window_f, fill = neg_log10_p)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = p_num, colour = text_colour_for(neg_log10_p, 3.78)), size = FIG_GEOM_TEXT_SIZE, family = "Arial") +
  scale_colour_identity() +
  scale_fill_gradientn(
    colours = pal_sig,
    name = "-log10(P)",
    limits = c(2.9, 4.4),
    breaks = c(3.0, 3.5, 4.0),
    values = scales::rescale(c(2.9, 3.2, 3.5, 3.8, 4.1, 4.4)),
    guide = guide_colourbar(frame.colour = NA, ticks.colour = "#333333")
  ) +
  scale_y_discrete(limits = rev(levels(df$window_f))) +
  labs(
    x = "PLV threshold percentile",
    y = "Window length",
    title = expression("Discriminative power of " * eta)
  ) +
  theme_heat()

p_b <- ggplot(df, aes(x = threshold_f, y = window_f, fill = mean_ictal_minus_pre)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = diff_num, colour = text_colour_for(mean_ictal_minus_pre, 0.014)), size = FIG_GEOM_TEXT_SIZE, family = "Arial") +
  scale_colour_identity() +
  scale_fill_gradientn(
    colours = pal_sig,
    name = expression(Delta * eta),
    limits = c(0.003, 0.0235),
    breaks = c(0.005, 0.010, 0.015, 0.020),
    values = scales::rescale(c(0.003, 0.006, 0.010, 0.014, 0.018, 0.0235)),
    guide = guide_colourbar(frame.colour = NA, ticks.colour = "#333333")
  ) +
  scale_y_discrete(limits = rev(levels(df$window_f))) +
  labs(
    x = "PLV threshold percentile",
    y = NULL,
    title = expression("Effect size: mean ictal " * eta * " - mean pre-ictal " * eta)
  ) +
  theme_heat()

fig <- (p_a | p_b) +
  plot_layout(widths = c(1, 1), guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = FIG_PANEL_PT, face = "bold", family = "Arial"))

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 82, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()
  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = 300)
  print(plot)
  dev.off()
}

out_base <- file.path(out_dir, "Fig6_eta_robustness_font6_label8")
save_pub_r(fig, out_base)

post_checks <- tibble::tibble(
  check = c("svg export", "pdf export", "tiff export", "png export", "source data", "stage means", "preplot consistency report"),
  path = c(
    paste0(out_base, ".svg"),
    paste0(out_base, ".pdf"),
    paste0(out_base, ".tiff"),
    paste0(out_base, ".png"),
    source_csv,
    stage_csv,
    file.path(out_dir, "Fig6_preplot_consistency_checks.csv")
  )
) %>%
  mutate(
    exists = file.exists(path),
    bytes = ifelse(exists, file.info(path)$size, NA_real_),
    status = ifelse(exists & bytes > 0, "PASS", "FAIL")
  )
write.csv(post_checks, file.path(out_dir, "Fig6_postplot_export_checks.csv"), row.names = FALSE)

best_p <- df %>% slice_max(order_by = neg_log10_p, n = 1, with_ties = FALSE)
max_effect <- df %>% slice_max(order_by = mean_ictal_minus_pre, n = 1, with_ties = FALSE)
min_effect <- df %>% slice_min(order_by = mean_ictal_minus_pre, n = 1, with_ties = FALSE)

qa_lines <- c(
  "Fig. 6 QA notes",
  "",
  "Figure contract:",
  "Core conclusion: the DMWA-derived effective refractive index eta remains higher during ictal than pre-ictal periods across the tested window lengths and PLV thresholds.",
  "Archetype: quantitative grid with paired heatmaps.",
  "Backend: R only; ggplot2 + patchwork + svglite/cairo_pdf/ragg.",
  "Export: double-column SVG/PDF/TIFF/PNG.",
  "Palette: low-saturation NC-style blue-grey gradient shared by both heatmaps and aligned with the preceding DMWA figures.",
  "",
  "Data consistency:",
  paste0("- Source grid: ", nrow(df), " parameter settings: window lengths ",
         paste(sort(unique(df$window_s)), collapse = ", "),
         " s and PLV thresholds ",
         paste(formatC(sort(unique(df$plv_threshold)), format = "f", digits = 2), collapse = ", "), "."),
  paste0("- All effects are positive: range ", formatC(min(df$mean_ictal_minus_pre), format = "f", digits = 6),
         " to ", formatC(max(df$mean_ictal_minus_pre), format = "f", digits = 6), "."),
  paste0("- All settings are significant at P < 0.01; ", sum(df$p_value < 0.001),
         " of ", nrow(df), " settings are significant at P < 0.001."),
  paste0("- Highest -log10(P): window ", best_p$window_s, " s, threshold ",
         formatC(best_p$plv_threshold, format = "f", digits = 2), ", value ",
         formatC(best_p$neg_log10_p, format = "f", digits = 3), "."),
  paste0("- Largest effect: window ", max_effect$window_s, " s, threshold ",
         formatC(max_effect$plv_threshold, format = "f", digits = 2), ", Delta eta ",
         formatC(max_effect$mean_ictal_minus_pre, format = "f", digits = 6), "."),
  paste0("- Smallest effect: window ", min_effect$window_s, " s, threshold ",
         formatC(min_effect$plv_threshold, format = "f", digits = 2), ", Delta eta ",
         formatC(min_effect$mean_ictal_minus_pre, format = "f", digits = 6), "."),
  "",
  "Interpretation boundary:",
  "This figure supports parameter robustness for the ictal-versus-pre-ictal eta increase. It does not test topology-specific controls, ablations, or low-order baselines."
)
writeLines(qa_lines, file.path(out_dir, "Fig6_QA_notes.txt"))

legend_lines <- c(
  "Fig. 6 | η dynamics are robust across window lengths and PLV thresholds.",
  "a, Discriminative strength of the ictal-versus-pre-ictal η difference across analysis parameters, shown as -log10(P). b, Corresponding effect size, computed as mean ictal η minus mean pre-ictal η. Rows denote window length and columns denote PLV threshold percentile. Values inside cells denote the plotted statistic. Source data are from n = 24 seizures; tests compare paired ictal and pre-ictal seizure summaries for each parameter setting."
)
writeLines(legend_lines, file.path(out_dir, "Fig6_legend_draft.txt"))

result_lines <- c(
  "η dynamics are robust across analysis parameters",
  "",
  "We next examined whether the seizure-associated increase in η depended on the analysis window length or the PLV threshold used to construct higher-order interaction units. We repeated the ictal-versus-pre-ictal comparison across 35 parameter settings, spanning window lengths from 1 to 5 s and PLV threshold percentiles from 0.40 to 0.70 (Fig. 6).",
  "",
  "The η increase was preserved throughout this parameter grid. The mean ictal-minus-pre-ictal difference was positive for every parameter setting, ranging from 0.00317 to 0.02308, and all paired comparisons were significant at P < 0.01. The strongest statistical evidence occurred at longer windows and higher PLV thresholds, with -log10(P) reaching 4.36, whereas larger effect sizes were observed at lower thresholds and shorter windows. This pattern indicates that η robustly separates ictal from pre-ictal states across the tested parameter range, while the balance between effect magnitude and statistical stability varies with window length and threshold choice.",
  "",
  "Boundary: this robustness analysis supports parameter stability of the η stage effect; it does not by itself establish the topology-specific mechanism tested in the control analyses."
)
writeLines(result_lines, file.path(out_dir, "Fig6_results_description_draft.txt"))

message("Fig. 6 complete: ", out_dir)
