#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(svglite)
  library(ragg)
  library(ggtext)
})

FIG_TEXT_PT <- 10.5
FIG_TITLE_PT <- 10.5
FIG_TICK_PT <- 10.5
FIG_X_TICK_PT <- 9.4
FIG_LEGEND_PT <- 9.2
FIG_CBAR_PT <- 8.2
FIG_PANEL_PT <- 12.5
FIG_HEATMAP_TEXT_PT <- 8.0
FIG_GEOM_TEXT_SIZE <- FIG_TEXT_PT / 2.845276
FIG_HEATMAP_TEXT_SIZE <- FIG_HEATMAP_TEXT_PT / 2.845276
FIG_X_ANGLE <- 35

out_dir <- file.path("figures", "Fig_6")
source_csv <- file.path(out_dir, "Fig6_retained_fraction_source_data.csv")
stage_csv <- file.path(out_dir, "Fig6_retained_fraction_stage_means.csv")

if (!file.exists(source_csv)) {
  stop("Missing retained-fraction source data: ", source_csv)
}
if (!file.exists(stage_csv)) {
  stop("Missing retained-fraction stage means: ", stage_csv)
}

df <- read.csv(source_csv, check.names = FALSE)
stage_df <- read.csv(stage_csv, check.names = FALSE)

expected_windows <- 1:5
expected_fraction <- seq(0.45, 0.90, by = 0.05)

df <- df %>%
  mutate(
    window_f = factor(WindowSec, levels = rev(expected_windows), labels = paste0(rev(expected_windows), " s")),
    window_line = factor(WindowSec, levels = expected_windows, labels = paste0(expected_windows, " s")),
    retained_f = factor(sprintf("%.2f", RetainedFraction), levels = sprintf("%.2f", expected_fraction)),
    p_text = sprintf("%.2f", NegLog10P),
    diff_text = sprintf("%.3f", MeanIctalMinusPre),
    sig_class = case_when(
      PValue < 0.01 ~ "P<0.01",
      PValue < 0.05 ~ "P<0.05",
      TRUE ~ "n.s."
    )
  )

stage_df <- stage_df %>%
  mutate(
    retained_f = factor(sprintf("%.2f", RetainedFraction), levels = sprintf("%.2f", expected_fraction)),
    window_line = factor(WindowSec, levels = expected_windows, labels = paste0(expected_windows, " s"))
  )

grid_ok <- nrow(df) == length(expected_windows) * length(expected_fraction) &&
  all(sort(unique(df$WindowSec)) == expected_windows) &&
  all(abs(sort(unique(df$RetainedFraction)) - expected_fraction) < 1e-9) &&
  all(df$NValidPairs == 24)

if (!grid_ok) {
  stop("Retained-fraction Fig. 6 source grid failed validation.")
}

pal <- c(
  ink = "#202124",
  grid = "#ECEFF3",
  neutral = "#7A828C",
  neutral_light = "#D9DEE7",
  blue_0 = "#F5F8FB",
  blue_1 = "#DDE7F0",
  blue_2 = "#B8CCDC",
  blue_3 = "#7FA8C9",
  blue_4 = "#5F8DBB",
  blue_5 = "#4778A8",
  summary_blue = "#4778A8",
  summary_band = "#C8D3E0",
  orange = "#E6A34A",
  red = "#D96661",
  purple = "#8D7BB8"
)

window_cols <- c(
  "1 s" = "#5D83B5",
  "2 s" = "#6EAD67",
  "3 s" = "#E6A34A",
  "4 s" = "#D96661",
  "5 s" = "#8D7BB8"
)

theme_fig <- function() {
  theme_classic(base_size = FIG_TEXT_PT, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = pal["ink"]),
      axis.ticks = element_line(linewidth = 0.35, colour = pal["ink"]),
      axis.title = ggtext::element_markdown(size = FIG_TEXT_PT, colour = pal["ink"]),
      axis.title.x = ggtext::element_markdown(margin = margin(t = 2.0)),
      axis.title.y = ggtext::element_markdown(margin = margin(r = 1.5)),
      axis.text = element_text(size = FIG_TICK_PT, colour = pal["ink"]),
      plot.title = ggtext::element_markdown(size = FIG_TITLE_PT, face = "bold", hjust = 0, colour = pal["ink"], margin = margin(b = 3.0)),
      legend.title = element_text(size = FIG_LEGEND_PT, colour = pal["ink"]),
      legend.text = element_text(size = FIG_LEGEND_PT, colour = pal["ink"]),
      legend.key.height = unit(2.4, "mm"),
      legend.key.width = unit(3.6, "mm"),
      legend.background = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(0.35, "mm"),
      legend.spacing.x = unit(0.8, "mm"),
      legend.spacing.y = unit(0.4, "mm"),
      plot.margin = margin(3.2, 3.2, 3.2, 3.2),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

theme_heat <- function() {
  theme_fig() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = FIG_X_TICK_PT, angle = FIG_X_ANGLE, hjust = 1, vjust = 1),
      legend.position = "right",
      legend.justification = "center",
      legend.key.height = unit(2.0, "mm")
    )
}

text_colour_for <- function(x, threshold) {
  ifelse(x >= threshold, "white", pal["ink"])
}

sig_points <- df %>%
  mutate(sig_class = factor(sig_class, levels = c("P<0.01", "P<0.05", "n.s.")))

profile <- df %>%
  group_by(RetainedFraction, retained_f) %>%
  summarise(
    mean_delta = mean(MeanIctalMinusPre),
    min_delta = min(MeanIctalMinusPre),
    max_delta = max(MeanIctalMinusPre),
    min_neglog10p = min(NegLog10P),
    n_p01 = sum(PValue < 0.01),
    n_p05 = sum(PValue < 0.05),
    .groups = "drop"
  )

best_p <- df %>% slice_max(order_by = NegLog10P, n = 1, with_ties = FALSE)
max_effect <- df %>% slice_max(order_by = MeanIctalMinusPre, n = 1, with_ties = FALSE)
min_effect <- df %>% slice_min(order_by = MeanIctalMinusPre, n = 1, with_ties = FALSE)

p_a <- ggplot(df, aes(retained_f, window_f, fill = NegLog10P)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  scale_fill_gradientn(
    colours = c(pal["blue_0"], pal["blue_1"], pal["blue_2"], pal["blue_3"], pal["blue_4"], pal["blue_5"]),
    name = "\u2212log\u2081\u2080(\U0001D443)",
    limits = c(0, 4.5),
    breaks = c(1, 2, 3, 4),
    values = rescale(c(0, 1, 2, 3, 4, 4.5)),
    guide = guide_colourbar(
      frame.colour = NA,
      ticks.colour = pal["ink"],
      barwidth = unit(1.45, "mm"),
      barheight = unit(34, "mm"),
      title.position = "right",
      label.theme = element_text(size = FIG_CBAR_PT),
      title.theme = element_text(size = FIG_CBAR_PT, angle = 90, hjust = 0.5, margin = margin(l = 0.6))
    )
  ) +
  labs(
    title = "Statistical support for <i>\u03b7</i> increase",
    x = "Retained PLV edge fraction",
    y = "Window length"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed(ratio = 1, clip = "off") +
  theme_heat()

p_b <- ggplot(df, aes(retained_f, window_f, fill = MeanIctalMinusPre)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  scale_fill_gradientn(
    colours = c(pal["blue_0"], pal["blue_1"], pal["blue_2"], pal["blue_3"], pal["blue_4"], pal["blue_5"]),
    name = "\U0001D6E5\U0001D702",
    limits = c(0, 0.018),
    breaks = c(0.004, 0.008, 0.012, 0.016),
    labels = number_format(accuracy = 0.001),
    values = rescale(c(0, 0.0025, 0.006, 0.010, 0.014, 0.018)),
    guide = guide_colourbar(
      frame.colour = NA,
      ticks.colour = pal["ink"],
      barwidth = unit(1.45, "mm"),
      barheight = unit(34, "mm"),
      title.position = "right",
      label.theme = element_text(size = FIG_CBAR_PT),
      title.theme = element_text(size = FIG_CBAR_PT, angle = 90, hjust = 0.5, margin = margin(l = 0.6))
    )
  ) +
  labs(
    title = "Effect size, <i>\u0394\u03b7</i>",
    x = "Retained PLV edge fraction",
    y = NULL
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed(ratio = 1, clip = "off") +
  theme_heat()

p_c <- ggplot() +
  annotate("rect", xmin = 0.425, xmax = 0.825, ymin = -Inf, ymax = Inf,
           fill = pal["blue_1"], alpha = 0.22, colour = NA) +
  annotate("rect", xmin = 0.825, xmax = 0.875, ymin = -Inf, ymax = Inf,
           fill = pal["orange"], alpha = 0.14, colour = NA) +
  geom_line(data = df, aes(RetainedFraction, MeanIctalMinusPre, group = window_line, colour = window_line),
            linewidth = 0.45, alpha = 0.58) +
  geom_point(data = df, aes(RetainedFraction, MeanIctalMinusPre, colour = window_line),
             size = 1.45, alpha = 0.78) +
  geom_ribbon(data = profile, aes(RetainedFraction, ymin = min_delta, ymax = max_delta),
              inherit.aes = FALSE, fill = pal["summary_band"], alpha = 0.40) +
  geom_line(data = profile, aes(RetainedFraction, mean_delta),
            inherit.aes = FALSE, linewidth = 0.85, colour = pal["summary_blue"]) +
  geom_point(data = profile, aes(RetainedFraction, mean_delta),
             inherit.aes = FALSE, size = 1.8, shape = 21, stroke = 0.35,
             fill = "white", colour = pal["summary_blue"]) +
  scale_colour_manual(
    values = window_cols,
    name = "Window",
    guide = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(linewidth = 0.45, size = 1.8)
    )
  ) +
  scale_x_continuous(breaks = expected_fraction, labels = number_format(accuracy = 0.01), expand = expansion(mult = c(0.025, 0.025))) +
  scale_y_continuous(labels = number_format(accuracy = 0.001), expand = expansion(mult = c(0.04, 0.08))) +
  labs(
    title = "<i>\u0394\u03b7</i> attenuation with scaffold density",
    x = "Retained PLV edge fraction",
    y = "Ictal - pre-ictal <i>\u03b7</i>"
  ) +
  theme_fig() +
  theme(
    axis.text.x = element_text(size = FIG_X_TICK_PT, angle = FIG_X_ANGLE, hjust = 1, vjust = 1),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title.position = "left",
    legend.box.spacing = unit(0.1, "mm")
  )

p_d <- ggplot(sig_points, aes(retained_f, window_f)) +
  geom_tile(fill = "white", colour = pal["grid"], linewidth = 0.32) +
  geom_point(aes(colour = sig_class, shape = sig_class), size = 2.45, stroke = 0.72) +
  scale_colour_manual(
    values = c("P<0.01" = unname(pal["summary_blue"]), "P<0.05" = unname(pal["orange"]), "n.s." = unname(pal["neutral"])),
    breaks = c("P<0.01", "P<0.05", "n.s."),
    labels = c("\U0001D443 < 0.01", "\U0001D443 < 0.05", "n.s."),
    name = "Paired test",
    guide = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(shape = c(16, 1, 4), size = 2.5, stroke = 0.72)
    )
  ) +
  scale_shape_manual(
    values = c("P<0.01" = 16, "P<0.05" = 1, "n.s." = 4),
    breaks = c("P<0.01", "P<0.05", "n.s."),
    labels = c("\U0001D443 < 0.01", "\U0001D443 < 0.05", "n.s."),
    name = "Paired test",
    guide = "none"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = "Significance tier",
    x = "Retained PLV edge fraction",
    y = "Window length"
  ) +
  theme_fig() +
  theme(
    axis.text.x = element_text(size = FIG_X_TICK_PT, angle = FIG_X_ANGLE, hjust = 1, vjust = 1),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title.position = "left",
    legend.box.spacing = unit(0.1, "mm"),
    plot.margin = margin(3.2, 0.4, 3.2, 3.2)
  )

p_a <- p_a + labs(tag = "a")
p_b <- p_b + labs(tag = "b")
p_c <- p_c + labs(tag = "c")
p_d <- p_d + labs(tag = "d")

top_row <- (p_a | p_b) +
  plot_layout(widths = c(1, 1), guides = "keep")

bottom_row <- (p_c | plot_spacer() | p_d) +
  plot_layout(widths = c(1, 0.136, 1), guides = "keep")

fig <- ((top_row / bottom_row) +
  plot_layout(heights = c(1, 1), guides = "keep")) &
  theme(
    plot.tag = element_text(size = FIG_PANEL_PT, face = "bold", family = "Arial", colour = pal["ink"], hjust = 0, vjust = 1),
    plot.tag.position = c(0, 1)
  )

save_pub_r <- function(plot, filename, width_mm = 178, height_mm = 126, dpi = 600) {
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

out_base <- file.path(out_dir, "Fig6_eta_robustness_retained_fraction")
save_pub_r(fig, out_base)

pre_checks <- data.frame(
  check = c(
    "source rows",
    "parameter grid",
    "valid paired seizures",
    "positive effects",
    "P < 0.01 settings",
    "P < 0.05 settings",
    "highest statistical support",
    "largest effect",
    "smallest effect"
  ),
  value = c(
    nrow(df),
    paste0("windows=", paste(expected_windows, collapse = ","), "; retained fractions=", paste(sprintf("%.2f", expected_fraction), collapse = ",")),
    paste0(sum(df$NValidPairs == 24), "/", nrow(df)),
    paste0(sum(df$MeanIctalMinusPre > 0), "/", nrow(df)),
    paste0(sum(df$PValue < 0.01), "/", nrow(df)),
    paste0(sum(df$PValue < 0.05), "/", nrow(df)),
    paste0(best_p$WindowSec, " s, retained fraction ", sprintf("%.2f", best_p$RetainedFraction), ", -log10(P)=", sprintf("%.3f", best_p$NegLog10P)),
    paste0(max_effect$WindowSec, " s, retained fraction ", sprintf("%.2f", max_effect$RetainedFraction), ", Delta eta=", sprintf("%.6f", max_effect$MeanIctalMinusPre)),
    paste0(min_effect$WindowSec, " s, retained fraction ", sprintf("%.2f", min_effect$RetainedFraction), ", Delta eta=", sprintf("%.6f", min_effect$MeanIctalMinusPre))
  )
)
write.csv(pre_checks, file.path(out_dir, "Fig6_retained_fraction_preplot_checks.csv"), row.names = FALSE)

post_checks <- data.frame(
  check = c("svg export", "pdf export", "tiff export", "png export", "source data", "stage means"),
  path = c(paste0(out_base, ".svg"), paste0(out_base, ".pdf"), paste0(out_base, ".tiff"), paste0(out_base, ".png"), source_csv, stage_csv)
) %>%
  mutate(
    exists = file.exists(path),
    bytes = ifelse(exists, file.info(path)$size, NA_real_),
    status = ifelse(exists & bytes > 0, "PASS", "FAIL")
  )
write.csv(post_checks, file.path(out_dir, "Fig6_retained_fraction_postplot_export_checks.csv"), row.names = FALSE)

legend_text <- paste0(
  "Fig. 6 | ", "\u03b7", " dynamics are robust across window lengths and retained PLV edge fractions. ",
  "a, Discriminative strength of the ictal-versus-pre-ictal ", "\u03b7", " difference across analysis parameters, shown as \u2212log10(P). ",
  "b, Corresponding effect size, computed as mean ictal ", "\u03b7", " minus mean pre-ictal ", "\u03b7", ". ",
  "Rows denote window length and columns denote the fraction of strongest PLV edges retained before hyperedge extraction. ",
  "c, Mean effect profile across retained fractions. Coloured points and thin lines denote individual window lengths; the dark line shows the mean across window lengths and the shaded band shows the range. ",
  "d, Significance tier for each window length and retained fraction. ",
  "Source data are from n = 24 seizures; tests compare paired ictal and pre-ictal seizure summaries using two-sided Wilcoxon signed-rank tests."
)
writeLines(legend_text, file.path(out_dir, "Fig6_retained_fraction_legend_draft.txt"), useBytes = TRUE)

result_text <- paste0(
  "\u03b7", " dynamics are robust across retained PLV edge fractions\n\n",
  "We next examined whether the seizure-associated increase in ", "\u03b7", " depended on the analysis window length or on the density of the PLV scaffold used for hyperedge extraction. ",
  "Using the retained-edge implementation of the thresholding routine, we repeated the ictal-versus-pre-ictal comparison across 50 parameter settings, spanning window lengths from 1 to 5 s and retained PLV edge fractions from 0.45 to 0.90 (Fig. 6). ",
  "The ictal-minus-pre-ictal difference remained positive across all settings (", sprintf("%.6f", min(df$MeanIctalMinusPre)), " to ", sprintf("%.6f", max(df$MeanIctalMinusPre)), "). ",
  "The effect was significant at P < 0.01 in ", sum(df$PValue < 0.01), " of 50 settings and at P < 0.05 in ", sum(df$PValue < 0.05), " of 50 settings. ",
  "Across retained fractions from 0.45 to 0.80, all five window lengths remained significant at P < 0.01; at 0.85, all window lengths remained significant at P < 0.05 but only two reached P < 0.01; at 0.90, the effect was positive but no longer significant. ",
  "Thus, ", "\u03b7", " distinguishes ictal from pre-ictal states across a broad range of window lengths and scaffold densities, while very dense PLV scaffolds attenuate the effect."
)
writeLines(result_text, file.path(out_dir, "Fig6_retained_fraction_results_description_draft.txt"), useBytes = TRUE)

qa_lines <- c(
  "Fig. 6 retained-fraction QA notes",
  "Core conclusion: the DMWA-derived effective refractive index eta remains higher during ictal than pre-ictal periods across broad retained PLV edge fractions, but the effect attenuates at very dense scaffolds.",
  "Archetype: quantitative grid with paired heatmaps plus summary panels.",
  "Backend: R only; ggplot2 + patchwork + svglite/cairo_pdf/ragg.",
  "Source grid: 50 settings; window lengths 1, 2, 3, 4, 5 s; retained PLV edge fractions 0.45 to 0.90.",
  paste0("Valid pairs: ", sum(df$NValidPairs == 24), " of ", nrow(df), " settings have n = 24 paired seizures."),
  paste0("Positive effects: ", sum(df$MeanIctalMinusPre > 0), " of ", nrow(df), "."),
  paste0("P < 0.01: ", sum(df$PValue < 0.01), " of ", nrow(df), "; P < 0.05: ", sum(df$PValue < 0.05), " of ", nrow(df), "."),
  paste0("Effect range: ", sprintf("%.6f", min(df$MeanIctalMinusPre)), " to ", sprintf("%.6f", max(df$MeanIctalMinusPre)), "."),
  "Statistics: paired two-sided Wilcoxon signed-rank tests comparing seizure-level ictal and pre-ictal eta summaries.",
  "Interpretation boundary: the figure supports parameter robustness and shows attenuation at very dense scaffolds; it does not establish topology specificity by itself."
)
writeLines(qa_lines, file.path(out_dir, "Fig6_retained_fraction_QA_notes.txt"), useBytes = TRUE)

message("Fig. 6 retained-fraction figure exported to: ", out_dir)
