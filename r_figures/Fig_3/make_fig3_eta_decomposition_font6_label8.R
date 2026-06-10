suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggbeeswarm)
  library(svglite)
  library(ragg)
  library(scales)
  library(ggrepel)
  library(grid)
})

FIG_TEXT_PT <- 6
FIG_PANEL_PT <- 8
FIG_GEOM_TEXT_SIZE <- FIG_TEXT_PT / 2.845276


base_candidates <- Sys.glob("F:/6/*0514-/nature_fig/Fig_3")
if (length(base_candidates) < 1) {
  stop("Cannot locate Fig_3 directory under F:/6/*0514-/nature_fig/Fig_3")
}
base_dir <- normalizePath(base_candidates[1], winslash = "/", mustWork = TRUE)
out_base <- file.path(base_dir, "Fig3_eta_decomposition_font6_label8")

phase_levels <- c("pre-ictal", "early", "mid", "late", "post-ictal")
phase_labels <- c("Pre", "Early", "Mid", "Late", "Post")
phase_cols <- c(
  "pre-ictal" = "#5D83B5",
  "early" = "#6EAD67",
  "mid" = "#E6A34A",
  "late" = "#D96661",
  "post-ictal" = "#8D7BB8"
)
phase_change_levels <- phase_levels
phase_change_labels <- phase_labels
eta_sym <- "\u03b7"
lambda1_sym <- "\u03bb1"
lambda2_sym <- "\u03bb2"
times_sym <- "\u00d7"

theme_nature <- function(base_size = FIG_TEXT_PT) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#202124"),
      axis.ticks = element_line(linewidth = 0.30, colour = "#202124"),
      axis.text = element_text(size = FIG_TEXT_PT, colour = "#202124"),
      axis.title = element_text(size = FIG_TEXT_PT, colour = "#202124"),
      plot.title = element_text(size = FIG_TEXT_PT, face = "bold", hjust = 0),
      legend.title = element_text(size = FIG_TEXT_PT),
      legend.text = element_text(size = FIG_TEXT_PT),
      legend.key.height = unit(3.0, "mm"),
      legend.key.width = unit(4.5, "mm"),
      legend.background = element_blank(),
      panel.grid.major = element_line(linewidth = 0.18, colour = "#ECEFF3"),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = FIG_TEXT_PT, face = "bold", colour = "#202124"),
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

metric_spec <- tibble::tribble(
  ~metric_id,       ~column,                 ~metric_plot,              ~metric_plain,        ~family,
  "eta",            "eta_mean",              "eta",                     "eta",                "index",
  "resource_R",     "resource_R_mean",       "R",                       "R",                  "primary",
  "lambda1",        "lambda1_mean",          "u03bb1",               "lambda1",            "primary",
  "spectral_gap",   "spectral_gap_mean",     "u03bb1-u03bb2",     "lambda1-lambda2",    "secondary",
  "lambda2_ratio",  "lambda2_ratio_mean",    "u03bb2/u03bb1",     "lambda2/lambda1",    "secondary",
  "spectral_PR",    "spectral_PR_mean",      "Spectral~PR",             "Spectral PR",        "secondary"
)

metric_plot_levels <- metric_spec$metric_plot

stage <- read_csv(
  file.path(base_dir, "ALL_CASES_stage5_eta_decomposition.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    phase5 = factor(phase5, levels = phase_levels),
    phase_label = factor(phase_labels[as.integer(phase5)], levels = phase_labels)
  )

window_level <- read_csv(
  file.path(base_dir, "ALL_CASES_window_level_eta_decomposition.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    phase5 = factor(phase5, levels = phase_levels),
    phase_label = factor(phase_labels[as.integer(phase5)], levels = phase_labels)
  )

stage_long <- stage %>%
  select(case_id, phase5, phase_label, all_of(metric_spec$column)) %>%
  pivot_longer(
    cols = all_of(metric_spec$column),
    names_to = "column",
    values_to = "value"
  ) %>%
  left_join(metric_spec, by = "column") %>%
  mutate(
    metric_plot = factor(metric_plot, levels = rev(metric_plot_levels)),
    metric_plain = factor(metric_plain, levels = metric_spec$metric_plain)
  )

stage_summary <- stage_long %>%
  group_by(metric_id, metric_plain, metric_plot, phase5, phase_label) %>%
  summarise(
    n = sum(is.finite(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    se = sd / sqrt(n),
    ci = qt(0.975, df = n - 1) * se,
    .groups = "drop"
  )

pre_values <- stage_long %>%
  filter(phase5 == "pre-ictal") %>%
  select(case_id, metric_id, pre_value = value)

change <- stage_long %>%
  filter(phase5 != "pre-ictal") %>%
  left_join(pre_values, by = c("case_id", "metric_id")) %>%
  filter(is.finite(value), is.finite(pre_value), value > 0, pre_value > 0) %>%
  mutate(
    dlog = log(value / pre_value),
    pct_change = 100 * (value / pre_value - 1),
    phase_label = factor(as.character(phase_label), levels = phase_labels[-1])
  )

change_stats <- change %>%
  group_by(metric_id, metric_plain, metric_plot, phase5, phase_label) %>%
  summarise(
    n = n(),
    mean_dlog = mean(dlog),
    median_dlog = median(dlog),
    ci_low = ci_one_sample(dlog)[1],
    ci_high = ci_one_sample(dlog)[2],
    mean_pct = 100 * (exp(mean_dlog) - 1),
    median_pct = 100 * (exp(median_dlog) - 1),
    p_wilcox = suppressWarnings(wilcox.test(dlog, mu = 0, exact = FALSE)$p.value),
    dz = ifelse(sd(dlog) > 0, mean(dlog) / sd(dlog), NA_real_),
    .groups = "drop"
  ) %>%
  group_by(metric_id) %>%
  mutate(
    p_holm = p.adjust(p_wilcox, method = "holm"),
    stars = stars_from_p(p_holm)
  ) %>%
  ungroup() %>%
  mutate(
    metric_plot = factor(as.character(metric_plot), levels = rev(metric_plot_levels)),
    phase_label = factor(as.character(phase_label), levels = phase_labels[-1])
  )

write_csv(stage_summary, file.path(base_dir, "Fig3_stage_summary.csv"))
write_csv(change_stats, file.path(base_dir, "Fig3_paired_logchange_statistics.csv"))

pre_heat <- metric_spec %>%
  transmute(
    metric_id,
    metric_plain,
    metric_plot = factor(metric_plot, levels = rev(metric_plot_levels)),
    phase5 = factor("pre-ictal", levels = phase_change_levels),
    phase_label = factor("Pre", levels = phase_change_labels),
    n = n_distinct(stage$case_id),
    mean_pct = 0,
    p_holm = NA_real_,
    stars = "",
    is_pre = TRUE
  )

heat_df <- change_stats %>%
  mutate(
    phase5 = factor(as.character(phase5), levels = phase_change_levels),
    phase_label = factor(as.character(phase_label), levels = phase_change_labels),
    is_pre = FALSE
  ) %>%
  bind_rows(pre_heat) %>%
  mutate(
    label = ifelse(is_pre, "0", paste0(ifelse(mean_pct >= 0, "+", ""), sprintf("%.1f", mean_pct), stars)),
    text_col = ifelse(abs(mean_pct) > 16, "white", "#202124")
  )

p_a <- ggplot(heat_df, aes(phase_label, metric_plot)) +
  geom_tile(
    data = heat_df %>% filter(is_pre),
    fill = "#F1F3F5", linewidth = 0.28, colour = "white", width = 0.92, height = 0.88
  ) +
  geom_tile(
    data = heat_df %>% filter(!is_pre),
    aes(fill = mean_pct),
    linewidth = 0.28, colour = "white", width = 0.92, height = 0.88
  ) +
  geom_text(aes(label = label, colour = text_col), size = FIG_GEOM_TEXT_SIZE, family = "Arial", fontface = "bold") +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = "#4778A8", mid = "white", high = "#C95C50",
    midpoint = 0, limits = c(-22, 22), oob = scales::squish,
    name = "Mean change\nfrom pre (%)"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete() +
  labs(
    title = "Paired component changes",
    x = NULL,
    y = NULL
  ) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

balance_wide <- change %>%
  filter(metric_id %in% c("eta", "resource_R", "lambda1")) %>%
  select(case_id, phase5, phase_label, metric_id, dlog) %>%
  pivot_wider(names_from = metric_id, values_from = dlog) %>%
  mutate(
    resource_contribution = 100 * resource_R,
    dominance_contribution = -100 * lambda1,
    eta_net = 100 * eta
  )

pre_balance_wide <- stage %>%
  distinct(case_id) %>%
  mutate(
    phase5 = factor("pre-ictal", levels = phase_change_levels),
    phase_label = factor("Pre", levels = phase_change_labels),
    eta = 0,
    resource_R = 0,
    lambda1 = 0,
    resource_contribution = 0,
    dominance_contribution = 0,
    eta_net = 0
  )

balance_wide_plot <- bind_rows(
  pre_balance_wide,
  balance_wide %>%
    mutate(
      phase5 = factor(as.character(phase5), levels = phase_change_levels),
      phase_label = factor(as.character(phase_label), levels = phase_change_labels)
    )
)

balance_long <- balance_wide_plot %>%
  select(case_id, phase5, phase_label, resource_contribution, dominance_contribution) %>%
  pivot_longer(
    cols = c(resource_contribution, dominance_contribution),
    names_to = "component",
    values_to = "contribution"
  ) %>%
  mutate(
    component = recode(
      component,
      resource_contribution = "Resource term (+dlog R)",
      dominance_contribution = paste0("Dominant term (-dlog ", lambda1_sym, ")")
    )
  )

balance_bar <- balance_long %>%
  group_by(phase5, phase_label, component) %>%
  summarise(mean = mean(contribution, na.rm = TRUE), .groups = "drop")

balance_fill_cols <- c("#4778A8", "#C95C50")
names(balance_fill_cols) <- c(
  paste0("Dominant term (-dlog ", lambda1_sym, ")"),
  "Resource term (+dlog R)"
)

net_summary <- balance_wide_plot %>%
  group_by(phase5, phase_label) %>%
  summarise(
    n = n(),
    mean = mean(eta_net, na.rm = TRUE),
    ci_low = ci_one_sample(eta_net)[1],
    ci_high = ci_one_sample(eta_net)[2],
    .groups = "drop"
  )

write_csv(
  balance_wide %>% arrange(case_id, phase5),
  file.path(base_dir, "Fig3_eta_logratio_component_balance.csv")
)

p_b <- ggplot() +
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.28, colour = "#7A828C") +
  geom_col(
    data = balance_bar,
    aes(phase_label, mean, fill = component),
    width = 0.58, colour = "white", linewidth = 0.20
  ) +
  geom_errorbar(
    data = net_summary,
    aes(phase_label, ymin = ci_low, ymax = ci_high),
    width = 0.12, linewidth = 0.35, colour = "#202124"
  ) +
  geom_line(
    data = net_summary,
    aes(phase_label, mean, group = 1),
    linewidth = 0.38, colour = "#202124"
  ) +
  geom_point(
    data = net_summary,
    aes(phase_label, mean),
    shape = 21, size = 2.35, stroke = 0.35, fill = "#202124", colour = "white"
  ) +
  scale_fill_manual(values = balance_fill_cols, name = NULL) +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = paste0("Ratio balance underlying ", eta_sym),
    x = NULL,
    y = paste0("Contribution to 100 ", times_sym, " dlog ", eta_sym)
  ) +
  theme(legend.position = "bottom")

net_eta <- balance_wide_plot %>%
  mutate(
    eta_net = 100 * eta,
    resource_excess = 100 * (resource_R - lambda1)
  )

net_stats <- net_eta %>%
  group_by(phase5, phase_label) %>%
  summarise(
    n = n(),
    mean = mean(resource_excess, na.rm = TRUE),
    median = median(resource_excess, na.rm = TRUE),
    ci_low = ci_one_sample(resource_excess)[1],
    ci_high = ci_one_sample(resource_excess)[2],
    p_wilcox = suppressWarnings(wilcox.test(resource_excess, mu = 0, exact = FALSE)$p.value),
    .groups = "drop"
  ) %>%
  mutate(
    p_holm = ifelse(as.character(phase5) == "pre-ictal", NA_real_, p.adjust(ifelse(as.character(phase5) == "pre-ictal", NA_real_, p_wilcox), method = "holm")),
    stars = ifelse(as.character(phase5) == "pre-ictal", "", stars_from_p(p_holm)),
    ypos = max(net_eta$resource_excess, na.rm = TRUE) + 0.12 * diff(range(net_eta$resource_excess, na.rm = TRUE))
  )

write_csv(net_stats, file.path(base_dir, "Fig3_resource_excess_statistics.csv"))

p_c_data <- net_eta %>% filter(as.character(phase5) != "pre-ictal")
p_c_pre <- tibble::tibble(
  phase_label = factor("Pre", levels = phase_change_labels),
  resource_excess = 0
)
p_c_stats <- net_stats %>% filter(as.character(phase5) != "pre-ictal")

p_c <- ggplot(p_c_data, aes(phase_label, resource_excess, fill = phase5)) +
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.28, colour = "#7A828C") +
  geom_violin(width = 0.82, trim = FALSE, alpha = 0.46, linewidth = 0.25, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.30, colour = "#202124", fill = "white", alpha = 0.88) +
  ggbeeswarm::geom_quasirandom(
    width = 0.12, shape = 21, size = 1.18, stroke = 0.10,
    colour = "white", alpha = 0.76
  ) +
  geom_errorbar(
    data = p_c_stats,
    aes(phase_label, ymin = ci_low, ymax = ci_high),
    inherit.aes = FALSE,
    width = 0.13, linewidth = 0.35, colour = "#202124"
  ) +
  geom_point(
    data = p_c_stats,
    aes(phase_label, mean),
    inherit.aes = FALSE,
    shape = 23, size = 2.15, fill = "#202124", colour = "white", stroke = 0.18
  ) +
  geom_point(
    data = p_c_pre,
    aes(phase_label, resource_excess),
    inherit.aes = FALSE,
    shape = 23, size = 2.35, fill = "#F1F3F5", colour = "#202124", stroke = 0.32
  ) +
  geom_text(
    data = p_c_stats,
    aes(phase_label, ypos, label = stars),
    inherit.aes = FALSE,
    family = "Arial", size = FIG_GEOM_TEXT_SIZE, fontface = "bold", colour = "#202124"
  ) +
  scale_fill_manual(values = phase_cols, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  coord_cartesian(
    ylim = c(
      min(net_eta$resource_excess, na.rm = TRUE) - 0.15 * diff(range(net_eta$resource_excess, na.rm = TRUE)),
      max(net_stats$ypos, na.rm = TRUE) + 0.06 * diff(range(net_eta$resource_excess, na.rm = TRUE))
    ),
    clip = "off"
  ) +
  labs(
    title = "Net resource excess over spectral dominance",
    x = NULL,
    y = paste0("100 ", times_sym, " dlog(R/", lambda1_sym, ")")
  ) +
  theme(
    axis.title.y = element_text(margin = margin(r = 1.5)),
    plot.margin = margin(3, 2, 3, 1)
  )

secondary_stats <- change_stats %>%
  filter(metric_id %in% c("lambda2_ratio", "spectral_PR")) %>%
  select(metric_id, phase_label, stars)

secondary_pre <- stage %>%
  distinct(case_id) %>%
  tidyr::crossing(
    metric_id = c("lambda2_ratio", "spectral_PR")
  ) %>%
  left_join(metric_spec %>% select(metric_id, metric_plot, metric_plain), by = "metric_id") %>%
  mutate(
    phase5 = factor("pre-ictal", levels = phase_change_levels),
    phase_label = factor("Pre", levels = phase_change_labels),
    dlog = 0,
    dlog100 = 0,
    metric_plot = factor(as.character(metric_plot), levels = c("u03bb2/u03bb1", "Spectral~PR"))
  )

secondary <- change %>%
  filter(metric_id %in% c("lambda2_ratio", "spectral_PR")) %>%
  mutate(
    dlog100 = 100 * dlog,
    phase5 = factor(as.character(phase5), levels = phase_change_levels),
    phase_label = factor(as.character(phase_label), levels = phase_change_labels),
    metric_plot = factor(as.character(metric_plot), levels = c("u03bb2/u03bb1", "Spectral~PR"))
  ) %>%
  bind_rows(secondary_pre)

secondary_plot <- secondary %>% filter(as.character(phase5) != "pre-ictal")
secondary_pre_marker <- secondary_pre %>%
  distinct(metric_id, metric_plot, phase_label, dlog100)

secondary_y <- secondary %>%
  group_by(metric_id) %>%
  summarise(
    ypos = max(dlog100, na.rm = TRUE) + 0.10 * diff(range(dlog100, na.rm = TRUE)),
    .groups = "drop"
  )

secondary_annot <- secondary_stats %>%
  left_join(secondary_y, by = "metric_id") %>%
  left_join(metric_spec %>% select(metric_id, metric_plot), by = "metric_id") %>%
  mutate(
    phase_label = factor(as.character(phase_label), levels = phase_change_labels),
    metric_plot = factor(metric_plot, levels = c("u03bb2/u03bb1", "Spectral~PR"))
  )

p_d <- ggplot(secondary_plot, aes(phase_label, dlog100, fill = phase5)) +
  geom_hline(yintercept = 0, linetype = "22", linewidth = 0.28, colour = "#7A828C") +
  geom_violin(width = 0.70, trim = FALSE, alpha = 0.42, linewidth = 0.22, colour = NA) +
  geom_boxplot(width = 0.13, outlier.shape = NA, linewidth = 0.26, colour = "#202124", fill = "white", alpha = 0.88) +
  ggbeeswarm::geom_quasirandom(
    width = 0.08, shape = 21, size = 0.82, stroke = 0.08,
    colour = "white", alpha = 0.72
  ) +
  stat_summary(
    fun = mean, geom = "point", shape = 23, size = 1.55,
    fill = "#202124", colour = "white", stroke = 0.15
  ) +
  geom_text(
    data = secondary_annot,
    aes(x = phase_label, y = ypos, label = stars),
    inherit.aes = FALSE,
    family = "Arial", size = FIG_GEOM_TEXT_SIZE, fontface = "bold", colour = "#202124"
  ) +
  geom_point(
    data = secondary_pre_marker,
    aes(phase_label, dlog100),
    inherit.aes = FALSE,
    shape = 23, size = 1.55, fill = "#F1F3F5", colour = "#202124", stroke = 0.25
  ) +
  facet_wrap(~ metric_plot, ncol = 1, scales = "free_y", labeller = label_parsed) +
  scale_fill_manual(values = phase_cols, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  labs(
    title = "Secondary spectral reorganization",
    x = NULL,
    y = paste0("Change from pre, 100 ", times_sym, " dlog")
  ) +
  theme(
    panel.spacing.y = unit(3.0, "mm"),
    strip.text = element_text(size = FIG_TEXT_PT, face = "bold"),
    axis.text.x = element_text(size = FIG_TEXT_PT),
    axis.title.y = element_text(margin = margin(r = 1.5)),
    plot.margin = margin(3, 2, 3, 1)
  )

fig <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(widths = c(1.08, 1.07), heights = c(1.00, 1.12), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = FIG_PANEL_PT, face = "bold", colour = "#202124"),
    plot.tag.position = c(0.055, 0.985),
    legend.position = "bottom"
  )

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 142, dpi = 600) {
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

save_pub_r(fig, out_base)

eta_stats <- change_stats %>%
  filter(metric_id == "eta") %>%
  mutate(summary = paste0(
    as.character(phase_label), ": ",
    sprintf("%+.2f%%", mean_pct),
    ", Holm P=", formatC(p_holm, format = "e", digits = 2)
  )) %>%
  pull(summary)

legend_text <- paste0(
  "Fig. 3 | ", eta_sym, " reflects proportional reorganization of propagation resources and spectral dominance.\n",
  "a, Paired stage changes in ", eta_sym, " and its decomposition metrics relative to the pre-ictal baseline, with the pre-ictal stage shown as the zero reference. Values show mean percentage changes from pre-ictal levels; asterisks denote paired Wilcoxon signed-rank tests with Holm correction within each metric. ",
  "b, Log-ratio decomposition of ", eta_sym, " = R/", lambda1_sym, ". Positive bars show the resource contribution (+dlog R), negative bars show the denominator contribution (-dlog ", lambda1_sym, "), and black points show the net ", eta_sym, " change; the pre-ictal stage is fixed at zero by definition. ",
  "c, Distribution of the net resource excess over spectral dominance, computed as 100 x dlog(R/", lambda1_sym, "), showing the individual seizure-level contribution that directly yields the proportional increase in ", eta_sym, ". ",
  "d, Stage-wise changes in secondary spectral organization, quantified by ", lambda2_sym, "/", lambda1_sym, " and spectral participation ratio. Source data are from n = 24 seizures. Tests are two-sided paired Wilcoxon signed-rank tests with Holm correction; *P < 0.05, **P < 0.01 and ***P < 0.001."
)
writeLines(legend_text, file.path(base_dir, "Fig3_legend_draft.txt"), useBytes = TRUE)

qa_lines <- c(
  "Core conclusion: eta increases because propagation resources R increase proportionally relative to lambda1; the effect is not caused by a decrease in lambda1.",
  "Interpretation guardrail: R and lambda1 both increase from pre-ictal baseline, while lambda2/lambda1 and spectral PR decrease, indicating stronger dominant spectral organization with reduced secondary-mode relative contribution.",
  "Figure archetype: quantitative mechanism grid with paired change heatmap, log-ratio decomposition, net resource-excess distribution and secondary spectral distributions.",
  "Backend: R only; ggplot2/patchwork/ggbeeswarm/ggrepel/svglite/cairo_pdf/ragg.",
  "Export: SVG/PDF/TIFF/PNG at double-column width.",
  "Statistics: paired log changes from pre-ictal baseline; two-sided paired Wilcoxon signed-rank tests with Holm correction within each metric.",
  paste("Eta paired changes:", paste(eta_stats, collapse = "; "))
)
writeLines(qa_lines, file.path(base_dir, "Fig3_QA_notes.txt"), useBytes = TRUE)

message("Fig. 3 outputs written to: ", base_dir)
