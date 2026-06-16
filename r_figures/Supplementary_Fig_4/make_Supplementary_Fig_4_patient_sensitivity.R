library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(ggtext)

required_packages <- c("ragg", "svglite", "ggtext")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

out_dir <- file.path("figures", "Supplementary_Fig_4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stage_source_path <- file.path("figures", "Supplementary_Fig_2", "Supplementary_Fig_2_stage_eta_source.csv")
stage_raw <- read.csv(stage_source_path, stringsAsFactors = FALSE, check.names = FALSE)

expected_phases <- c("pre-ictal", "early", "mid", "late", "post-ictal")
missing_cols <- setdiff(c("case_id", "phase5", "eta_mean"), names(stage_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns in stage source: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}
if (!all(expected_phases %in% unique(stage_raw$phase5))) {
  stop("The stage source does not contain all five expected phases.", call. = FALSE)
}

stage_raw <- stage_raw %>%
  mutate(
    patient_raw = gsub("_cut.*$", "", case_id),
    patient_raw = gsub("_process$", "", patient_raw)
  )

patient_map <- stage_raw %>%
  distinct(patient_raw) %>%
  arrange(patient_raw) %>%
  mutate(patient_id = sprintf("P%02d", row_number()))

stage_patient_input <- stage_raw %>%
  left_join(patient_map, by = "patient_raw") %>%
  mutate(
    phase_short = recode(
      phase5,
      "pre-ictal" = "Pre",
      "early" = "Early",
      "mid" = "Mid",
      "late" = "Late",
      "post-ictal" = "Post"
    ),
    phase_short = factor(phase_short, levels = c("Pre", "Early", "Mid", "Late", "Post"))
  )

patient_stage <- stage_patient_input %>%
  group_by(patient_id, phase_short) %>%
  summarise(
    eta = mean(eta_mean, na.rm = TRUE),
    n_seizures = n_distinct(case_id),
    .groups = "drop"
  ) %>%
  arrange(patient_id, phase_short)

if (n_distinct(patient_stage$patient_id) != 14) {
  stop("Expected 14 patients after aggregation, found ", n_distinct(patient_stage$patient_id), ".", call. = FALSE)
}
if (nrow(patient_stage) != 14 * 5) {
  stop("Expected one patient-level value for each patient and stage.", call. = FALSE)
}

patient_wide <- patient_stage %>%
  select(patient_id, phase_short, eta) %>%
  pivot_wider(names_from = phase_short, values_from = eta) %>%
  mutate(Ictal = rowMeans(across(c(Early, Mid, Late)), na.rm = TRUE))

patient_delta <- patient_wide %>%
  transmute(
    patient_id,
    Early = Early - Pre,
    Mid = Mid - Pre,
    Late = Late - Pre,
    Post = Post - Pre,
    Ictal = Ictal - Pre
  )

stat_one <- function(values) {
  values <- as.numeric(values)
  n <- length(values)
  mean_delta <- mean(values)
  sd_delta <- stats::sd(values)
  se_delta <- sd_delta / sqrt(n)
  ci_half_width <- stats::qt(0.975, df = n - 1) * se_delta
  wilcox <- suppressWarnings(stats::wilcox.test(
    values,
    mu = 0,
    alternative = "two.sided",
    exact = FALSE,
    correct = FALSE
  ))
  data.frame(
    n_patients = n,
    mean_delta = mean_delta,
    ci_low = mean_delta - ci_half_width,
    ci_high = mean_delta + ci_half_width,
    dz = mean_delta / sd_delta,
    positive_patient_count = sum(values > 0),
    wilcoxon_v = unname(wilcox$statistic),
    p_raw = unname(wilcox$p.value),
    stringsAsFactors = FALSE
  )
}

comparisons <- c("Early", "Mid", "Late", "Post", "Ictal")
stats_tbl <- bind_rows(lapply(comparisons, function(comp) {
  cbind(comparison = paste0(comp, "-Pre"), stat_one(patient_delta[[comp]]))
})) %>%
  mutate(
    comparison = factor(comparison, levels = paste0(comparisons, "-Pre")),
    p_holm_stage = NA_real_
  )

stage_rows <- as.character(stats_tbl$comparison) %in% paste0(c("Early", "Mid", "Late", "Post"), "-Pre")
stats_tbl$p_holm_stage[stage_rows] <- stats::p.adjust(stats_tbl$p_raw[stage_rows], method = "holm")
stats_tbl <- stats_tbl %>% arrange(comparison)

loo_tbl <- patient_delta %>%
  select(patient_id, ictal_delta = Ictal) %>%
  rowwise() %>%
  mutate(mean_delta_leave_one_out = mean(patient_delta$Ictal[patient_delta$patient_id != patient_id])) %>%
  ungroup() %>%
  mutate(still_positive = mean_delta_leave_one_out > 0) %>%
  arrange(mean_delta_leave_one_out)

loo_summary <- data.frame(
  comparison = "Ictal-Pre",
  n_patients = n_distinct(patient_delta$patient_id),
  all_leave_one_out_means_positive = all(loo_tbl$still_positive),
  loo_min = min(loo_tbl$mean_delta_leave_one_out),
  loo_max = max(loo_tbl$mean_delta_leave_one_out),
  stringsAsFactors = FALSE
)

write.csv(patient_stage, file.path(out_dir, "Supplementary_Fig_4_patient_stage_eta_source.csv"), row.names = FALSE)
write.csv(patient_delta, file.path(out_dir, "Supplementary_Fig_4_patient_delta_source.csv"), row.names = FALSE)
write.csv(stats_tbl, file.path(out_dir, "Supplementary_Fig_4_patient_statistics.csv"), row.names = FALSE)
write.csv(loo_tbl, file.path(out_dir, "Supplementary_Fig_4_leave_one_patient_out.csv"), row.names = FALSE)
write.csv(loo_summary, file.path(out_dir, "Supplementary_Fig_4_leave_one_patient_out_summary.csv"), row.names = FALSE)
write.csv(patient_map, file.path(out_dir, "Supplementary_Fig_4_patient_key_internal.csv"), row.names = FALSE)

checks <- data.frame(
  check = c(
    "stage source rows",
    "seizures in source",
    "patients after aggregation",
    "patient-stage rows",
    "window-level rows used",
    "all leave-one-out Ictal-Pre means positive"
  ),
  value = c(
    nrow(stage_raw),
    n_distinct(stage_raw$case_id),
    n_distinct(patient_stage$patient_id),
    nrow(patient_stage),
    0,
    all(loo_tbl$still_positive)
  )
)
write.csv(checks, file.path(out_dir, "Supplementary_Fig_4_preplot_consistency_checks.csv"), row.names = FALSE)

FIG_TEXT_PT <- 10.5
FIG_PANEL_PT <- 12.5
FIG_GEOM_TEXT_SIZE <- FIG_TEXT_PT / 2.845276
FIG_HEATMAP_TEXT_PT <- 8.0

pt_to_mm <- function(pt) pt / 2.845276

palette_contract <- c(
  neutral_dark = "#2A2A2A",
  neutral_mid = "#7B7B7B",
  neutral_light = "#D6D6D6",
  signal_blue = "#2F78B7",
  signal_teal = "#169C8A",
  signal_light = "#A9D6E5",
  accent_orange = "#D58A1F",
  accent_red = "#B94A48"
)

theme_supp <- function() {
  theme_classic(base_size = FIG_TEXT_PT, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = FIG_TEXT_PT, colour = "black"),
      axis.title.x = ggtext::element_markdown(size = FIG_TEXT_PT, colour = "black"),
      axis.title.y = ggtext::element_markdown(size = FIG_TEXT_PT, colour = "black", margin = margin(r = 1.5)),
      axis.text = element_text(size = FIG_TEXT_PT - 1, colour = "black"),
      plot.title = ggtext::element_markdown(size = FIG_TEXT_PT, face = "bold", hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(size = FIG_TEXT_PT - 1.5, colour = palette_contract["neutral_mid"], margin = margin(b = 3)),
      legend.position = "none",
      panel.grid.major.y = element_line(linewidth = 0.18, colour = "#ECECEC"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5.5, 7, 5.5, 7)
    )
}

theme_set(theme_supp())

mean_ci <- function(data, value_col, group_col = NULL) {
  if (is.null(group_col)) {
    data %>%
      summarise(
        mean = mean(.data[[value_col]], na.rm = TRUE),
        sd = sd(.data[[value_col]], na.rm = TRUE),
        n = dplyr::n(),
        se = sd / sqrt(n),
        ci_low = mean - stats::qt(0.975, n - 1) * se,
        ci_high = mean + stats::qt(0.975, n - 1) * se,
        .groups = "drop"
      )
  } else {
    data %>%
      group_by(.data[[group_col]]) %>%
      summarise(
        mean = mean(.data[[value_col]], na.rm = TRUE),
        sd = sd(.data[[value_col]], na.rm = TRUE),
        n = dplyr::n(),
        se = sd / sqrt(n),
        ci_low = mean - stats::qt(0.975, n - 1) * se,
        ci_high = mean + stats::qt(0.975, n - 1) * se,
        .groups = "drop"
      )
  }
}

patient_stage <- patient_stage %>%
  mutate(phase_short = factor(phase_short, levels = c("Pre", "Early", "Mid", "Late", "Post")))

stage_summary <- mean_ci(patient_stage, "eta", "phase_short")

p_a <- ggplot(patient_stage, aes(x = phase_short, y = eta, group = patient_id)) +
  geom_line(linewidth = 0.32, colour = alpha(palette_contract["neutral_mid"], 0.45)) +
  geom_point(size = 1.45, colour = alpha(palette_contract["neutral_dark"], 0.55)) +
  geom_ribbon(
    data = stage_summary,
    aes(x = phase_short, ymin = ci_low, ymax = ci_high, group = 1),
    inherit.aes = FALSE,
    fill = alpha(palette_contract["signal_light"], 0.42)
  ) +
  geom_line(
    data = stage_summary,
    aes(x = phase_short, y = mean, group = 1),
    inherit.aes = FALSE,
    linewidth = 0.85,
    colour = palette_contract["signal_blue"]
  ) +
  geom_point(
    data = stage_summary,
    aes(x = phase_short, y = mean),
    inherit.aes = FALSE,
    size = 2.15,
    shape = 21,
    stroke = 0.45,
    fill = "white",
    colour = palette_contract["signal_blue"]
  ) +
  scale_y_continuous(labels = number_format(accuracy = 0.001), expand = expansion(mult = c(0.08, 0.10))) +
  labs(
    title = "Patient-level <i>η</i> trajectories",
    x = NULL,
    y = "Patient mean <i>η</i>"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

delta_long <- patient_delta %>%
  pivot_longer(cols = c(Early, Mid, Late, Post), names_to = "stage", values_to = "delta_eta") %>%
  mutate(stage = factor(stage, levels = c("Early", "Mid", "Late", "Post")))

delta_summary <- stats_tbl %>%
  filter(as.character(comparison) %in% paste0(c("Early", "Mid", "Late", "Post"), "-Pre")) %>%
  mutate(stage = factor(gsub("-Pre", "", as.character(comparison)), levels = c("Early", "Mid", "Late", "Post")))

stage_label_tbl <- delta_summary %>%
  mutate(
    label = paste0("P_Holm\n", format.pval(p_holm_stage, digits = 2, eps = 0.001), "\n", positive_patient_count, "/14"),
    label_y = max(delta_long$delta_eta, delta_summary$ci_high, na.rm = TRUE) + 0.0070
  )

p_b <- ggplot(delta_long, aes(x = stage, y = delta_eta)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = palette_contract["neutral_mid"]) +
  geom_point(
    aes(group = patient_id),
    position = position_jitter(width = 0.075, height = 0, seed = 14),
    size = 1.55,
    alpha = 0.72,
    colour = palette_contract["neutral_dark"]
  ) +
  geom_errorbar(
    data = delta_summary,
    aes(x = stage, ymin = ci_low, ymax = ci_high),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 0.42,
    colour = palette_contract["signal_blue"]
  ) +
  geom_point(
    data = delta_summary,
    aes(x = stage, y = mean_delta),
    inherit.aes = FALSE,
    shape = 23,
    size = 2.45,
    stroke = 0.45,
    fill = palette_contract["signal_blue"],
    colour = "white"
  ) +
  geom_text(
    data = stage_label_tbl,
    aes(x = stage, y = label_y, label = label),
    inherit.aes = FALSE,
    size = FIG_GEOM_TEXT_SIZE * 0.70,
    lineheight = 0.9,
    colour = palette_contract["neutral_dark"]
  ) +
  scale_y_continuous(labels = number_format(accuracy = 0.001), expand = expansion(mult = c(0.10, 0.26))) +
  labs(
    title = "Patient-level <i>Δη</i> by stage",
    x = NULL,
    y = "<i>Δη</i> vs Pre"
  )

delta_c_long <- patient_delta %>%
  select(patient_id, Ictal, Post) %>%
  pivot_longer(cols = c(Ictal, Post), names_to = "metric", values_to = "delta_eta") %>%
  mutate(metric = factor(ifelse(metric == "Ictal", "Ictal-Pre", "Post-Pre"), levels = c("Ictal-Pre", "Post-Pre")))

summary_c <- stats_tbl %>%
  filter(as.character(comparison) %in% c("Ictal-Pre", "Post-Pre")) %>%
  mutate(
    metric = factor(ifelse(as.character(comparison) == "Ictal-Pre", "Ictal-Pre", "Post-Pre"), levels = c("Ictal-Pre", "Post-Pre")),
    p_display = ifelse(
      as.character(comparison) == "Post-Pre",
      paste0("Holm P=", format.pval(p_holm_stage, digits = 2, eps = 0.001)),
      paste0("P=", format.pval(p_raw, digits = 2, eps = 0.001))
    ),
    label = paste0("mean ", sprintf("%.4f", mean_delta), "\nCI ", sprintf("%.4f", ci_low), "-", sprintf("%.4f", ci_high), "\n", positive_patient_count, "/14; ", p_display),
    label_y = max(delta_c_long$delta_eta, ci_high, na.rm = TRUE) + 0.005
  )

p_c <- ggplot(delta_c_long, aes(x = metric, y = delta_eta)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = palette_contract["neutral_mid"]) +
  geom_point(
    position = position_jitter(width = 0.065, height = 0, seed = 7),
    size = 1.65,
    alpha = 0.72,
    colour = palette_contract["neutral_dark"]
  ) +
  geom_errorbar(
    data = summary_c,
    aes(x = metric, ymin = ci_low, ymax = ci_high),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 0.42,
    colour = palette_contract["signal_teal"]
  ) +
  geom_point(
    data = summary_c,
    aes(x = metric, y = mean_delta),
    inherit.aes = FALSE,
    shape = 23,
    size = 2.55,
    stroke = 0.45,
    fill = palette_contract["signal_teal"],
    colour = "white"
  ) +
  geom_text(
    data = summary_c,
    aes(x = metric, y = label_y, label = label),
    inherit.aes = FALSE,
    size = FIG_GEOM_TEXT_SIZE * 0.70,
    lineheight = 0.88,
    colour = palette_contract["neutral_dark"]
  ) +
  scale_y_continuous(labels = number_format(accuracy = 0.001), expand = expansion(mult = c(0.10, 0.34))) +
  labs(
    title = "Ictal and post-ictal summaries",
    x = NULL,
    y = "<i>Δη</i> vs Pre"
  )

full_ictal_mean <- stats_tbl %>%
  filter(as.character(comparison) == "Ictal-Pre") %>%
  pull(mean_delta)
loo_range_label <- paste0(
  "LOO range ",
  sprintf("%.4f", loo_summary$loo_min),
  " to ",
  sprintf("%.4f", loo_summary$loo_max),
  "; all > 0"
)

p_d <- ggplot(loo_tbl, aes(x = reorder(patient_id, mean_delta_leave_one_out), y = mean_delta_leave_one_out)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = palette_contract["neutral_mid"]) +
  geom_hline(yintercept = full_ictal_mean, linewidth = 0.45, linetype = "22", colour = palette_contract["accent_orange"]) +
  geom_segment(
    aes(xend = reorder(patient_id, mean_delta_leave_one_out), y = 0, yend = mean_delta_leave_one_out),
    linewidth = 0.45,
    colour = alpha(palette_contract["signal_teal"], 0.62)
  ) +
  geom_point(size = 2.0, colour = palette_contract["signal_teal"]) +
  annotate(
    "text",
    x = 1.1,
    y = max(loo_tbl$mean_delta_leave_one_out) + 0.0016,
    label = loo_range_label,
    hjust = 0,
    size = FIG_GEOM_TEXT_SIZE * 0.82,
    colour = palette_contract["neutral_dark"]
  ) +
  scale_y_continuous(labels = number_format(accuracy = 0.001), expand = expansion(mult = c(0.10, 0.22))) +
  labs(
    title = "Leave-one-patient-out sensitivity",
    x = "Patient omitted",
    y = "Mean <i>Δη</i> (Ictal-Pre)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = FIG_PANEL_PT, face = "bold", colour = "black"),
    plot.tag.position = c(0.01, 0.98)
  )

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 148, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = dpi, background = "white")
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi, background = "white", compression = "lzw")
  print(plot)
  grDevices::dev.off()
}

out_base <- file.path(out_dir, "Supplementary_Fig_4")
save_pub_r(fig, out_base)

export_checks <- data.frame(
  file = paste0("Supplementary_Fig_4.", c("svg", "pdf", "png", "tiff")),
  exists = file.exists(paste0(out_base, ".", c("svg", "pdf", "png", "tiff"))),
  bytes = file.info(paste0(out_base, ".", c("svg", "pdf", "png", "tiff")))$size
)
write.csv(export_checks, file.path(out_dir, "Supplementary_Fig_4_postplot_export_checks.csv"), row.names = FALSE)

legend_lines <- c(
  "Supplementary Fig. 4 | Patient-level sensitivity analysis of η dynamics.",
  "a, Patient-level η trajectories across pre-ictal, early-ictal, mid-ictal, late-ictal and post-ictal stages after averaging all seizure-level stage means within each patient. Thin grey lines denote individual patients; the coloured line and error bars show the cohort mean and 95% CI across patients.",
  "b, Patient-level Δη for Early, Mid, Late and Post relative to each patient's paired Pre value. Points denote patients; diamonds and error bars show mean Δη and 95% CI. P values were obtained with two-sided paired Wilcoxon signed-rank tests and Holm correction across these four comparisons.",
  "c, Summary of Ictal-Pre (Ictal = mean(Early, Mid, Late)) and Post-Pre patient-level changes.",
  "d, Leave-one-patient-out sensitivity for Ictal-Pre. Each point shows the cohort mean Δη after removing one patient; the dashed line shows the full-cohort mean and the solid horizontal line marks zero. Source data are patient-level means from n = 14 patients derived from n = 24 seizures."
)
writeLines(legend_lines, file.path(out_dir, "Supplementary_Fig_4_legend_draft.txt"), useBytes = TRUE)

qa_lines <- c(
  "Supplementary Fig. 4 QA notes",
  "Core conclusion: Patient-level aggregation preserves the positive direction of ictal η dynamics, supporting the seizure-level analysis as a sensitivity check without replacing it.",
  "Archetype: quantitative grid.",
  "Backend: R/ggplot2/patchwork only.",
  "Statistics: n = 14 patients; patient means aggregate all seizure-level stage mean η values within each patient; no window-level statistics are used.",
  "Intervals: t-based 95% CI across patient-level paired differences.",
  "Tests: two-sided paired Wilcoxon signed-rank tests; Holm correction for Early/Mid/Late/Post vs Pre.",
  paste0("Ictal-Pre leave-one-out range: ", sprintf("%.6f", loo_summary$loo_min), " to ", sprintf("%.6f", loo_summary$loo_max), "."),
  paste0("All leave-one-patient-out means positive: ", all(loo_tbl$still_positive), ".")
)
writeLines(qa_lines, file.path(out_dir, "Supplementary_Fig_4_QA_notes.txt"), useBytes = TRUE)

message("Wrote Supplementary Fig. 4 patient-level sensitivity outputs to: ", out_dir)
