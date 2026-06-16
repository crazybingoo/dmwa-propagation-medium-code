#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(svglite)
  library(ragg)
})

FIG_TEXT_PT <- 6
FIG_PANEL_PT <- 8
FIG_GEOM_TEXT_SIZE <- FIG_TEXT_PT / 2.845276


phase_order <- c("pre-ictal", "early", "mid", "late", "post-ictal")
phase_label <- c(
  "pre-ictal" = "Pre",
  "early" = "Early",
  "mid" = "Mid",
  "late" = "Late",
  "post-ictal" = "Post"
)

region_order <- c("SOZ", "PZ", "NIZ")
macro_order <- c("SOZ_only", "PZ_only", "NIZ_only", "SOZ_PZ", "SOZ_NIZ", "PZ_NIZ", "SOZ_PZ_NIZ")
role_order <- c("source-like", "balanced", "sink-like")

input_dir <- file.path("data", "Fig_5_size_adjusted")
out_dir_candidates <- Sys.glob(file.path("example_project", "*0514-", "nature_fig", "Fig_5"))
if (length(out_dir_candidates) < 1) {
  stop("Cannot locate Fig_5 directory under example_project/*0514-/nature_fig")
}
out_dir <- normalizePath(out_dir_candidates[1], winslash = "/", mustWork = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

region_case_csv <- file.path(input_dir, "ALL_CASES_exp4_stage_region_contribution_size_adjusted.csv")
macro_case_csv <- file.path(input_dir, "ALL_CASES_exp4_stage_macro_role_contribution.csv")
region_stats_csv <- file.path(input_dir, "GROUP_exp4_region_phase_stats_size_adjusted.csv")
macro_stats_csv <- file.path(input_dir, "GROUP_exp4_macro_role_phase_stats.csv")

required_files <- c(region_case_csv, macro_case_csv, region_stats_csv, macro_stats_csv)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

region_case <- read.csv(region_case_csv, check.names = FALSE)
macro_case <- read.csv(macro_case_csv, check.names = FALSE)
region_stats <- read.csv(region_stats_csv, check.names = FALSE)
macro_stats <- read.csv(macro_stats_csv, check.names = FALSE)

as_phase <- function(x) factor(x, levels = phase_order, labels = unname(phase_label[phase_order]))

mean_sem <- function(x) {
  x <- x[is.finite(x)]
  tibble(
    n = length(x),
    mean = ifelse(length(x) == 0, NA_real_, mean(x)),
    sd = ifelse(length(x) <= 1, NA_real_, sd(x)),
    sem = ifelse(length(x) <= 1, NA_real_, sd(x) / sqrt(length(x)))
  )
}

format_num <- function(x, digits = 5) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

# Figure contract:
# Core conclusion: eta burden is distributed beyond the seizure onset zone; after
# regional size adjustment, PZ/NIZ retain comparable per-node burdens, while total
# hyperedge-level burden is dominated by PZ-NIZ combinations and balanced roles.
# Archetype: quantitative grid with three linked line-summary panels.
# Backend: R only for plot generation, export and visual QA.

validate_levels <- function(df, var, expected, label) {
  observed <- sort(unique(as.character(df[[var]])))
  missing <- setdiff(expected, observed)
  extra <- setdiff(observed, expected)
  tibble(
    check = label,
    status = ifelse(length(missing) == 0 && length(extra) == 0, "PASS", "WARN"),
    detail = paste0(
      "observed=", paste(observed, collapse = "|"),
      "; missing=", ifelse(length(missing) == 0, "none", paste(missing, collapse = "|")),
      "; extra=", ifelse(length(extra) == 0, "none", paste(extra, collapse = "|"))
    )
  )
}

pre_checks <- bind_rows(
  validate_levels(region_case, "phase5", phase_order, "region phase levels"),
  validate_levels(region_case, "region", region_order, "region labels"),
  validate_levels(macro_case, "phase5", phase_order, "macro phase levels"),
  validate_levels(macro_case, "macro_state", macro_order, "macro-state labels"),
  validate_levels(macro_case, "role", role_order, "role labels"),
  tibble(
    check = "case count",
    status = ifelse(length(unique(region_case$case_id)) == 24 && length(unique(macro_case$case_id)) == 24, "PASS", "WARN"),
    detail = paste0(
      "region cases=", length(unique(region_case$case_id)),
      "; macro-role cases=", length(unique(macro_case$case_id))
    )
  ),
  tibble(
    check = "panel A metric",
    status = ifelse("eta_abs_burden_per_node_mean" %in% names(region_case), "PASS", "FAIL"),
    detail = "Panel A uses eta_abs_burden_per_node_mean to match size-adjusted screenshot."
  ),
  tibble(
    check = "panel B/C metric",
    status = ifelse("eta_abs_burden_mean" %in% names(macro_case), "PASS", "FAIL"),
    detail = "Panels B/C use eta_abs_burden_mean, matching total absolute-proxy burden in the screenshot."
  )
)

if (any(pre_checks$status == "FAIL")) {
  write.csv(pre_checks, file.path(out_dir, "Fig5_preplot_consistency_checks.csv"), row.names = FALSE)
  stop("Pre-plot consistency checks failed. See Fig5_preplot_consistency_checks.csv.")
}

panel_a <- region_case %>%
  filter(region %in% region_order, phase5 %in% phase_order) %>%
  mutate(
    region = factor(region, levels = region_order),
    phase5 = factor(phase5, levels = phase_order),
    phase = as_phase(as.character(phase5))
  ) %>%
  group_by(region, phase5, phase) %>%
  summarise(mean_sem(eta_abs_burden_per_node_mean), .groups = "drop")

panel_b_case <- macro_case %>%
  filter(macro_state %in% macro_order, role %in% role_order, phase5 %in% phase_order) %>%
  mutate(
    macro_state = factor(macro_state, levels = macro_order),
    role = factor(role, levels = role_order),
    phase5 = factor(phase5, levels = phase_order)
  ) %>%
  group_by(case_id, macro_state, phase5) %>%
  summarise(value = sum(eta_abs_burden_mean, na.rm = TRUE), .groups = "drop")

panel_b <- panel_b_case %>%
  mutate(phase = as_phase(as.character(phase5))) %>%
  group_by(macro_state, phase5, phase) %>%
  summarise(mean_sem(value), .groups = "drop")

panel_c_case <- macro_case %>%
  filter(macro_state %in% macro_order, role %in% role_order, phase5 %in% phase_order) %>%
  mutate(
    macro_state = factor(macro_state, levels = macro_order),
    role = factor(role, levels = role_order),
    phase5 = factor(phase5, levels = phase_order)
  ) %>%
  group_by(case_id, role, phase5) %>%
  summarise(value = sum(eta_abs_burden_mean, na.rm = TRUE), .groups = "drop")

panel_c <- panel_c_case %>%
  mutate(phase = as_phase(as.character(phase5))) %>%
  group_by(role, phase5, phase) %>%
  summarise(mean_sem(value), .groups = "drop")

source_data <- bind_rows(
  panel_a %>% transmute(panel = "a", group = as.character(region), phase = as.character(phase), n, mean, sd, sem),
  panel_b %>% transmute(panel = "b", group = as.character(macro_state), phase = as.character(phase), n, mean, sd, sem),
  panel_c %>% transmute(panel = "c", group = as.character(role), phase = as.character(phase), n, mean, sd, sem)
)
write.csv(source_data, file.path(out_dir, "Fig5_source_data.csv"), row.names = FALSE)

expected_summary <- tribble(
  ~panel, ~expected_pattern,
  "a", "SOZ/PZ/NIZ start close in pre-ictal stage after per-node adjustment; SOZ drops most strongly in early ictal stage; PZ and NIZ remain mutually close.",
  "b", "PZ_NIZ is the largest macro-state combination across stages; SOZ_only and SOZ_PZ_NIZ are low; mixed-region categories contribute more than single SOZ-only burden.",
  "c", "Balanced roles dominate total eta burden across all stages; source-like and sink-like burdens remain near zero."
)

pattern_checks <- bind_rows(
  panel_a %>%
    select(group = region, phase, mean) %>%
    pivot_wider(names_from = group, values_from = mean) %>%
    arrange(factor(phase, levels = unname(phase_label[phase_order]))) %>%
    summarise(
      check = "screenshot pattern panel A",
      status = ifelse(
        abs(SOZ[1] - PZ[1]) < 0.003 &&
          abs(SOZ[1] - NIZ[1]) < 0.003 &&
          SOZ[2] < SOZ[1] &&
          abs(PZ[5] - NIZ[5]) < 0.002,
        "PASS", "WARN"
      ),
      detail = paste0(
        "Pre SOZ/PZ/NIZ=", format_num(SOZ[1]), "/", format_num(PZ[1]), "/", format_num(NIZ[1]),
        "; Early SOZ=", format_num(SOZ[2]),
        "; Post PZ/NIZ=", format_num(PZ[5]), "/", format_num(NIZ[5])
      )
    ),
  panel_b %>%
    select(group = macro_state, phase, mean) %>%
    group_by(phase) %>%
    summarise(top_group = as.character(group[which.max(mean)]), top_mean = max(mean), .groups = "drop") %>%
    summarise(
      check = "screenshot pattern panel B",
      status = ifelse(all(top_group == "PZ_NIZ"), "PASS", "WARN"),
      detail = paste0("top macro-state by phase=", paste(paste(phase, top_group, format_num(top_mean), sep = ":"), collapse = "; "))
    ),
  panel_c %>%
    select(group = role, phase, mean) %>%
    group_by(phase) %>%
    summarise(top_group = as.character(group[which.max(mean)]), top_mean = max(mean), .groups = "drop") %>%
    summarise(
      check = "screenshot pattern panel C",
      status = ifelse(all(top_group == "balanced"), "PASS", "WARN"),
      detail = paste0("top role by phase=", paste(paste(phase, top_group, format_num(top_mean), sep = ":"), collapse = "; "))
    )
)

consistency_report <- bind_rows(pre_checks, pattern_checks)
write.csv(consistency_report, file.path(out_dir, "Fig5_preplot_consistency_checks.csv"), row.names = FALSE)

if (any(consistency_report$status == "WARN")) {
  warning("Some consistency checks emitted WARN; inspect Fig5_preplot_consistency_checks.csv.")
}

region_palette <- c(SOZ = "#D55E00", PZ = "#0072B2", NIZ = "#009E73")
macro_palette <- c(
  SOZ_only = "#202020",
  PZ_only = "#686868",
  NIZ_only = "#9A9A9A",
  SOZ_PZ = "#D6603D",
  SOZ_NIZ = "#4C78A8",
  PZ_NIZ = "#4E9A4B",
  SOZ_PZ_NIZ = "#8056B3"
)
role_palette <- c(`source-like` = "#CC79A7", balanced = "#5A5A5A", `sink-like` = "#56B4E9")

theme_nature <- function(base_size = FIG_TEXT_PT, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#1E1E1E"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#1E1E1E"),
      axis.title = element_text(size = FIG_TEXT_PT),
      axis.text = element_text(size = FIG_TEXT_PT, colour = "#252525"),
      plot.title = element_text(size = FIG_TEXT_PT, face = "bold", hjust = 0),
      legend.title = element_blank(),
      legend.text = element_text(size = FIG_TEXT_PT),
      legend.key.height = unit(3.1, "mm"),
      legend.key.width = unit(5.0, "mm"),
      legend.spacing.y = unit(0.3, "mm"),
      panel.grid.major.y = element_line(linewidth = 0.18, colour = "#E8EDF2"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
}

line_panel <- function(df, x_group, y_label, palette, title, y_limits = NULL, y_breaks = waiver(), legend_position = "inside") {
  p <- ggplot(df, aes(x = phase, y = mean, group = {{ x_group }}, colour = {{ x_group }})) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.08, linewidth = 0.38, alpha = 0.95) +
    geom_line(linewidth = 0.58) +
    geom_point(shape = 21, size = 1.8, stroke = 0.35, fill = "white") +
    scale_colour_manual(values = palette, labels = function(x) gsub("_", "-", x, fixed = TRUE)) +
    scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = expansion(mult = c(0.04, 0.08))) +
    labs(x = NULL, y = y_label, title = title) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
      legend.position = legend_position
    )
  if (identical(legend_position, "inside")) {
    p <- p + theme(
      legend.position.inside = c(0.78, 0.82),
      legend.background = element_rect(fill = scales::alpha("white", 0.76), colour = NA),
      legend.margin = margin(0, 0, 0, 0)
    )
  }
  p
}

p_a <- line_panel(
  panel_a,
  region,
  expression(eta~"burden per node"),
  region_palette,
  "Regional eta burden\n(size-adjusted)",
  y_limits = c(0.009, 0.0192),
  y_breaks = seq(0.009, 0.019, 0.002)
)

p_b <- line_panel(
  panel_b,
  macro_state,
  expression(eta~"burden, absolute proxy"),
  macro_palette,
  "Macro-state combination",
  y_limits = c(0, 0.115),
  y_breaks = seq(0, 0.10, 0.025),
  legend_position = "inside"
) +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(linewidth = 0.65, size = 1.7))) +
  theme(
    legend.position.inside = c(0.99, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.86), colour = NA),
    legend.margin = margin(0, 0, 0, 0)
  )

p_c <- line_panel(
  panel_c,
  role,
  expression(eta~"burden, absolute proxy"),
  role_palette,
  "Source/sink role",
  y_limits = c(0, 0.25),
  y_breaks = seq(0, 0.25, 0.05)
)

fig <- (p_a | p_b | p_c) +
  plot_layout(widths = c(1.06, 1.16, 0.96), guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = FIG_PANEL_PT, face = "bold", family = "Arial"))

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 68, dpi = 600) {
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

out_base <- file.path(out_dir, "Fig5_eta_burden_distribution_font6_label8")
save_pub_r(fig, out_base)

post_checks <- tibble(
  check = c("svg export", "pdf export", "tiff export", "png export", "source data", "preplot consistency report"),
  path = c(
    paste0(out_base, ".svg"),
    paste0(out_base, ".pdf"),
    paste0(out_base, ".tiff"),
    paste0(out_base, ".png"),
    file.path(out_dir, "Fig5_source_data.csv"),
    file.path(out_dir, "Fig5_preplot_consistency_checks.csv")
  )
) %>%
  mutate(
    exists = file.exists(path),
    bytes = ifelse(exists, file.info(path)$size, NA_real_),
    status = ifelse(exists & bytes > 0, "PASS", "FAIL")
  )
write.csv(post_checks, file.path(out_dir, "Fig5_postplot_export_checks.csv"), row.names = FALSE)

panel_a_values <- panel_a %>%
  select(group = region, phase, mean, sem) %>%
  mutate(value = paste0(format_num(mean, 4), " +/- ", format_num(sem, 4))) %>%
  select(group, phase, value) %>%
  pivot_wider(names_from = phase, values_from = value)

panel_b_top <- panel_b %>%
  group_by(phase) %>%
  slice_max(order_by = mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(phase, top_macro = as.character(macro_state), mean = format_num(mean, 4))

panel_c_values <- panel_c %>%
  select(group = role, phase, mean, sem) %>%
  mutate(value = paste0(format_num(mean, 4), " +/- ", format_num(sem, 4))) %>%
  select(group, phase, value) %>%
  pivot_wider(names_from = phase, values_from = value)

qa_lines <- c(
  "Fig. 5 QA notes",
  "",
  "Figure contract:",
  "Core conclusion: eta burden is distributed beyond the SOZ; after regional size adjustment, PZ and NIZ retain comparable per-node burdens, while total burden is dominated by PZ-NIZ macro-state combinations and balanced source/sink roles.",
  "Archetype: quantitative grid.",
  "Backend: R only; ggplot2 + patchwork + svglite/cairo_pdf/ragg.",
  "Export: double-column SVG/PDF/TIFF/PNG.",
  "",
  "Data and screenshot consistency:",
  paste0("- Panel A uses eta_abs_burden_per_node_mean from the size-adjusted regional table; this matches the screenshot label 'Regional eta burden (size-adjusted)' and y-axis 'eta burden per node'."),
  paste0("- Panels B and C use eta_abs_burden_mean from the macro-role table; this matches the screenshot note that only panel A is size adjusted, whereas panels B/C show total absolute-proxy burden."),
  paste0("- n = ", length(unique(region_case$case_id)), " seizures for regional summaries and ", length(unique(macro_case$case_id)), " seizures for macro-role summaries."),
  paste0("- Pre-plot consistency check status: ", paste(unique(consistency_report$status), collapse = ", "), "."),
  "",
  "Panel A mean +/- SEM:",
  paste(capture.output(print(panel_a_values, n = Inf)), collapse = "\n"),
  "",
  "Panel B largest macro-state by phase:",
  paste(capture.output(print(panel_b_top, n = Inf)), collapse = "\n"),
  "",
  "Panel C mean +/- SEM:",
  paste(capture.output(print(panel_c_values, n = Inf)), collapse = "\n"),
  "",
  "Interpretation boundary:",
  "Do not describe panels B/C as size-adjusted per-node effects. Only panel A supports size-adjusted regional per-node claims; panels B/C support total burden claims."
)
writeLines(qa_lines, file.path(out_dir, "Fig5_QA_notes.txt"))

legend_lines <- c(
  "Fig. 5 | Regional and hyperedge-level distribution of eta burden.",
  "a, Size-adjusted regional eta burden across seizure stages for hyperedges involving the seizure onset zone (SOZ), propagation zone (PZ) and non-involved zone (NIZ). Values are normalized per regional node. b, Total eta burden grouped by macro-state combinations of SOZ, PZ and NIZ involvement. c, Total eta burden grouped by source-like, balanced and sink-like hyperedge roles. Points denote seizure-level cohort means and error bars denote s.e.m.; source data are from n = 24 seizures. Panel a uses per-node size-adjusted burden, whereas b and c show total absolute-proxy burden after summing across roles or macro-state combinations, respectively."
)
writeLines(legend_lines, file.path(out_dir, "Fig5_legend_draft.txt"))

result_lines <- c(
  "η burden is distributed across regional and hyperedge-level contributors",
  "",
  "We next asked which regional and hyperedge categories carried the DMWA-derived η burden. Because SOZ, PZ and NIZ categories differed in their available node counts, regional burden was first examined after per-node size adjustment. Under this normalization, SOZ, PZ and NIZ showed comparable pre-ictal unit burdens, but SOZ burden decreased more prominently after seizure onset, whereas PZ and NIZ remained closer to one another across ictal and post-ictal stages (Fig. 5a). This indicates that the apparent regional contribution of SOZ is sensitive to regional sampling size and should not be interpreted as a stable SOZ-dominant effect.",
  "",
  "At the hyperedge-combination level, total η burden was not concentrated in single-region categories. PZ-NIZ combinations carried the largest absolute-proxy burden across all stages, whereas SOZ-only and three-region SOZ-PZ-NIZ combinations remained comparatively small (Fig. 5b). When the same burden was grouped by source/sink role, balanced hyperedges dominated the total burden, while source-like and sink-like categories contributed only weakly (Fig. 5c). These results support a distributed propagation-medium interpretation: η burden is carried primarily by cross-regional and relatively balanced higher-order interactions, rather than by a single onset-zone compartment alone.",
  "",
  "Boundary: Panel a supports size-adjusted per-node regional comparisons; panels b and c support total burden comparisons and should not be described as per-node normalized effects."
)
writeLines(result_lines, file.path(out_dir, "Fig5_results_description_draft.txt"))

message("Fig. 5 complete: ", out_dir)
