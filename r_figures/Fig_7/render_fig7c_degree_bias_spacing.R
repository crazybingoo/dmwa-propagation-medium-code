#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(".", winslash = "/", mustWork = TRUE)
}
source_script <- file.path(script_dir, "make_fig7_eta_ablation_final.R")
source_lines <- readLines(source_script, warn = FALSE, encoding = "UTF-8")

replacements <- c(
  'c("DMW-HLG", "Coverage", "Degree bias")[match(as.character(model), model_levels)],' =
    'c("DMW-HLG", "Coverage", "Degree-bias")[match(as.character(model), model_levels)],',
  'y = c(0.34, 0.20, 0.06)' =
    'y = c(0.25, 0.16, 0.07)',
  'base_file <- file.path(out_dir, "Fig7_eta_ablation_components_method_highlight_softer_delta_label_clear_dmwhlg_label_titlebold_font10p5_label12p5_microadjust_eta_entity_label_higher")' =
    'base_file <- file.path(out_dir, "Fig7_eta_ablation_components_fig7c_degree_bias_spacing")'
)

for (old in names(replacements)) {
  matched <- trimws(source_lines) == old
  if (sum(matched) != 1L) {
    stop("Expected one exact match for: ", old, "; found ", sum(matched), ".")
  }
  source_lines[matched] <- unname(replacements[[old]])
}

export_line <- grep("^save_pub\\(fig, base_file\\)$", source_lines)
if (length(export_line) != 1L) {
  stop("Could not identify the unique Fig. 7 export line.")
}

write.csv <- function(...) invisible(NULL)
invisible(Sys.setlocale("LC_CTYPE", "Chinese_China.utf8"))
eval(parse(text = source_lines[seq_len(export_line)]), envir = globalenv())
