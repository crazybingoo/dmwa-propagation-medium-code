#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(".", winslash = "/", mustWork = TRUE)
}
source_script <- file.path(script_dir, "make_fig4_topology_controls_font10p5_label12p5_no_c_legend_direct_labels.R")
source_lines <- readLines(source_script, warn = FALSE, encoding = "UTF-8")

replacements <- c(
  'out_base <- file.path(base_dir, "Fig4_topology_controls_font10p5_label12p5_no_c_legend_direct_labels")' =
    'out_base <- file.path(base_dir, "Fig4_topology_controls_control_label_nouns")',
  'model_labels <- c("Original", "Density matched", "Weight shuffled", "Topology shuffled")' =
    'model_labels <- c("Empirical", "density matching", "weight shuffling", "topology shuffling")'
)

for (old in names(replacements)) {
  matched <- trimws(source_lines) == old
  if (sum(matched) != 1L) {
    stop("Expected one exact match for: ", old, "; found ", sum(matched), ".")
  }
  source_lines[matched] <- unname(replacements[[old]])
}

export_line <- grep("^save_pub_r\\(fig, out_base\\)$", source_lines)
if (length(export_line) != 1L) {
  stop("Could not identify the unique Fig. 4 export line.")
}

write_csv <- function(...) invisible(NULL)
eval(parse(text = source_lines[seq_len(export_line)]), envir = globalenv())
