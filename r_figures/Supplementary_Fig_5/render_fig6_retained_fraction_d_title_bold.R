#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(".", winslash = "/", mustWork = TRUE)
}
source_script <- file.path(script_dir, "make_fig6_eta_robustness_retained_fraction.R")
source_lines <- readLines(source_script, warn = FALSE, encoding = "UTF-8")

source_lines <- sub(
  'title = "Significance tier"',
  'title = "<b>Significance tier</b>"',
  source_lines,
  fixed = TRUE
)
source_lines <- sub(
  'out_base <- file.path(out_dir, "Fig6_eta_robustness_retained_fraction")',
  'out_base <- file.path(out_dir, "Fig6_eta_robustness_retained_fraction_d_title_bold")',
  source_lines,
  fixed = TRUE
)

export_line <- grep("^save_pub_r\\(fig, out_base\\)$", source_lines)
if (length(export_line) != 1L) {
  stop("Could not identify the unique Fig. 6 export line.")
}

eval(parse(text = source_lines[seq_len(export_line)]), envir = globalenv())
