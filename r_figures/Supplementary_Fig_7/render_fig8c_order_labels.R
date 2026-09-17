#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(".", winslash = "/", mustWork = TRUE)
}
source_script <- file.path(script_dir, "make_fig8_eta_order_comparison_final.R")
source_lines <- readLines(source_script, warn = FALSE, encoding = "UTF-8")

replacements <- c(
  '"High-Order(3-nodes)" = "High order",' =
    '"High-Order(3-nodes)" = "High-order",',
  '"Low-Order(2-nodes)" = "Pairwise"' =
    '"Low-Order(2-nodes)" = "Low-order"',
  'base_file <- file.path(out_dir, "Fig8_eta_order_representation_comparison_method_first_label_clear_titlebold_font10p5_label12p5_microadjust")' =
    'base_file <- file.path(out_dir, "Fig8_eta_order_representation_comparison_fig8c_order_labels")'
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
  stop("Could not identify the unique Fig. 8 export line.")
}

write.csv <- function(...) invisible(NULL)
invisible(Sys.setlocale("LC_CTYPE", "Chinese_China.utf8"))
eval(parse(text = source_lines[seq_len(export_line)]), envir = globalenv())
